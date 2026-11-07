"""Finite-fuel local CTMC for production and renewal of protective contacts.

Only physical site species and fuel states enter transition rates. Instance tags,
connected components, ancestry, and any measure of integrity are readouts. This
is a supplied chemistry, not an atomistic derivation or a membrane theorem.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Sequence
import math

import numpy as np

D, R0, R1, L = 0, 1, 2, 3
PASSIVE, COPY, PRODUCE_L, MATERIAL_SWAP, FUEL_SWAP = range(5)
KIND_NAMES = ("passive", "copy_R", "produce_L", "material_swap", "fuel_swap")
SPECIES_NAMES = ("D", "R0", "R1", "L")


@dataclass(frozen=True)
class Params:
    side: int = 10
    theta: float = 1.0
    epsilon_R: float = 8.0
    epsilon_L: float = 5.0
    E_F: float = 12.0
    J_LL: float = 1.5
    J_RL: float = 1.0
    J_RR: float = 0.0
    copy_rate: float = 0.03
    production_rate: float = 0.04
    passive_R: float = 0.08
    passive_L: float = 0.05
    material_mobility: float = 0.10
    fuel_mobility: float = 0.20
    protection_strength: float = 0.50
    transport_protection_strength: float = 0.25
    mutation: tuple[tuple[float, float], tuple[float, float]] = ((0.95, 0.05), (0.05, 0.95))
    allele_passive_multipliers: tuple[float, float] = (1.0, 1.0)
    initial_templates: int = 12
    initial_R1_fraction: float = 0.0
    initial_fuel: int | None = None


@dataclass(frozen=True)
class Flags:
    produce_L: bool = True
    protect_R: bool = True
    copy_R: bool = True


@dataclass(frozen=True)
class Event:
    kind: int
    i: int
    j: int = -1
    before: int = -1
    after: int = -1


@dataclass
class State:
    material: np.ndarray
    fuel: np.ndarray
    instance_tags: np.ndarray
    next_tag: int
    time: float = 0.0
    cumulative_heat: float = 0.0
    external_work: float = 0.0
    event_count: int = 0


class Model:
    """Periodic square-lattice jump process with a fixed reversible catalogue.

    ``rates(state)`` returns ``(model.events, rate_array)``. Zero-rate events
    stay in the catalogue, so event indices remain stable across trajectories
    and factorial flags. ``reverse_index[k]`` is the inverse event index.
    """

    def __init__(self, params: Params | None = None, flags: Flags | dict[str, bool] | None = None):
        self.params = params or Params()
        self.flags = flags if isinstance(flags, Flags) else Flags(**(flags or {}))
        p = self.params
        self._validate_params()
        self.n = p.side * p.side
        site = np.arange(self.n).reshape(p.side, p.side)
        self.neighbors = np.column_stack([
            np.roll(site, 1, axis=0).ravel(), np.roll(site, -1, axis=0).ravel(),
            np.roll(site, 1, axis=1).ravel(), np.roll(site, -1, axis=1).ravel(),
        ])
        self.edge_i = np.repeat(np.arange(self.n), 2)
        self.edge_j = np.column_stack([
            np.roll(site, -1, axis=0).ravel(), np.roll(site, -1, axis=1).ravel()
        ]).ravel()
        self.edge_unions = tuple(np.unique(np.r_[i, j, self.neighbors[i], self.neighbors[j]])
                                 for i, j in zip(self.edge_i, self.edge_j))
        width = max(map(len, self.edge_unions))
        self._union_index = np.zeros((len(self.edge_unions), width), dtype=np.int64)
        self._union_mask = np.zeros_like(self._union_index, dtype=bool)
        for k, u in enumerate(self.edge_unions):
            self._union_index[k, :len(u)] = u
            self._union_mask[k, :len(u)] = True
        self.onsite = np.array([0.0, p.epsilon_R, p.epsilon_R, p.epsilon_L])
        self.pair_energy = np.zeros((4, 4))
        self.pair_energy[L, L] = -p.J_LL
        self.pair_energy[1:3, L] = -p.J_RL
        self.pair_energy[L, 1:3] = -p.J_RL
        self.pair_energy[1:3, 1:3] = -p.J_RR
        self.mutation = np.asarray(p.mutation, dtype=float)
        events: list[Event] = []
        for i in range(self.n):
            for a in (R0, R1, L):
                events.extend((Event(PASSIVE, i, -1, D, a), Event(PASSIVE, i, -1, a, D)))
            for j in self.neighbors[i]:
                for a in (R0, R1):
                    events.extend((Event(COPY, i, int(j), D, a), Event(COPY, i, int(j), a, D)))
                events.extend((Event(PRODUCE_L, i, int(j), D, L),
                               Event(PRODUCE_L, i, int(j), L, D)))
        for i, j in zip(self.edge_i, self.edge_j):
            events.extend((Event(MATERIAL_SWAP, int(i), int(j)), Event(FUEL_SWAP, int(i), int(j))))
        self.events = tuple(events)
        self._event_lookup = {event: k for k, event in enumerate(self.events)}
        self.reverse_index = np.array([
            self._event_lookup[Event(e.kind, e.i, e.j, e.after, e.before)]
            if e.kind <= PRODUCE_L else k for k, e in enumerate(self.events)
        ], dtype=np.int64)
        self._kind = np.array([e.kind for e in self.events])
        self._i = np.array([e.i for e in self.events])
        self._j = np.array([e.j for e in self.events])
        self._before = np.array([e.before for e in self.events])
        self._after = np.array([e.after for e in self.events])
        self._by_kind = tuple(np.flatnonzero(self._kind == kind) for kind in range(5))

    def _validate_params(self) -> None:
        p = self.params
        if p.side < 3 or p.theta <= 0:
            raise ValueError("side must be at least 3 and theta positive")
        if not 0 <= p.initial_templates <= p.side ** 2 or not 0 <= p.initial_R1_fraction <= 1:
            raise ValueError("invalid initial template inventory")
        if p.initial_fuel is not None and not 0 <= p.initial_fuel <= p.side ** 2:
            raise ValueError("invalid finite fuel inventory")
        for name in ("copy_rate", "production_rate", "passive_R", "passive_L", "material_mobility",
                     "fuel_mobility", "protection_strength", "transport_protection_strength"):
            if not math.isfinite(getattr(p, name)) or getattr(p, name) < 0:
                raise ValueError(f"{name} must be finite and nonnegative")
        mu = np.asarray(p.mutation, dtype=float)
        if mu.shape != (2, 2) or not np.all(np.isfinite(mu)) or np.any(mu < 0) or not np.allclose(mu.sum(axis=1), 1):
            raise ValueError("mutation must be a nonnegative 2x2 row-stochastic kernel")
        if len(p.allele_passive_multipliers) != 2 or any(not math.isfinite(x) or x < 0 for x in p.allele_passive_multipliers):
            raise ValueError("allele passive multipliers must be nonnegative")
        for name in ("theta", "epsilon_R", "epsilon_L", "E_F", "J_LL", "J_RL", "J_RR"):
            if not math.isfinite(getattr(p, name)):
                raise ValueError(f"{name} must be finite")

    def initial(self, seed: int = 0) -> State:
        """Dispersed templates, zero L; material/fuel choices do not use flags."""
        rng = np.random.default_rng(seed)
        material = np.zeros(self.n, dtype=np.uint8)
        founders = rng.choice(self.n, self.params.initial_templates, replace=False)
        material[founders] = R0
        n1 = int(round(self.params.initial_templates * self.params.initial_R1_fraction))
        material[founders[:n1]] = R1
        fuel = np.zeros(self.n, dtype=np.uint8)
        k = self.n if self.params.initial_fuel is None else self.params.initial_fuel
        fuel[rng.choice(self.n, k, replace=False)] = 1
        tags = np.zeros(self.n, dtype=np.int64)
        tags[founders] = np.arange(1, len(founders) + 1)
        return State(material, fuel, tags, len(founders) + 1)

    @staticmethod
    def copy_state(state: State) -> State:
        return State(state.material.copy(), state.fuel.copy(), state.instance_tags.copy(),
                     state.next_tag, state.time, state.cumulative_heat, state.external_work, state.event_count)

    def energy(self, state: State) -> float:
        m = state.material
        return float(self.onsite[m].sum() + self.params.E_F * state.fuel.sum()
                     + self.pair_energy[m[self.edge_i], m[self.edge_j]].sum())

    def _potentials(self, material: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        adjacent = material[self.neighbors]
        nr = ((adjacent == R0) | (adjacent == R1)).sum(axis=1)
        nl = (adjacent == L).sum(axis=1)
        p = self.params
        potential = np.zeros((self.n, 4))
        potential[:, R0] = potential[:, R1] = p.epsilon_R - p.J_RR * nr - p.J_RL * nl
        potential[:, L] = p.epsilon_L - p.J_RL * nr - p.J_LL * nl
        return potential, nl

    @staticmethod
    def _environment(environment: Sequence[float] | None) -> np.ndarray:
        a = np.asarray((1.0, 1.0) if environment is None else environment, dtype=float)
        if a.shape != (2,) or not np.all(np.isfinite(a)) or np.any(a < 0):
            raise ValueError("environment must contain two finite nonnegative barrier multipliers")
        return a

    def rates(self, state: State, environment: Sequence[float] | None = None) -> tuple[tuple[Event, ...], np.ndarray]:
        """All rates; environment changes paired R passive mobilities only.

        A constant environment needs no energy protocol. A prescribed changing
        environment is not a modeled physical apparatus and must be declared.
        """
        p, m, f = self.params, state.material, state.fuel
        env = self._environment(environment)
        potential, nl = self._potentials(m)
        out = np.zeros(len(self.events))
        for kind in (PASSIVE, COPY, PRODUCE_L):
            idx = self._by_kind[kind]
            i, j, before, after = self._i[idx], self._j[idx], self._before[idx], self._after[idx]
            valid = m[i] == before
            product = np.maximum(before, after)
            delta = potential[i, after] - potential[i, before]
            if kind == PASSIVE:
                is_r = product != L
                allele = np.clip(product - 1, 0, 1)
                mobility = np.where(is_r, p.passive_R * np.asarray(p.allele_passive_multipliers)[allele] * env[allele], p.passive_L)
                if self.flags.protect_R:
                    mobility = mobility * np.where(is_r, np.exp(-p.protection_strength * nl[i]), 1.0)
            else:
                forward = before == D
                valid &= ((m[j] == R0) | (m[j] == R1)) & (f[i] == forward)
                delta = delta + np.where(forward, -p.E_F, p.E_F)
                if kind == COPY:
                    parent = np.clip(m[j].astype(int) - 1, 0, 1)
                    mobility = p.copy_rate * self.mutation[parent, product - 1] * self.flags.copy_R
                else:
                    mobility = p.production_rate * self.flags.produce_L
            out[idx] = valid * mobility * np.exp(-np.maximum(delta, 0) / p.theta)
        idx = self._by_kind[MATERIAL_SWAP]
        i, j = self._i[idx], self._j[idx]
        a, b = m[i], m[j]
        delta = (potential[i, b] - potential[i, a] + potential[j, a] - potential[j, b]
                 - self.pair_energy[a, a] - self.pair_energy[b, b] + 2 * self.pair_energy[a, b])
        mobility = np.full(len(idx), p.material_mobility)
        if self.flags.protect_R and p.transport_protection_strength:
            union_l = ((m[self._union_index] == L) & self._union_mask).sum(axis=1)
            has_r = (a == R0) | (a == R1) | (b == R0) | (b == R1)
            mobility *= np.where(has_r, np.exp(-p.transport_protection_strength * union_l), 1.0)
        out[idx] = (a != b) * mobility * np.exp(-np.maximum(delta, 0) / p.theta)
        idx = self._by_kind[FUEL_SWAP]
        out[idx] = (f[self._i[idx]] != f[self._j[idx]]) * p.fuel_mobility
        return self.events, out

    def event_index(self, event: Event | int) -> int:
        return int(event) if isinstance(event, (int, np.integer)) else self._event_lookup[event]

    def reverse(self, event: Event | int) -> Event:
        return self.events[int(self.reverse_index[self.event_index(event)])]

    def delta_energy(self, state: State, event: Event | int) -> float:
        e = self.events[self.event_index(event)]
        m, p = state.material, self.params
        if e.kind == FUEL_SWAP:
            return 0.0
        if e.kind == MATERIAL_SWAP:
            a, b = int(m[e.i]), int(m[e.j])
            di = self.onsite[b] - self.onsite[a] + (self.pair_energy[b, m[self.neighbors[e.i]]] - self.pair_energy[a, m[self.neighbors[e.i]]]).sum()
            dj = self.onsite[a] - self.onsite[b] + (self.pair_energy[a, m[self.neighbors[e.j]]] - self.pair_energy[b, m[self.neighbors[e.j]]]).sum()
            return float(di + dj - self.pair_energy[a, a] - self.pair_energy[b, b] + 2 * self.pair_energy[a, b])
        delta = (self.onsite[e.after] - self.onsite[e.before]
                 + (self.pair_energy[e.after, m[self.neighbors[e.i]]] - self.pair_energy[e.before, m[self.neighbors[e.i]]]).sum())
        if e.kind in (COPY, PRODUCE_L):
            delta += -p.E_F if e.before == D else p.E_F
        return float(delta)

    def _check_event(self, state: State, e: Event) -> None:
        if e.kind <= PRODUCE_L:
            if state.material[e.i] != e.before:
                raise ValueError("event source species does not match state")
            if e.kind in (COPY, PRODUCE_L):
                if state.material[e.j] not in (R0, R1) or bool(state.fuel[e.i]) != (e.before == D):
                    raise ValueError("catalytic reaction lacks helper or fuel/waste")
        elif e.kind == MATERIAL_SWAP and state.material[e.i] == state.material[e.j]:
            raise ValueError("identity material swap is not an event")
        elif e.kind == FUEL_SWAP and state.fuel[e.i] == state.fuel[e.j]:
            raise ValueError("identity fuel swap is not an event")

    def apply(self, state: State, event: Event | int) -> dict[str, Any]:
        """Mutate state and return one complete thermodynamic event record.

        This primitive can also apply an energetically allowed catalogue edge
        for audits, even when a factorial flag sets its mobility to zero. The
        SSA selects only positive-rate edges from ``rates``.
        """
        index = self.event_index(event)
        e = self.events[index]
        self._check_event(state, e)
        delta = self.delta_energy(state, e)
        fuel_before = int(state.fuel.sum())
        created = removed = parent = catalyst = 0
        if e.kind <= PRODUCE_L:
            if e.before != D:
                removed = int(state.instance_tags[e.i])
                state.instance_tags[e.i] = 0
            if e.after != D:
                created = state.next_tag
                state.next_tag += 1
                state.instance_tags[e.i] = created
                if e.kind in (COPY, PRODUCE_L):
                    catalyst = int(state.instance_tags[e.j])
                    if e.kind == COPY:
                        parent = catalyst
            state.material[e.i] = e.after
            if e.kind in (COPY, PRODUCE_L):
                state.fuel[e.i] = 0 if e.before == D else 1
        elif e.kind == MATERIAL_SWAP:
            state.material[e.i], state.material[e.j] = state.material[e.j], state.material[e.i]
            state.instance_tags[e.i], state.instance_tags[e.j] = state.instance_tags[e.j], state.instance_tags[e.i]
        else:
            state.fuel[e.i], state.fuel[e.j] = state.fuel[e.j], state.fuel[e.i]
        heat = -delta
        state.cumulative_heat += heat
        state.event_count += 1
        fuel_after = int(state.fuel.sum())
        fuel_energy_change = self.params.E_F * (fuel_after - fuel_before)
        return {"time": float(state.time), "event_index": index, "kind": KIND_NAMES[e.kind],
                "i": e.i, "j": e.j, "before": e.before, "after": e.after,
                "delta_H": delta, "heat_to_bath": heat, "external_work": 0.0,
                "fuel_before": fuel_before, "fuel_after": fuel_after,
                "fuel_standard_energy_change": fuel_energy_change,
                "material_energy_change": delta - fuel_energy_change,
                "instance_created": created, "instance_removed": removed,
                "parent_instance": parent, "catalyst_instance": catalyst}

    def summary(self, state: State) -> dict[str, Any]:
        m, a, b = state.material, state.material[self.edge_i], state.material[self.edge_j]
        ar, br = (a == R0) | (a == R1), (b == R0) | (b == R1)
        return {"time": float(state.time), "D": int((m == D).sum()), "R0": int((m == R0).sum()),
                "R1": int((m == R1).sum()), "L": int((m == L).sum()), "fuel": int(state.fuel.sum()),
                "H": self.energy(state), "heat_to_bath": float(state.cumulative_heat),
                "external_work": float(state.external_work), "RL_contacts": int(((ar & (b == L)) | (br & (a == L))).sum()),
                "LL_contacts": int(((a == L) & (b == L)).sum()), "events": state.event_count}

    @staticmethod
    def state_to_dict(state: State) -> dict[str, Any]:
        return {"material": state.material.tolist(), "fuel": state.fuel.tolist(),
                "instance_tags": state.instance_tags.tolist(), "next_tag": state.next_tag,
                "time": state.time, "cumulative_heat": state.cumulative_heat,
                "external_work": state.external_work, "event_count": state.event_count}

    @staticmethod
    def state_from_dict(data: dict[str, Any]) -> State:
        return State(np.asarray(data["material"], dtype=np.uint8), np.asarray(data["fuel"], dtype=np.uint8),
                     np.asarray(data["instance_tags"], dtype=np.int64), int(data["next_tag"]),
                     float(data.get("time", 0)), float(data.get("cumulative_heat", 0)),
                     float(data.get("external_work", 0)), int(data.get("event_count", 0)))


def simulate(params: Params | None, seed: int, horizon: float,
             flags: Flags | dict[str, bool] | None = None, *, sample_dt: float = 1.0,
             environment: Sequence[float] | None = None, initial_state: State | None = None,
             max_events: int = 2_000_000) -> dict[str, Any]:
    """Exact direct SSA; no snapshots, components, or thresholds affect rates."""
    if horizon < 0 or sample_dt <= 0:
        raise ValueError("horizon must be nonnegative and sample_dt positive")
    model = Model(params, flags)
    state = model.initial(seed) if initial_state is None else model.copy_state(initial_state)
    if state.time != 0:
        raise ValueError("initial state for simulate must start at time zero")
    initial = model.state_to_dict(state)
    rng = np.random.default_rng(np.random.SeedSequence([seed, 914173]))
    samples: list[dict[str, Any]] = []
    ledger: list[dict[str, Any]] = []
    sample_times = np.r_[np.arange(0.0, horizon, sample_dt), horizon]
    sample_index = 0
    while state.time < horizon:
        _, rates = model.rates(state, environment)
        total = float(rates.sum())
        next_time = math.inf if total == 0 else state.time + float(rng.exponential(1 / total))
        while sample_index < len(sample_times) and sample_times[sample_index] <= min(next_time, horizon):
            sample = model.summary(state)
            sample["time"] = float(sample_times[sample_index])
            samples.append(sample)
            sample_index += 1
        if next_time > horizon:
            state.time = float(horizon)
            break
        state.time = next_time
        event_index = int(np.searchsorted(np.cumsum(rates), rng.random() * total, side="right"))
        ledger.append(model.apply(state, event_index))
        if len(ledger) >= max_events:
            raise RuntimeError("max_events exceeded; no truncated result is returned")
    if horizon == 0:
        samples.append(model.summary(state))
    return {"params": asdict(model.params), "flags": asdict(model.flags), "seed": seed,
            "environment": model._environment(environment).tolist(), "horizon": horizon,
            "initial_state": initial, "final_state": model.state_to_dict(state),
            "samples": samples, "ledger": ledger,
            "final_summary": model.summary(state),
            "scope": "Physical protective contacts; no enclosing membrane, aging law, functional memory, or group individuality is inferred. Instance tags track molecular instances, not conserved atoms."}


if __name__ == "__main__":
    import json
    result = simulate(Params(side=6, initial_templates=6), seed=7, horizon=5.0)
    print(json.dumps({"params": result["params"], "flags": result["flags"],
                      "initial_H": Model(Params(**result["params"])).energy(Model.state_from_dict(result["initial_state"])),
                      "final": result["final_summary"], "events": len(result["ledger"])}, indent=2))
