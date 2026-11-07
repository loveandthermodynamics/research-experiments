"""Spatial chemical CTMC with explicit reservoir work and conservative mixing.

The frozen sibling chemistry is unchanged. Added channels obey local detailed
balance with declared reservoir chemical work or with a single thermal bath.
The all-pairs transport control is a supplied kinetic graph, not a stirring
apparatus derived or energetically constructed by the chemistry.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np

PRIOR = Path(__file__).resolve().parent.parent / "interface-maintenance"
sys.path.insert(0, str(PRIOR))
from interface_model import D, L, R0, R1, Event, Flags, Model as BaseModel, Params, State, MATERIAL_SWAP


@dataclass(frozen=True)
class DriveParams:
    exchange_rate: float = 1.0
    fuel_probability: float = 0.9
    mu: float | None = None


@dataclass(frozen=True)
class ExtraEvent:
    kind: str
    i: int
    j: int = -1
    before: int = -1
    after: int = -1


class DrivenModel:
    def __init__(self, params: Params | None = None, flags: Flags | dict[str, bool] | None = None,
                 drive: DriveParams | None = None, transport_mode: str = "local"):
        self.base = BaseModel(params, flags)
        self.params, self.flags, self.n = self.base.params, self.base.flags, self.base.n
        self.drive = drive or DriveParams()
        if transport_mode not in ("local", "mixed_L"):
            raise ValueError("transport_mode must be local or mixed_L")
        self.transport_mode = transport_mode
        d = self.drive
        if not math.isfinite(d.exchange_rate) or d.exchange_rate <= 0:
            raise ValueError("exchange_rate must be finite and positive")
        if d.mu is None:
            if not math.isfinite(d.fuel_probability) or not 0 < d.fuel_probability < 1:
                raise ValueError("fuel_probability must lie strictly between zero and one")
            self.fuel_probability = d.fuel_probability
            self.mu = self.params.E_F + self.params.theta * math.log(d.fuel_probability / (1 - d.fuel_probability))
        else:
            if not math.isfinite(d.mu):
                raise ValueError("mu must be finite")
            self.mu = float(d.mu)
            x = (self.mu - self.params.E_F) / self.params.theta
            self.fuel_probability = 1 / (1 + math.exp(-x)) if x >= 0 else math.exp(x) / (1 + math.exp(x))
            if not 0 < self.fuel_probability < 1:
                raise ValueError("chemical potential exceeds floating-point rate range")
        self.k_up = d.exchange_rate * self.fuel_probability
        self.k_down = d.exchange_rate * (1 - self.fuel_probability)
        self.base_event_count = len(self.base.events)
        reservoirs = tuple(ExtraEvent("reservoir_exchange", i, -1, before, 1 - before)
                           for i in range(self.n) for before in (0, 1))
        self.reservoir_event_count = len(reservoirs)
        self.mix_i, self.mix_j = np.triu_indices(self.n, 1)
        mixing = tuple(ExtraEvent("mixed_L_swap", int(i), int(j)) for i, j in zip(self.mix_i, self.mix_j))
        self.events: tuple[Event | ExtraEvent, ...] = self.base.events + reservoirs + mixing
        self._lookup = {event: k for k, event in enumerate(self.events)}
        reservoir_reverse = self.base_event_count + (np.arange(2 * self.n) ^ 1)
        mix_start = self.base_event_count + self.reservoir_event_count
        self.reverse_index = np.r_[self.base.reverse_index, reservoir_reverse, np.arange(mix_start, len(self.events))]
        self._reservoir_i = np.repeat(np.arange(self.n), 2)
        self._reservoir_before = np.tile([0, 1], self.n)
        self._base_swap = np.flatnonzero(np.array([e.kind for e in self.base.events]) == MATERIAL_SWAP)
        self._mix_adjacent = np.array([int(j) in self.base.neighbors[int(i)] for i, j in zip(self.mix_i, self.mix_j)])
        unions = [np.unique(np.r_[i, j, self.base.neighbors[i], self.base.neighbors[j]])
                  for i, j in zip(self.mix_i, self.mix_j)]
        width = max(map(len, unions))
        self._union_i = np.zeros((len(unions), width), dtype=int)
        self._union_mask = np.zeros_like(self._union_i, dtype=bool)
        for k, union in enumerate(unions):
            self._union_i[k, :len(union)] = union
            self._union_mask[k, :len(union)] = True
        self.mix_mobility = 4 * self.params.material_mobility / (self.n - 1)

    def initial(self, seed: int = 0) -> State:
        return self.base.initial(seed)

    copy_state = staticmethod(BaseModel.copy_state)
    state_to_dict = staticmethod(BaseModel.state_to_dict)
    state_from_dict = staticmethod(BaseModel.state_from_dict)

    def energy(self, state: State) -> float:
        return self.base.energy(state)

    def rates(self, state: State) -> tuple[tuple[Event | ExtraEvent, ...], np.ndarray]:
        _, original = self.base.rates(state)
        reservoir = (state.fuel[self._reservoir_i] == self._reservoir_before) * np.where(self._reservoir_before == 0, self.k_up, self.k_down)
        mixed = np.zeros(len(self.mix_i))
        if self.transport_mode == "mixed_L":
            m = state.material
            bi, bj = self.base._i[self._base_swap], self.base._j[self._base_swap]
            original[self._base_swap[(m[bi] == L) | (m[bj] == L)]] = 0.0
            a, b = m[self.mix_i], m[self.mix_j]
            valid = ((a == L) | (b == L)) & (a != b)
            potential, _ = self.base._potentials(m)
            pair = self.base.pair_energy
            correction = pair[a, a] + pair[b, b] - 2 * pair[a, b]
            delta = (potential[self.mix_i, b] - potential[self.mix_i, a]
                     + potential[self.mix_j, a] - potential[self.mix_j, b]
                     - self._mix_adjacent * correction)
            mobility = np.full(len(a), self.mix_mobility)
            if self.flags.protect_R and self.params.transport_protection_strength:
                nl = ((m[self._union_i] == L) & self._union_mask).sum(axis=1)
                has_r = np.isin(a, (R0, R1)) | np.isin(b, (R0, R1))
                mobility *= np.where(has_r, np.exp(-self.params.transport_protection_strength * nl), 1)
            mixed = valid * mobility * np.exp(-np.maximum(delta, 0) / self.params.theta)
        return self.events, np.r_[original, reservoir, mixed]

    def event_index(self, event: Event | ExtraEvent | int) -> int:
        return int(event) if isinstance(event, (int, np.integer)) else self._lookup[event]

    def reverse(self, event: Event | ExtraEvent | int) -> Event | ExtraEvent:
        return self.events[int(self.reverse_index[self.event_index(event)])]

    def delta_energy(self, state: State, event: Event | ExtraEvent | int) -> float:
        k = self.event_index(event)
        if k < self.base_event_count:
            return self.base.delta_energy(state, k)
        e = self.events[k]
        assert isinstance(e, ExtraEvent)
        if e.kind == "reservoir_exchange":
            return self.params.E_F * (e.after - e.before)
        m, pair, onsite = state.material, self.base.pair_energy, self.base.onsite
        a, b = int(m[e.i]), int(m[e.j])
        delta = (onsite[b] - onsite[a] + (pair[b, m[self.base.neighbors[e.i]]] - pair[a, m[self.base.neighbors[e.i]]]).sum()
                 + onsite[a] - onsite[b] + (pair[a, m[self.base.neighbors[e.j]]] - pair[b, m[self.base.neighbors[e.j]]]).sum())
        if e.j in self.base.neighbors[e.i]:
            delta -= pair[a, a] + pair[b, b] - 2 * pair[a, b]
        return float(delta)

    def apply(self, state: State, event: Event | ExtraEvent | int) -> dict[str, Any]:
        k = self.event_index(event)
        if k < self.base_event_count:
            record = self.base.apply(state, k)
            record.update(chemical_work=0.0, subsystem="internal")
            return record
        e = self.events[k]
        assert isinstance(e, ExtraEvent)
        old_fuel = int(state.fuel.sum())
        delta = self.delta_energy(state, k)
        work = 0.0
        if e.kind == "reservoir_exchange":
            if int(state.fuel[e.i]) != e.before:
                raise ValueError("reservoir edge source does not match fuel state")
            state.fuel[e.i] = e.after
            work = self.mu * (e.after - e.before)
        else:
            a, b = int(state.material[e.i]), int(state.material[e.j])
            if a == b or (a != L and b != L):
                raise ValueError("mixed transport requires unequal species including L")
            state.material[e.i], state.material[e.j] = state.material[e.j], state.material[e.i]
            state.instance_tags[e.i], state.instance_tags[e.j] = state.instance_tags[e.j], state.instance_tags[e.i]
        heat = work - delta
        state.cumulative_heat += heat
        state.external_work += work
        state.event_count += 1
        fuel_after = int(state.fuel.sum())
        fuel_delta = self.params.E_F * (fuel_after - old_fuel)
        return {"time": float(state.time), "event_index": k, "kind": e.kind, "i": e.i, "j": e.j,
                "before": e.before, "after": e.after, "delta_H": delta, "heat_to_bath": heat,
                "external_work": work, "chemical_work": work,
                "fuel_before": old_fuel, "fuel_after": fuel_after,
                "fuel_standard_energy_change": fuel_delta, "material_energy_change": delta - fuel_delta,
                "instance_created": 0, "instance_removed": 0, "parent_instance": 0, "catalyst_instance": 0,
                "subsystem": "reservoir" if e.kind == "reservoir_exchange" else "transport"}

    def summary(self, state: State, include_state: bool = False) -> dict[str, Any]:
        out = self.base.summary(state)
        out["chemical_work"] = state.external_work
        if include_state:
            out["material"] = state.material.tolist()
            out["instance_tags"] = state.instance_tags.tolist()
        return out


def simulate_driven(params: Params | None, seed: int, horizon: float,
                    flags: Flags | dict[str, bool] | None = None, drive: DriveParams | None = None,
                    transport_mode: str = "local", *, sample_dt: float = 2.0,
                    initial_state: State | None = None, max_events: int = 5_000_000) -> dict[str, Any]:
    if horizon < 0 or sample_dt <= 0:
        raise ValueError("invalid observation times")
    model = DrivenModel(params, flags, drive, transport_mode)
    state = model.initial(seed) if initial_state is None else model.copy_state(initial_state)
    if state.time != 0:
        raise ValueError("initial time must be zero")
    initial = model.state_to_dict(state)
    rng = np.random.default_rng(np.random.SeedSequence([seed, 54719]))
    grid = np.r_[np.arange(0.0, horizon, sample_dt), horizon]
    index = 0
    ledger, samples = [], []
    while state.time < horizon:
        _, rates = model.rates(state)
        total = float(rates.sum())
        jump_time = math.inf if total == 0 else state.time + float(rng.exponential(1 / total))
        while index < len(grid) and grid[index] <= min(jump_time, horizon):
            sample = model.summary(state, include_state=True)
            sample["time"] = float(grid[index])
            samples.append(sample)
            index += 1
        if jump_time > horizon:
            state.time = float(horizon)
            break
        state.time = jump_time
        k = int(np.searchsorted(np.cumsum(rates), rng.random() * total, side="right"))
        ledger.append(model.apply(state, k))
        if len(ledger) >= max_events:
            raise RuntimeError("event safeguard exceeded; no partial result returned")
    if horizon == 0:
        samples.append(model.summary(state, include_state=True))
    return {"params": asdict(model.params), "flags": asdict(model.flags), "drive": asdict(model.drive),
            "effective_mu": model.mu, "effective_reservoir_fuel_probability": model.fuel_probability,
            "transport_mode": transport_mode, "seed": seed, "horizon": horizon,
            "initial_state": initial, "final_state": model.state_to_dict(state), "samples": samples,
            "ledger": ledger, "final_summary": model.summary(state),
            "scope": "Supplied uniform chemical reservoirs and transport graph. Local contacts and molecular lineage readouts do not establish a membrane, aging, functional memory, autonomous apparatus, or individuality. Bath heat is not total irreversible entropy production."}
