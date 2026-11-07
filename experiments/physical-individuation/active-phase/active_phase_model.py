"""Local reversible phase chemistry with explicit ideal-chemostat work.

Natural rates inspect only local occupancies and one Hamiltonian. No cluster,
interface, lineage or preferred-size observable enters the generator.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Any

import numpy as np


S, A, B = 0, 1, 2
MATERIAL_SWAP, ENZYME_SWAP = "material_swap", "enzyme_swap"
PASSIVE, CATALYZED = "passive", "catalyzed"


@dataclass(frozen=True)
class Params:
    side: int = 16
    theta: float = 1.0
    epsilon_B: float = 2.0
    J_BB: float = 2.4
    J_EB: float = 2.5
    mu: float = 12.0
    material_mobility: float = 1.0
    enzyme_mobility: float = 0.5
    passive_rate: float = 0.02
    catalytic_rate: float = 0.2
    material_count: int = 77
    enzyme_count: int = 8


@dataclass
class State:
    material: np.ndarray
    enzyme: np.ndarray
    material_tags: np.ndarray
    time: float = 0.0
    cumulative_work: float = 0.0
    cumulative_heat: float = 0.0
    net_fuel_to_waste: int = 0
    event_count: int = 0


@dataclass(frozen=True)
class Event:
    kind: str
    i: int
    j: int = -1
    before: int = -1
    after: int = -1


class Model:
    def __init__(self, params: Params | None = None):
        self.params = params or Params()
        p = self.params
        if not isinstance(p.side, int) or p.side < 4:
            raise ValueError("side must be an integer >=4")
        self.n = p.side ** 2
        for name in ("material_count", "enzyme_count"):
            value = getattr(p, name)
            if not isinstance(value, int) or not 0 <= value <= self.n:
                raise ValueError(f"{name} must be an integer in [0,N]")
        for name in ("theta", "epsilon_B", "J_BB", "J_EB", "mu", "material_mobility",
                     "enzyme_mobility", "passive_rate", "catalytic_rate"):
            if not math.isfinite(getattr(p, name)):
                raise ValueError(f"{name} must be finite")
        if p.theta <= 0 or min(p.material_mobility, p.enzyme_mobility,
                              p.passive_rate, p.catalytic_rate) < 0:
            raise ValueError("temperature must be positive and mobilities nonnegative")
        neighbors, bonds = [], []
        for i in range(self.n):
            row, col = divmod(i, p.side)
            right = row * p.side + (col + 1) % p.side
            down = ((row + 1) % p.side) * p.side + col
            left = row * p.side + (col - 1) % p.side
            up = ((row - 1) % p.side) * p.side + col
            neighbors.append((right, down, left, up))
            bonds.extend(((i, right), (i, down)))
        self.neighbors = np.asarray(neighbors, dtype=np.int64)
        self.bonds = np.asarray(bonds, dtype=np.int64)
        self.bond_i, self.bond_j = self.bonds.T
        self.edges_per_family = 2 * self.n
        events = [Event(MATERIAL_SWAP, int(i), int(j)) for i, j in bonds]
        events += [Event(ENZYME_SWAP, int(i), int(j)) for i, j in bonds]
        events += [Event(PASSIVE, i, before=before, after=after)
                   for i in range(self.n) for before, after in ((A, B), (B, A))]
        events += [Event(CATALYZED, i, before=before, after=after)
                   for i in range(self.n) for before, after in ((B, A), (A, B))]
        self.events = tuple(events)
        self._event_lookup = {event: index for index, event in enumerate(events)}
        swaps = np.arange(4 * self.n)
        conversions = np.arange(4 * self.n, 8 * self.n) ^ 1
        self.reverse_index = np.r_[swaps, conversions]
        self.proposal_mobilities = np.repeat(
            np.array([p.material_mobility, p.enzyme_mobility, p.passive_rate, p.catalytic_rate]),
            self.edges_per_family,
        )
        self._proposal_family_cumulative = np.cumsum(
            np.array([p.material_mobility, p.enzyme_mobility, p.passive_rate, p.catalytic_rate])
        )
        self.total_proposal_rate = float(self.edges_per_family * self._proposal_family_cumulative[-1])

    def initial(self, seed: int = 0) -> State:
        material_rng = np.random.default_rng(np.random.SeedSequence([seed, 74501]))
        enzyme_rng = np.random.default_rng(np.random.SeedSequence([seed, 74502]))
        material = np.zeros(self.n, dtype=np.int8)
        tags = np.zeros(self.n, dtype=np.int64)
        sites = material_rng.choice(self.n, self.params.material_count, replace=False)
        material[sites] = A
        tags[sites] = np.arange(1, self.params.material_count + 1)
        enzyme = np.zeros(self.n, dtype=np.int8)
        enzyme[enzyme_rng.choice(self.n, self.params.enzyme_count, replace=False)] = 1
        return State(material, enzyme, tags)

    @staticmethod
    def copy_state(state: State) -> State:
        return State(state.material.copy(), state.enzyme.copy(), state.material_tags.copy(), state.time,
                     state.cumulative_work, state.cumulative_heat, state.net_fuel_to_waste, state.event_count)

    @staticmethod
    def state_to_dict(state: State) -> dict[str, Any]:
        return dict(material=state.material.tolist(), enzyme=state.enzyme.tolist(),
                    material_tags=state.material_tags.tolist(), time=float(state.time),
                    cumulative_work=float(state.cumulative_work), cumulative_heat=float(state.cumulative_heat),
                    net_fuel_to_waste=int(state.net_fuel_to_waste), event_count=int(state.event_count))

    @staticmethod
    def state_from_dict(data: dict[str, Any]) -> State:
        return State(np.asarray(data["material"], dtype=np.int8), np.asarray(data["enzyme"], dtype=np.int8),
                     np.asarray(data["material_tags"], dtype=np.int64), float(data.get("time", 0.)),
                     float(data.get("cumulative_work", 0.)), float(data.get("cumulative_heat", 0.)),
                     int(data.get("net_fuel_to_waste", 0)), int(data.get("event_count", 0)))

    def validate_state(self, state: State) -> None:
        if any(array.shape != (self.n,) for array in (state.material, state.enzyme, state.material_tags)):
            raise ValueError("state arrays must have exactly N entries")
        if not np.isin(state.material, (S, A, B)).all() or not np.isin(state.enzyme, (0, 1)).all():
            raise ValueError("invalid material/enzyme occupancy")
        present = state.material != S
        tags = state.material_tags[present]
        if int(present.sum()) != self.params.material_count or int(state.enzyme.sum()) != self.params.enzyme_count:
            raise ValueError("state inventory differs from supplied inventory")
        if np.any(state.material_tags[~present] != 0) or np.any(tags <= 0) or len(np.unique(tags)) != len(tags):
            raise ValueError("conserved material tags must be unique, positive and absent on solvent")
        if not math.isfinite(state.time) or state.time < 0:
            raise ValueError("state time must be finite and nonnegative")

    def energy(self, state: State) -> float:
        b = state.material == B
        p = self.params
        return float(p.epsilon_B * b.sum() - p.J_BB * np.sum(b[self.bond_i] & b[self.bond_j])
                     - p.J_EB * np.dot(state.enzyme, b.astype(np.int64)))

    def event_index(self, event: Event | int) -> int:
        return int(event) if isinstance(event, (int, np.integer)) else self._event_lookup[event]

    def reverse(self, event: Event | int) -> Event:
        return self.events[int(self.reverse_index[self.event_index(event)])]

    def work(self, event: Event | int) -> float:
        edge = self.events[self.event_index(event)]
        if edge.kind != CATALYZED:
            return 0.
        return self.params.mu if edge.before == B else -self.params.mu

    def delta_energy(self, state: State, event: Event | int) -> float:
        edge = self.events[self.event_index(event)]
        p, i, j = self.params, edge.i, edge.j
        b_i = int(state.material[i] == B)
        if edge.kind == MATERIAL_SWAP:
            b_j = int(state.material[j] == B)
            change = b_j - b_i
            near_i = int(np.count_nonzero(state.material[self.neighbors[i]] == B))
            near_j = int(np.count_nonzero(state.material[self.neighbors[j]] == B))
            # Endpoints are neighbors. Their mutual BB product is unchanged.
            return float(-p.J_BB * change * (near_i - near_j - change)
                         - p.J_EB * change * (int(state.enzyme[i]) - int(state.enzyme[j])))
        if edge.kind == ENZYME_SWAP:
            return float(-p.J_EB * (int(state.enzyme[j]) - int(state.enzyme[i]))
                         * (b_i - int(state.material[j] == B)))
        if int(state.material[i]) != edge.before:
            raise ValueError("conversion source does not match microscopic state")
        change = int(edge.after == B) - b_i
        field = p.epsilon_B - p.J_BB * np.count_nonzero(state.material[self.neighbors[i]] == B)
        field -= p.J_EB * int(state.enzyme[i])
        return float(change * field)

    def event_rate(self, state: State, event: Event | int) -> float:
        k = self.event_index(event)
        edge, mobility = self.events[k], float(self.proposal_mobilities[k])
        if mobility == 0:
            return 0.
        if edge.kind == MATERIAL_SWAP:
            enabled = state.material[edge.i] != state.material[edge.j]
        elif edge.kind == ENZYME_SWAP:
            enabled = state.enzyme[edge.i] != state.enzyme[edge.j]
        else:
            enabled = state.material[edge.i] == edge.before
            if edge.kind == CATALYZED:
                enabled = enabled and bool(state.enzyme[edge.i])
        if not enabled:
            return 0.
        exponent = min(0., (self.work(k) - self.delta_energy(state, k)) / self.params.theta)
        value = mobility * math.exp(exponent)
        if value == 0:
            raise FloatingPointError("enabled reverse rate underflow; parameter scale is not representable")
        return value

    def rates(self, state: State):
        """Vectorized public reference API; the simulator uses local thinning."""
        p = self.params
        b = (state.material == B).astype(np.int64)
        neighbor_b = b[self.neighbors].sum(axis=1)
        i, j = self.bond_i, self.bond_j
        change = b[j] - b[i]
        delta_material = -p.J_BB * change * (neighbor_b[i] - neighbor_b[j] - change)
        delta_material -= p.J_EB * change * (state.enzyme[i].astype(float) - state.enzyme[j])
        delta_enzyme = -p.J_EB * (state.enzyme[j].astype(float) - state.enzyme[i]) * (b[i] - b[j])
        field = p.epsilon_B - p.J_BB * neighbor_b - p.J_EB * state.enzyme
        delta_passive = np.stack((field, -field), axis=1).ravel()
        delta_catalyzed = np.stack((-field, field), axis=1).ravel()
        deltas = np.r_[delta_material, delta_enzyme, delta_passive, delta_catalyzed]
        works = np.r_[np.zeros(6 * self.n), np.tile((p.mu, -p.mu), self.n)]
        passive_enabled = np.stack((state.material == A, state.material == B), axis=1).ravel()
        catalyzed_enabled = np.stack(((state.material == B) & (state.enzyme == 1),
                                      (state.material == A) & (state.enzyme == 1)), axis=1).ravel()
        enabled = np.r_[state.material[i] != state.material[j], state.enzyme[i] != state.enzyme[j],
                        passive_enabled, catalyzed_enabled]
        rates = self.proposal_mobilities * enabled * np.exp(np.minimum(0., (works - deltas) / p.theta))
        if np.any(enabled & (self.proposal_mobilities > 0) & (rates == 0)):
            raise FloatingPointError("enabled rate underflow; parameter scale is not representable")
        return self.events, rates

    @staticmethod
    def _local_snapshot(state: State, sites: list[int]) -> dict[str, list[int]]:
        return dict(material=[int(state.material[i]) for i in sites],
                    enzyme=[int(state.enzyme[i]) for i in sites],
                    material_tags=[int(state.material_tags[i]) for i in sites])

    def apply(self, state: State, event: Event | int) -> dict[str, Any]:
        k = self.event_index(event)
        edge = self.events[k]
        if self.event_rate(state, k) <= 0:
            raise ValueError("cannot apply disabled microscopic event")
        sites = [edge.i] if edge.j == -1 else [edge.i, edge.j]
        before = self._local_snapshot(state, sites)
        delta, work = self.delta_energy(state, k), self.work(k)
        fuel = (1 if edge.before == B else -1) if edge.kind == CATALYZED else 0
        if edge.kind == MATERIAL_SWAP:
            state.material[edge.i], state.material[edge.j] = state.material[edge.j], state.material[edge.i]
            state.material_tags[edge.i], state.material_tags[edge.j] = state.material_tags[edge.j], state.material_tags[edge.i]
        elif edge.kind == ENZYME_SWAP:
            state.enzyme[edge.i], state.enzyme[edge.j] = state.enzyme[edge.j], state.enzyme[edge.i]
        else:
            state.material[edge.i] = edge.after
        heat = work - delta
        state.cumulative_work += work
        state.cumulative_heat += heat
        state.net_fuel_to_waste += fuel
        state.event_count += 1
        return dict(time=float(state.time), event_index=k, kind=edge.kind, i=edge.i, j=edge.j,
                    before=before, after=self._local_snapshot(state, sites), delta_H=delta,
                    chemical_work=work, heat_to_bath=heat, fuel_to_waste=fuel,
                    cumulative_work=float(state.cumulative_work), cumulative_heat=float(state.cumulative_heat),
                    net_fuel_to_waste=int(state.net_fuel_to_waste), event_count=state.event_count)

    def summary(self, state: State, include_state: bool = False) -> dict[str, Any]:
        b = state.material == B
        result = dict(time=float(state.time), H=self.energy(state), S=int((state.material == S).sum()),
                      A=int((state.material == A).sum()), B=int(b.sum()), enzymes=int(state.enzyme.sum()),
                      enzymes_on_B=int(np.dot(state.enzyme, b.astype(np.int64))),
                      BB_contacts=int(np.sum(b[self.bond_i] & b[self.bond_j])),
                      chemical_work=float(state.cumulative_work), heat_to_bath=float(state.cumulative_heat),
                      net_fuel_to_waste=int(state.net_fuel_to_waste), events=state.event_count)
        if include_state:
            result.update(material=state.material.tolist(), enzyme=state.enzyme.tolist(),
                          material_tags=state.material_tags.tolist())
        return result


def simulate(params: Params | None, seed: int, horizon: float, *, initial_state: State | None = None,
             sample_dt: float = 2., max_proposals: int = 10000000, max_events: int = 5000000):
    """Exact constant-bound proposal thinning, including rejected-event time.

    Clock/event/acceptance draws use one documented stream independent of
    initialization. Rejected proposals cause no microscopic change. Accepted
    timestamps, proposal numbers, counts and terminal censoring are retained.
    """
    if not math.isfinite(horizon) or horizon < 0 or not math.isfinite(sample_dt) or sample_dt <= 0:
        raise ValueError("invalid observation times")
    if max_proposals <= 0 or max_events <= 0:
        raise ValueError("safeguards must be positive")
    model = Model(params)
    state = model.initial(seed) if initial_state is None else model.copy_state(initial_state)
    model.validate_state(state)
    if state.time != 0 or state.event_count != 0 or state.cumulative_work != 0 or state.cumulative_heat != 0 or state.net_fuel_to_waste != 0:
        raise ValueError("initial observer counters must be zero; physical inventories are preserved")
    initial = model.state_to_dict(state)
    initial_h = model.energy(state)
    rng = np.random.default_rng(np.random.SeedSequence([seed, 74503]))
    grid = np.r_[np.arange(0., horizon, sample_dt), horizon]
    samples, ledger = [], []
    index, proposals, rejected_since_last, accepted = 0, 0, 0, 0
    proposed_by_family = np.zeros(4, dtype=np.int64)
    accepted_by_family = np.zeros(4, dtype=np.int64)
    final_censored_proposal_time = None
    while state.time < horizon:
        if model.total_proposal_rate == 0:
            next_time = math.inf
        else:
            next_time = state.time + float(rng.exponential(1. / model.total_proposal_rate))
        while index < len(grid) and grid[index] <= min(next_time, horizon):
            sample = model.summary(state, include_state=True)
            sample["time"] = float(grid[index])
            samples.append(sample)
            index += 1
        if next_time > horizon:
            final_censored_proposal_time = next_time if math.isfinite(next_time) else None
            state.time = float(horizon)
            break
        state.time = next_time
        proposals += 1
        family = int(np.searchsorted(model._proposal_family_cumulative,
                                     rng.random() * model._proposal_family_cumulative[-1], side="right"))
        offset = int(rng.integers(model.edges_per_family))
        k = family * model.edges_per_family + offset
        proposed_by_family[family] += 1
        probability = model.event_rate(state, k) / model.proposal_mobilities[k]
        if rng.random() < probability:
            record = model.apply(state, k)
            record.update(proposal_number=proposals, rejected_since_last_accepted=rejected_since_last)
            ledger.append(record)
            accepted += 1
            accepted_by_family[family] += 1
            rejected_since_last = 0
        else:
            rejected_since_last += 1
        if proposals >= max_proposals or accepted >= max_events:
            raise RuntimeError("event safeguard exceeded; no truncated result returned")
    if horizon == 0:
        samples.append(model.summary(state, include_state=True))
    model.validate_state(state)
    first_law_error = abs(model.energy(state) - initial_h - state.cumulative_work + state.cumulative_heat)
    if first_law_error > 1e-6:
        raise RuntimeError("cumulative first law failed")
    families = (MATERIAL_SWAP, ENZYME_SWAP, PASSIVE, CATALYZED)
    return dict(
        params=asdict(model.params), seed=seed, horizon=float(horizon), initial_state=initial,
        final_state=model.state_to_dict(state), ledger=ledger, samples=samples,
        final_summary=model.summary(state), first_law_error=first_law_error,
        proposal_statistics=dict(
            total_bound=model.total_proposal_rate, proposals=proposals, accepted=accepted,
            rejected=proposals - accepted, terminal_rejected_since_last_accepted=rejected_since_last,
            proposed_by_family=dict(zip(families, map(int, proposed_by_family))),
            accepted_by_family=dict(zip(families, map(int, accepted_by_family))),
            terminal_censored_proposal_time=final_censored_proposal_time,
        ),
        algorithm="Exact constant-bound proposal thinning: uniform channel within family; family weighted by fixed mobility; exponential total-bound clock; all rejected waiting times included.",
        random_streams={"initial_material": [seed, 74501], "initial_enzyme": [seed, 74502],
                        "proposal_clock_event_acceptance": [seed, 74503]},
        supplied_resources="Randomized A carriers, conserved mobile enzymes, lattice/affinities and reaction grammar, thermal bath, homogeneous ideal F/W chemostats with fixed mu. Chemical work counts fuel-to-waste conversion, not apparatus construction or initial molecular preparation.",
        scope="No droplet seed, inside/outside mask, group or diversity reward. Formation, selected finite size, arrested coarsening, turnover and individuation are distinct empirical questions. Conserved enzyme is not renewed metabolism; copying and variation are absent.",
    )
