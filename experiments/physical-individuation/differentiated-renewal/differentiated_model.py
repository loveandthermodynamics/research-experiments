"""Supplied opposed/aligned response phenotypes in the frozen driven chemistry.

A stationary, independent symmetric telegraph gate changes only paired passive
template mobilities. Site energies, copying, reservoirs and neutral instance
tags retain their prior meanings. No group or diversity measurement enters rates.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from pathlib import Path
import sys
from typing import Any, Sequence

import numpy as np

DRIVEN = Path(__file__).resolve().parent.parent / "driven-spatial"
sys.path.insert(0, str(DRIVEN))
from driven_model import DriveParams, DrivenModel, Flags, Params, State
from interface_model import COPY, PRODUCE_L, PASSIVE, Event


@dataclass(frozen=True)
class ResponseParams:
    asymmetry: float = 0.5
    flip_rate: float = 0.025
    initial_gate: int | None = None


@dataclass(frozen=True)
class GateEvent:
    before: int
    after: int
    kind: str = "gate_flip"


@dataclass
class JointState:
    base_state: State
    gate: int
    gate_event_count: int = 0

    @property
    def time(self) -> float:
        return self.base_state.time

    @time.setter
    def time(self, value: float) -> None:
        self.base_state.time = float(value)


class DifferentiatedModel:
    def __init__(self, params: Params | None = None, flags: Flags | dict[str, bool] | None = None,
                 drive: DriveParams | None = None, response: ResponseParams | None = None,
                 response_mode: str = "opposed", transport_mode: str = "local",
                 production_site_mobility: Sequence[float] | None = None):
        self.driven = DrivenModel(params, flags, drive, transport_mode)
        self.params, self.flags, self.drive = self.driven.params, self.driven.flags, self.driven.drive
        self.n, self.transport_mode = self.driven.n, transport_mode
        self.response = response or ResponseParams()
        if response_mode not in ("opposed", "aligned"):
            raise ValueError("response_mode must be opposed or aligned")
        self.response_mode = response_mode
        r = self.response
        if not math.isfinite(r.asymmetry) or not 0 <= r.asymmetry < 1:
            raise ValueError("asymmetry must be finite and lie in [0,1)")
        if not math.isfinite(r.flip_rate) or r.flip_rate <= 0:
            raise ValueError("flip_rate must be finite and positive")
        if r.initial_gate is not None and r.initial_gate not in (-1, 1):
            raise ValueError("physical gate is -1/+1; use None for stationary initialization")
        mob = np.ones(self.n) if production_site_mobility is None else np.array(production_site_mobility, dtype=float, copy=True)
        if mob.shape != (self.n,) or not np.all(np.isfinite(mob)) or np.any(mob < 0):
            raise ValueError("production_site_mobility must have N finite nonnegative entries")
        self._production_site_mobility = tuple(float(x) for x in mob)
        self._mobility = np.asarray(self._production_site_mobility)
        self._mobility.flags.writeable = False
        self.driven_event_count = len(self.driven.events)
        self.events = self.driven.events + (GateEvent(-1, 1), GateEvent(1, -1))
        self.reverse_index = np.r_[self.driven.reverse_index, self.driven_event_count + np.array([1, 0])]
        self._lookup = {event: index for index, event in enumerate(self.events)}
        passive, alleles, production, helpers = [], [], [], []
        for k, event in enumerate(self.driven.events):
            if isinstance(event, Event) and event.kind == PASSIVE and max(event.before, event.after) in (1, 2):
                passive.append(k); alleles.append(max(event.before, event.after) - 1)
            if isinstance(event, Event) and event.kind == PRODUCE_L:
                production.append(k); helpers.append(event.j)
        self._passive_indices, self._alleles = np.array(passive), np.array(alleles)
        self._production_indices, self._producer_sites = np.array(production), np.array(helpers)

    @property
    def production_site_mobility(self) -> tuple[float, ...]:
        return self._production_site_mobility

    def initial(self, seed: int = 0) -> JointState:
        gate = self.response.initial_gate
        if gate is None:
            rng = np.random.default_rng(np.random.SeedSequence([seed, 62411]))
            gate = int(rng.choice((-1, 1)))
        return JointState(self.driven.initial(seed), gate)

    @staticmethod
    def copy_state(state: JointState) -> JointState:
        return JointState(DrivenModel.copy_state(state.base_state), state.gate, state.gate_event_count)

    @staticmethod
    def state_to_dict(state: JointState) -> dict[str, Any]:
        return {"base_state": DrivenModel.state_to_dict(state.base_state), "gate": state.gate,
                "gate_event_count": state.gate_event_count}

    @staticmethod
    def state_from_dict(data: dict[str, Any]) -> JointState:
        gate = int(data["gate"])
        if gate not in (-1, 1):
            raise ValueError("saved physical gate must be -1/+1")
        return JointState(DrivenModel.state_from_dict(data["base_state"]), gate, int(data.get("gate_event_count", 0)))

    def energy(self, state: JointState) -> float:
        return self.driven.energy(state.base_state)  # Both physical gate wells have energy zero.

    def passive_multipliers(self, state: JointState) -> tuple[float, float]:
        first = 1 + self.response.asymmetry * state.gate
        second = first if self.response_mode == "aligned" else 1 - self.response.asymmetry * state.gate
        return first, second

    def material_rates(self, state: JointState):
        events, rates = self.driven.rates(state.base_state)
        rates[self._passive_indices] *= np.asarray(self.passive_multipliers(state))[self._alleles]
        rates[self._production_indices] *= self._mobility[self._producer_sites]
        return events, rates

    def rates(self, state: JointState):
        _, rates = self.material_rates(state)
        gate_rates = np.zeros(2)
        gate_rates[0 if state.gate == -1 else 1] = self.response.flip_rate
        return self.events, np.r_[rates, gate_rates]

    def event_index(self, event) -> int:
        return int(event) if isinstance(event, (int, np.integer)) else self._lookup[event]

    def reverse(self, event):
        return self.events[int(self.reverse_index[self.event_index(event)])]

    def delta_energy(self, state: JointState, event) -> float:
        k = self.event_index(event)
        return self.driven.delta_energy(state.base_state, k) if k < self.driven_event_count else 0.0

    def apply(self, state: JointState, event) -> dict[str, Any]:
        k = self.event_index(event)
        before_gate = state.gate
        helper = -1
        if k < self.driven_event_count:
            edge = self.events[k]
            if isinstance(edge, Event) and edge.kind in (COPY, PRODUCE_L):
                helper = int(state.base_state.material[edge.j])
            record = self.driven.apply(state.base_state, k)
        else:
            edge = self.events[k]
            assert isinstance(edge, GateEvent)
            if state.gate != edge.before:
                raise ValueError("gate transition source does not match state")
            state.gate = edge.after
            state.gate_event_count += 1
            fuel = int(state.base_state.fuel.sum())
            record = {"time": state.time, "event_index": k, "kind": "gate_flip", "i": -1, "j": -1,
                      "before": edge.before, "after": edge.after, "delta_H": 0.0, "heat_to_bath": 0.0,
                      "external_work": 0.0, "chemical_work": 0.0, "fuel_before": fuel, "fuel_after": fuel,
                      "fuel_standard_energy_change": 0.0, "material_energy_change": 0.0,
                      "instance_created": 0, "instance_removed": 0, "parent_instance": 0, "catalyst_instance": 0,
                      "subsystem": "environment"}
        record.update(gate_before=before_gate, gate_after=state.gate, helper_species=helper)
        return record

    def summary(self, state: JointState, include_state=False):
        out = self.driven.summary(state.base_state, include_state)
        out["material_reservoir_transport_events"] = out["events"]
        out["events"] += state.gate_event_count
        out.update(gate=state.gate, gate_events=state.gate_event_count,
                   passive_multipliers=list(self.passive_multipliers(state)))
        return out


def simulate_differentiated(params: Params | None, seed: int, horizon: float,
                            flags: Flags | dict[str, bool] | None = None, drive: DriveParams | None = None,
                            response: ResponseParams | None = None, response_mode: str = "opposed",
                            transport_mode: str = "local", *, production_site_mobility: Sequence[float] | None = None,
                            initial_state: JointState | None = None, sample_dt=2.0, max_events=5_000_000,
                            max_gate_events=1_000_000):
    if horizon < 0 or sample_dt <= 0:
        raise ValueError("invalid observation times")
    model = DifferentiatedModel(params, flags, drive, response, response_mode, transport_mode, production_site_mobility)
    state = model.initial(seed) if initial_state is None else model.copy_state(initial_state)
    if state.time != 0:
        raise ValueError("initial time must be zero; preserve realized gate when resetting bookkeeping")
    initial = model.state_to_dict(state)
    gate_rng = np.random.default_rng(np.random.SeedSequence([seed, 62412]))
    material_rng = np.random.default_rng(np.random.SeedSequence([seed, 62413]))
    next_gate = float(gate_rng.exponential(1 / model.response.flip_rate))
    grid = np.r_[np.arange(0.0, horizon, sample_dt), horizon]
    index = 0
    ledger, samples = [], []
    while state.time < horizon:
        _, rates = model.material_rates(state)
        total = float(rates.sum())
        material_time = math.inf if total == 0 else state.time + float(material_rng.exponential(1 / total))
        next_time = min(material_time, next_gate)
        while index < len(grid) and grid[index] <= min(next_time, horizon):
            sample = model.summary(state, include_state=True); sample["time"] = float(grid[index])
            samples.append(sample); index += 1
        if next_time > horizon:
            state.time = float(horizon); break
        state.time = next_time
        if material_time < next_gate:
            k = int(np.searchsorted(np.cumsum(rates), material_rng.random() * total, side="right"))
        else:
            k = model.driven_event_count + (0 if state.gate == -1 else 1)
        ledger.append(model.apply(state, k))
        if k >= model.driven_event_count:
            next_gate = state.time + float(gate_rng.exponential(1 / model.response.flip_rate))
        if len(ledger) >= max_events or state.gate_event_count >= max_gate_events:
            raise RuntimeError("event safeguard exceeded; no truncated result returned")
    if horizon == 0:
        samples.append(model.summary(state, include_state=True))
    return {"params": asdict(model.params), "flags": asdict(model.flags), "drive": asdict(model.drive),
            "effective_mu": model.driven.mu, "response": asdict(model.response), "response_mode": response_mode,
            "transport_mode": transport_mode, "production_site_mobility": list(model.production_site_mobility),
            "seed": seed, "horizon": horizon, "initial_state": initial, "final_state": model.state_to_dict(state),
            "samples": samples, "ledger": ledger, "final_summary": model.summary(state),
            "algorithm": "Independent stationary initial gate, autonomous exponential gate clock, and conditional material/reservoir SSA; matched seeds share gate histories.",
            "scope": "Supplied finite response alphabet, copying grammar, reservoirs and kinetic gate. Immutable producer-site mobility is an external fixed-site kinetic intervention with excluded apparatus construction cost, never genealogy feedback. Tags track instances/copying provenance, not functional identity, aging, subjectivity or conserved atoms."}
