"""Autonomous, finite-state kinetic environment for the frozen local CTMC.

Two physical binary gates have declared well energies and thermal transition
rates. They modulate forward and reverse passive-template barriers together.
The construction is a coarse kinetic-gate idealization, not a derivation of a
weather apparatus. No damage reset or group-dependent rule is added.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any
import math

import numpy as np

from interface_model import Event, Flags, Model as BaseModel, Params, State


@dataclass(frozen=True)
class GateParams:
    p_bad: tuple[float, float] = (0.08, 0.08)
    gamma: tuple[float, float] = (0.25, 0.25)
    bad_multiplier: tuple[float, float] = (12.0, 12.0)
    initial_gates: tuple[int, int] = (0, 0)


@dataclass(frozen=True)
class GateEvent:
    gate: int
    before: int
    after: int
    kind: str = "environment_gate"


@dataclass
class EnvironmentState:
    base_state: State
    gates: np.ndarray
    cumulative_gate_heat: float = 0.0
    gate_event_count: int = 0

    @property
    def time(self) -> float:
        return self.base_state.time

    @time.setter
    def time(self, value: float) -> None:
        self.base_state.time = float(value)


class EnvironmentModel:
    """Joint local-material and two-gate Markov jump generator.

    Base event indices retain their original values; four directed gate edges
    are appended. ``reverse_index`` and ``energy`` refer to the full joint
    physical state. Observer-assigned instance IDs are not thermodynamic state.
    """

    def __init__(self, params: Params | None = None, flags: Flags | dict[str, bool] | None = None,
                 gateparams: GateParams | None = None, response_mode: str = "differentiated"):
        self.base = BaseModel(params, flags)
        self.params, self.flags = self.base.params, self.base.flags
        self.gateparams = gateparams or GateParams()
        if response_mode not in ("differentiated", "common"):
            raise ValueError("response_mode must be 'differentiated' or 'common'")
        self.response_mode = response_mode
        gp = self.gateparams
        for name in ("p_bad", "gamma", "bad_multiplier", "initial_gates"):
            if len(getattr(gp, name)) != 2:
                raise ValueError(f"{name} must have two entries")
        if any(not math.isfinite(p) or not 0 < p < 1 for p in gp.p_bad):
            raise ValueError("p_bad must lie strictly between zero and one")
        if any(not math.isfinite(g) or g <= 0 for g in gp.gamma):
            raise ValueError("gamma must be finite and strictly positive")
        if any(not math.isfinite(m) or m < 1 for m in gp.bad_multiplier):
            raise ValueError("bad_multiplier must be finite and at least one")
        if any(g not in (0, 1) for g in gp.initial_gates):
            raise ValueError("initial_gates must be binary")
        p = np.asarray(gp.p_bad, dtype=float)
        self.gate_energy_levels = -self.params.theta * np.log(p / (1 - p))
        self.gate_up_rates = np.asarray(gp.gamma) * p / (1 - p)
        self.gate_down_rates = np.asarray(gp.gamma, dtype=float)
        self.base_event_count = len(self.base.events)
        self.gate_events = tuple(GateEvent(g, b, 1 - b) for g in range(2) for b in range(2))
        self.events: tuple[Event | GateEvent, ...] = self.base.events + self.gate_events
        self._event_lookup = {event: k for k, event in enumerate(self.events)}
        self.reverse_index = np.r_[self.base.reverse_index,
                                   self.base_event_count + np.array([1, 0, 3, 2])]

    def initial(self, seed: int = 0) -> EnvironmentState:
        return EnvironmentState(self.base.initial(seed), np.asarray(self.gateparams.initial_gates, dtype=np.uint8))

    @staticmethod
    def copy_state(state: EnvironmentState) -> EnvironmentState:
        return EnvironmentState(BaseModel.copy_state(state.base_state), state.gates.copy(),
                                state.cumulative_gate_heat, state.gate_event_count)

    def environment(self, state: EnvironmentState) -> tuple[float, float]:
        source = (0, 1) if self.response_mode == "differentiated" else (0, 0)
        return tuple(self.gateparams.bad_multiplier[a] if state.gates[source[a]] else 1.0
                     for a in range(2))

    def gate_energy(self, state: EnvironmentState) -> float:
        return float(self.gate_energy_levels @ state.gates)

    def energy(self, state: EnvironmentState) -> float:
        return self.base.energy(state.base_state) + self.gate_energy(state)

    def gate_rate(self, state: EnvironmentState, gate: int) -> float:
        return float(self.gate_down_rates[gate] if state.gates[gate] else self.gate_up_rates[gate])

    def rates(self, state: EnvironmentState) -> tuple[tuple[Event | GateEvent, ...], np.ndarray]:
        _, material_rates = self.base.rates(state.base_state, self.environment(state))
        gates = np.zeros(4)
        for g in range(2):
            gates[2 * g + int(state.gates[g])] = self.gate_rate(state, g)
        return self.events, np.r_[material_rates, gates]

    def event_index(self, event: Event | GateEvent | int) -> int:
        return int(event) if isinstance(event, (int, np.integer)) else self._event_lookup[event]

    def reverse(self, event: Event | GateEvent | int) -> Event | GateEvent:
        return self.events[int(self.reverse_index[self.event_index(event)])]

    def delta_energy(self, state: EnvironmentState, event: Event | GateEvent | int) -> float:
        k = self.event_index(event)
        if k < self.base_event_count:
            return self.base.delta_energy(state.base_state, k)
        e = self.events[k]
        assert isinstance(e, GateEvent)
        return float(self.gate_energy_levels[e.gate] * (e.after - e.before))

    def apply(self, state: EnvironmentState, event: Event | GateEvent | int) -> dict[str, Any]:
        k = self.event_index(event)
        gates_before = state.gates.tolist()
        if k < self.base_event_count:
            record = self.base.apply(state.base_state, k)
            record.update({"subsystem": "material", "delta_H_material": record["delta_H"],
                           "delta_H_environment": 0.0})
        else:
            e = self.events[k]
            assert isinstance(e, GateEvent)
            if int(state.gates[e.gate]) != e.before:
                raise ValueError("gate edge source does not match state")
            delta = self.delta_energy(state, k)
            state.gates[e.gate] = e.after
            state.cumulative_gate_heat -= delta
            state.gate_event_count += 1
            fuel = int(state.base_state.fuel.sum())
            record = {"time": float(state.time), "event_index": k, "kind": "environment_gate",
                      "i": e.gate, "j": -1, "before": e.before, "after": e.after,
                      "delta_H": delta, "heat_to_bath": -delta, "external_work": 0.0,
                      "fuel_before": fuel, "fuel_after": fuel,
                      "fuel_standard_energy_change": 0.0, "material_energy_change": 0.0,
                      "instance_created": 0, "instance_removed": 0, "parent_instance": 0,
                      "catalyst_instance": 0, "subsystem": "environment",
                      "delta_H_material": 0.0, "delta_H_environment": delta}
        record["gates_before"] = gates_before
        record["gates_after"] = state.gates.tolist()
        return record

    def summary(self, state: EnvironmentState) -> dict[str, Any]:
        out = self.base.summary(state.base_state)
        out["H_material_and_fuel"] = out["H"]
        out["H_environment"] = self.gate_energy(state)
        out["H"] += out["H_environment"]
        out["material_heat_to_bath"] = out["heat_to_bath"]
        out["gate_heat_to_bath"] = float(state.cumulative_gate_heat)
        out["heat_to_bath"] += state.cumulative_gate_heat
        out["material_events"] = out["events"]
        out["gate_events"] = state.gate_event_count
        out["events"] += state.gate_event_count
        out["gates"] = state.gates.tolist()
        out["passive_environment_multipliers"] = list(self.environment(state))
        return out

    @staticmethod
    def state_to_dict(state: EnvironmentState) -> dict[str, Any]:
        return {"base_state": BaseModel.state_to_dict(state.base_state), "gates": state.gates.tolist(),
                "cumulative_gate_heat": state.cumulative_gate_heat, "gate_event_count": state.gate_event_count}

    @staticmethod
    def state_from_dict(data: dict[str, Any]) -> EnvironmentState:
        return EnvironmentState(BaseModel.state_from_dict(data["base_state"]),
                                np.asarray(data["gates"], dtype=np.uint8),
                                float(data.get("cumulative_gate_heat", 0)), int(data.get("gate_event_count", 0)))


def simulate_environment(params: Params | None, seed: int, horizon: float,
                         flags: Flags | dict[str, bool] | None = None,
                         gateparams: GateParams | None = None,
                         response_mode: str = "differentiated", *, sample_dt: float = 1.0,
                         initial_state: EnvironmentState | None = None,
                         max_events: int = 2_000_000,
                         max_gate_events: int = 1_000_000) -> dict[str, Any]:
    """Exact joint CTMC via competing autonomous-gate and material clocks.

    Each gate has an independent random stream and pending exponential clock.
    Consequently equal seeds and GateParams give identical autonomous gate
    histories across material controls. After a gate flip, the material waiting
    time is resampled for the new rates; exponential memorylessness makes this
    exactly equivalent to direct SSA on the combined generator.
    """
    if horizon < 0 or sample_dt <= 0:
        raise ValueError("horizon must be nonnegative and sample_dt positive")
    model = EnvironmentModel(params, flags, gateparams, response_mode)
    state = model.initial(seed) if initial_state is None else model.copy_state(initial_state)
    if state.time != 0:
        raise ValueError("initial joint state must start at time zero")
    material_rng = np.random.default_rng(np.random.SeedSequence([seed, 34017]))
    gate_rng = [np.random.default_rng(np.random.SeedSequence([seed, 34018, g])) for g in range(2)]
    next_gate = np.array([float(gate_rng[g].exponential(1 / model.gate_rate(state, g))) for g in range(2)])
    initial = model.state_to_dict(state)
    sample_times = np.r_[np.arange(0.0, horizon, sample_dt), horizon]
    sample_index = 0
    samples: list[dict[str, Any]] = []
    ledger: list[dict[str, Any]] = []
    while state.time < horizon:
        _, rates = model.base.rates(state.base_state, model.environment(state))
        total = float(rates.sum())
        material_time = math.inf if total == 0 else state.time + float(material_rng.exponential(1 / total))
        g = int(np.argmin(next_gate))
        next_time = min(material_time, float(next_gate[g]))
        while sample_index < len(sample_times) and sample_times[sample_index] <= min(next_time, horizon):
            sample = model.summary(state)
            sample["time"] = float(sample_times[sample_index])
            samples.append(sample)
            sample_index += 1
        if next_time > horizon:
            state.time = float(horizon)
            break
        state.time = next_time
        if material_time < next_gate[g]:
            k = int(np.searchsorted(np.cumsum(rates), material_rng.random() * total, side="right"))
        else:
            k = model.base_event_count + 2 * g + int(state.gates[g])
        ledger.append(model.apply(state, k))
        if k >= model.base_event_count:
            next_gate[g] = state.time + float(gate_rng[g].exponential(1 / model.gate_rate(state, g)))
        if len(ledger) >= max_events or state.gate_event_count >= max_gate_events:
            raise RuntimeError("event safeguard exceeded; no truncated result is returned")
    if horizon == 0:
        samples.append(model.summary(state))
    return {"params": asdict(model.params), "flags": asdict(model.flags), "gateparams": asdict(model.gateparams),
            "response_mode": response_mode, "seed": seed, "horizon": horizon,
            "initial_state": initial, "final_state": model.state_to_dict(state),
            "samples": samples, "ledger": ledger, "final_summary": model.summary(state),
            "algorithm": "Exact competing independent gate clocks and material SSA; paired seeds share gate histories.",
            "scope": "Finite thermal barrier-gate idealization; supplied response map and mutation alphabet. No state reset, derived weather apparatus, functional biodiversity advantage, or open-ended evolution is asserted."}


if __name__ == "__main__":
    import json
    result = simulate_environment(Params(side=6, initial_templates=6), seed=7, horizon=10)
    model = EnvironmentModel(Params(**result["params"]), gateparams=GateParams(**result["gateparams"]),
                             response_mode=result["response_mode"])
    initial_h = model.energy(model.state_from_dict(result["initial_state"]))
    print(json.dumps({"initial_H": initial_h, "final": result["final_summary"],
                      "events": len(result["ledger"]), "response_mode": result["response_mode"]}, indent=2))
