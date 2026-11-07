"""Independent additive-energy, gate-clock and event-replay wrapper audit.

Uses the independent scalar Reference in frozen verify_interface.py, never
the production apply() during trajectory replay. Prints JSON without writing
model or result files. Accepts JSON/GZ records, summaries, directories or globs.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import asdict
import glob
import gzip
from hashlib import sha256
import itertools
import json
import math
from pathlib import Path
import random
import sys

sys.dont_write_bytecode = True
import numpy as np
from interface_model import Flags, Params
from autonomous_environment import EnvironmentModel, GateParams
from verify_interface import Reference, close, compare_state


def compare(expected, actual, errors=None):
    for key, value in expected.items():
        assert key in actual, ("missing field", key)
        if isinstance(value, float):
            error = abs(value-float(actual[key]))
            close(value, actual[key])
            if errors is not None:
                errors[key] = max(errors.get(key, 0.), error)
        else:
            assert value == actual[key], (key, value, actual[key])


class JointReference:
    def __init__(self, params, flags, gateparams, mode):
        self.base = Reference(params, flags)
        self.p, self.gp, self.mode = params, gateparams, mode
        self.levels = [-params.theta*math.log(p/(1-p)) for p in gateparams.p_bad]
        self.up = [g*p/(1-p) for g, p in zip(gateparams.gamma, gateparams.p_bad)]
        self.down = list(gateparams.gamma)
        self.base_count = len(self.base.events)

    def environment(self, state):
        indices = (0, 1) if self.mode == "differentiated" else (0, 0)
        return [self.gp.bad_multiplier[a] if state["gates"][indices[a]] else 1. for a in range(2)]

    def energy(self, state):
        return self.base.energy(state["base_state"])+sum(a*b for a, b in zip(self.levels, state["gates"]))

    def rate(self, state, index):
        if index < self.base_count:
            self.base.env = self.environment(state)
            return self.base.rate(state["base_state"], self.base.events[index])
        gate, before = divmod(index-self.base_count, 2)
        return (self.down[gate] if before else self.up[gate]) if state["gates"][gate] == before else 0.

    def reverse_index(self, index):
        if index >= self.base_count:
            return index+1 if (index-self.base_count) % 2 == 0 else index-1
        e = self.base.events[index]
        return self.base.lookup[(e[0], e[1], e[2], e[4], e[3])] if e[0] <= 2 else index

    def apply(self, state, index):
        assert self.rate(state, index) > 0, ("disabled transition", index)
        before_gates = list(state["gates"])
        if index < self.base_count:
            record = self.base.apply(state["base_state"], index)
            record.update(subsystem="material", delta_H_material=record["delta_H"], delta_H_environment=0.)
        else:
            gate, before = divmod(index-self.base_count, 2)
            after = 1-before
            delta = self.levels[gate]*(after-before)
            state["gates"][gate] = after
            state["cumulative_gate_heat"] -= delta
            state["gate_event_count"] += 1
            fuel = sum(state["base_state"]["fuel"])
            record = {"time": float(state["base_state"]["time"]), "event_index": index,
                      "kind": "environment_gate", "i": gate, "j": -1,
                      "before": before, "after": after, "delta_H": delta,
                      "heat_to_bath": -delta, "external_work": 0., "fuel_before": fuel,
                      "fuel_after": fuel, "fuel_standard_energy_change": 0., "material_energy_change": 0.,
                      "instance_created": 0, "instance_removed": 0, "parent_instance": 0,
                      "catalyst_instance": 0, "subsystem": "environment",
                      "delta_H_material": 0., "delta_H_environment": delta}
        record.update(gates_before=before_gates, gates_after=list(state["gates"]))
        return record

    def summary(self, state):
        out = self.base.summary(state["base_state"])
        out["H_material_and_fuel"] = out["H"]
        out["H_environment"] = sum(a*b for a, b in zip(self.levels, state["gates"]))
        out["H"] += out["H_environment"]
        out["material_heat_to_bath"] = out["heat_to_bath"]
        out["gate_heat_to_bath"] = float(state["cumulative_gate_heat"])
        out["heat_to_bath"] += out["gate_heat_to_bath"]
        out["material_events"] = out["events"]
        out["gate_events"] = state["gate_event_count"]
        out["events"] += out["gate_events"]
        out["gates"] = list(state["gates"])
        out["passive_environment_multipliers"] = self.environment(state)
        return out


def compare_joint(expected, actual):
    compare_state(expected["base_state"], actual["base_state"])
    assert expected["gates"] == actual["gates"]
    assert expected["gate_event_count"] == actual["gate_event_count"]
    close(expected["cumulative_gate_heat"], actual["cumulative_gate_heat"])


def clock_history(seed, horizon, gateparams, initial_gates):
    """Reconstruct autonomous gate paths without consulting material events."""
    history = []
    for gate in range(2):
        rng = np.random.default_rng(np.random.SeedSequence([seed, 34018, gate]))
        bit, time = initial_gates[gate], 0.
        while True:
            p, gamma = gateparams.p_bad[gate], gateparams.gamma[gate]
            rate = gamma if bit else gamma*p/(1-p)
            time += float(rng.exponential(1/rate))
            if time > horizon:
                break
            history.append((time, gate, bit, 1-bit))
            bit = 1-bit
    return sorted(history)


def verify_simulation(result):
    p, flags, gp = Params(**result["params"]), Flags(**result["flags"]), GateParams(**result["gateparams"])
    ref = JointReference(p, flags, gp, result["response_mode"])
    state = copy.deepcopy(result["initial_state"])
    ref.base.validate(state["base_state"])
    initial_h = ref.energy(state)
    initial_q = state["base_state"]["cumulative_heat"]+state["cumulative_gate_heat"]
    initial_w = state["base_state"]["external_work"]
    samples, sample_index, errors = result.get("samples", []), 0, {}
    gate_history = []
    direct_h = initial_h

    def check_samples(until):
        nonlocal sample_index
        while sample_index < len(samples) and samples[sample_index]["time"] <= until:
            expected = ref.summary(state)
            expected["time"] = float(samples[sample_index]["time"])
            compare(expected, samples[sample_index], errors)
            sample_index += 1

    for k, row in enumerate(result["ledger"]):
        assert state["base_state"]["time"] <= row["time"] <= result["horizon"]
        check_samples(row["time"])
        state["base_state"]["time"] = row["time"]
        expected = ref.apply(state, int(row["event_index"]))
        compare(expected, row, errors)
        direct_h += expected["delta_H"]
        if expected["subsystem"] == "environment":
            gate_history.append((row["time"], row["i"], row["before"], row["after"]))
        total_q = state["base_state"]["cumulative_heat"]+state["cumulative_gate_heat"]
        close(initial_h+state["base_state"]["external_work"]-initial_w-direct_h-total_q+initial_q, 0.)
        if k % 97 == 0 or k == len(result["ledger"])-1:
            close(ref.energy(state), direct_h)
            ref.base.validate(state["base_state"])
    check_samples(result["horizon"])
    assert sample_index == len(samples)
    state["base_state"]["time"] = float(result["horizon"])
    compare_joint(state, result["final_state"])
    compare(ref.summary(state), result["final_summary"], errors)
    expected_history = clock_history(result["seed"], result["horizon"], gp, result["initial_state"]["gates"])
    assert len(expected_history) == len(gate_history)
    for expected, actual in zip(expected_history, gate_history):
        close(expected[0], actual[0])
        assert expected[1:] == actual[1:]
    return {"seed": result["seed"], "response_mode": result["response_mode"],
            "flags": asdict(flags), "events_replayed": len(result["ledger"]),
            "gate_events_replayed": len(gate_history), "samples_checked": len(samples),
            "initial_joint_H": initial_h, "final_joint_H": ref.energy(state),
            "maximum_record_errors": errors,
            "independent_gate_clock_history_matches": True,
            "gate_history_sha256": sha256(json.dumps(gate_history, separators=(",", ":")).encode()).hexdigest(),
            "initial_physical_state_sha256": sha256(json.dumps({"material": result["initial_state"]["base_state"]["material"],
                         "fuel": result["initial_state"]["base_state"]["fuel"], "gates": result["initial_state"]["gates"]},
                         sort_keys=True).encode()).hexdigest(),
            "common_vs_differentiated_gate_marginals_match":
                gp.p_bad[0] == gp.p_bad[1] and gp.gamma[0] == gp.gamma[1]
                and result["initial_state"]["gates"][0] == result["initial_state"]["gates"][1]}


def verify_model(configurations=6, seed=414812):
    rng = random.Random(seed)
    cases = pairs = propensities = 0
    max_h = max_delta = max_rate = max_log = 0.
    p = Params(side=4, initial_templates=5, theta=.7, J_RR=.6)
    # Unequal gates stress the code; matched-marginal experiment comparisons
    # must use equal p_bad/gamma, as the default gates do.
    gp_settings = [GateParams(), GateParams(p_bad=(.23, .71), gamma=(.2, .7), bad_multiplier=(3., 11.))]
    for gp, mode, production, protection in itertools.product(gp_settings, ("common", "differentiated"), (False, True), (False, True)):
        flags = Flags(production, protection, True)
        model, ref = EnvironmentModel(p, flags, gp, mode), JointReference(p, flags, gp, mode)
        for configuration in range(configurations):
            state = model.state_to_dict(model.initial(seed+configuration))
            b = state["base_state"]
            b["material"] = [rng.randrange(4) for _ in b["material"]]
            b["fuel"] = [rng.randrange(2) for _ in b["fuel"]]
            b["instance_tags"] = [i+1 if a else 0 for i, a in enumerate(b["material"])]
            b["next_tag"] = len(b["material"])+1
            state["gates"] = [configuration % 2, (configuration//2) % 2]
            observed_state = model.state_from_dict(state)
            events, rates = model.rates(observed_state)
            assert len(events) == ref.base_count+4
            expected = np.array([ref.rate(state, k) for k in range(len(events))])
            max_rate = max(max_rate, float(np.max(np.abs(rates-expected))))
            assert np.allclose(rates, expected, rtol=2e-12, atol=2e-13)
            propensities += len(events)
            max_h = max(max_h, abs(ref.energy(state)-model.energy(observed_state)))
            close(ref.energy(state), model.energy(observed_state))
            selected = []
            for kind in range(5):
                candidates = [k for k, e in enumerate(ref.base.events) if e[0] == kind and rates[k] > 0]
                selected.extend(rng.sample(candidates, min(2, len(candidates))))
            selected.extend(k for k in range(ref.base_count, len(events)) if rates[k] > 0)
            for index in selected:
                changed = copy.deepcopy(state)
                expected_record = ref.apply(changed, index)
                actual_state = model.state_from_dict(state)
                actual_record = model.apply(actual_state, index)
                compare(expected_record, actual_record)
                compare_joint(changed, model.state_to_dict(actual_state))
                delta = ref.energy(changed)-ref.energy(state)
                max_delta = max(max_delta, abs(delta-model.delta_energy(observed_state, index)))
                close(delta, expected_record["delta_H"])
                inverse = ref.reverse_index(index)
                assert int(model.reverse_index[index]) == inverse
                assert model.reverse(index) == model.events[inverse]
                reverse = model.rates(model.state_from_dict(changed))[1][inverse]
                assert reverse > 0
                error = abs(math.log(rates[index]/reverse)+delta/p.theta)
                max_log = max(max_log, error)
                assert error < 2e-10
                pairs += 1
            altered = model.state_from_dict(state)
            altered.time = 1e9
            altered.base_state.instance_tags[:] = np.arange(model.base.n)+10000
            altered.base_state.next_tag, altered.base_state.event_count = 999999, 99999
            altered.base_state.cumulative_heat, altered.cumulative_gate_heat = -1e8, 1e8
            altered.gate_event_count = 99999
            assert np.array_equal(model.rates(altered)[1], rates)
            close(model.energy(altered), model.energy(observed_state), 0.)
            # Response maps change kinetics, never the full energy or initial resources.
            other = EnvironmentModel(p, flags, gp, "common" if mode == "differentiated" else "differentiated")
            close(model.energy(observed_state), other.energy(other.state_from_dict(state)), 0.)
            assert model.state_to_dict(model.initial(818)) == other.state_to_dict(other.initial(818))
            cases += 1
    return {"random_state_mode_flag_cases": cases, "individual_propensities_compared": propensities,
            "enabled_reverse_pairs": pairs, "maximum_errors": {"joint_H": max_h, "delta_H": max_delta,
                                                                 "rate": max_rate, "LDB_log": max_log},
            "both_gate_energies_retained_in_both_modes": True,
            "same_mode_initial_resources": True, "observer_tags_and_age_metadata_not_in_rates": True}


def discover(path, seen):
    path = Path(path).resolve()
    if path in seen:
        return
    seen.add(path)
    with (gzip.open if path.suffix == ".gz" else open)(path, "rt", encoding="utf-8") as stream:
        data = json.load(stream)
    def walk(value, label="$"):
        if isinstance(value, dict):
            if {"params", "flags", "gateparams", "response_mode", "initial_state", "ledger"} <= value.keys():
                yield str(path)+":"+label, value
                return
            if isinstance(value.get("trajectory_file"), str):
                yield from discover(path.parent/value["trajectory_file"], seen)
            for key, item in value.items():
                yield from walk(item, label+"."+key)
        elif isinstance(value, list):
            for index, item in enumerate(value):
                yield from walk(item, label+"["+str(index)+"]")
    yield from walk(data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*")
    parser.add_argument("--configurations", type=int, default=6)
    args = parser.parse_args()
    assert args.configurations > 0
    directory = Path(__file__).resolve().parent
    sources = ("interface_model.py", "autonomous_environment.py", "verify_interface.py")
    hashes = {name: sha256((directory/name).read_bytes()).hexdigest() for name in sources}
    checks = verify_model(args.configurations)
    seen, replays = set(), []
    for pattern in args.paths:
        matches = sorted(glob.glob(pattern, recursive=True))
        assert matches, pattern
        for match in matches:
            path = Path(match)
            files = sorted(set(path.rglob("*.json")) | set(path.rglob("*.json.gz"))) if path.is_dir() else [path]
            for file in files:
                for label, result in discover(file, seen):
                    replays.append({"record": label, **verify_simulation(result)})
    assert hashes == {name: sha256((directory/name).read_bytes()).hexdigest() for name in sources}, "audit source changed"
    print(json.dumps({"checks_passed": True, "source_sha256": hashes,
                      "verifier_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
                      "model_checks": checks, "trajectory_records_replayed": len(replays),
                      "total_events_replayed": sum(r["events_replayed"] for r in replays),
                      "total_gate_events_replayed": sum(r["gate_events_replayed"] for r in replays),
                      "replays": replays,
                      "scope": "Additive energetic gate idealization with local reversible material kinetics. This audit does not establish maintained boundaries, functional information, aging, heavy-tail asymptotics, open-ended diversity or agency."}, indent=2))


if __name__ == "__main__":
    main()
