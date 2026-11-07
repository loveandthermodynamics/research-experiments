"""Independent scalar/replay audit of the differentiated driven CTMC.

Only the frozen scalar Reference applies material events during replay. The
production DifferentiatedModel is used only for challenge comparisons and tiny
clock smokes, never as the independent replay implementation.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import asdict, replace
import json
import math
from pathlib import Path
import random
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / "driven-spatial"))
from verify_driven import (DrivenReference, Params, Flags, DriveParams, close,
                           compare_record, compare_state, digest, load)
sys.path.insert(0, str(ROOT))
from differentiated_model import (DifferentiatedModel, ResponseParams,
                                  simulate_differentiated)


class JointReference:
    """Rebuild only the added paired mobility factors and physical gate."""

    def __init__(self, params, flags, drive, response, mode="opposed",
                 transport="local", mobility=None):
        self.base = DrivenReference(params, flags, drive, transport)
        self.response, self.mode = response, mode
        self.mobility = tuple([1.] * self.base.n if mobility is None else mobility)
        self.count = len(self.base.events)
        self.events = self.base.events + [("gate_flip", -1, 1), ("gate_flip", 1, -1)]
        self.reverse_index = self.base.reverse_index + [self.count + 1, self.count]

    def multipliers(self, state):
        first = 1. + self.response.asymmetry * state["gate"]
        return [first, first if self.mode == "aligned" else 2. - first]

    def rate(self, state, index):
        if index >= self.count:
            return self.response.flip_rate if state["gate"] == self.events[index][1] else 0.
        event = self.base.events[index]
        value = self.base.rate(state["base_state"], event)
        kind, _, helper, before, after = event
        if kind == 0 and max(before, after) in (1, 2):
            value *= self.multipliers(state)[max(before, after) - 1]
        if kind == 2:
            value *= self.mobility[helper]
        return value

    def delta(self, state, index):
        return self.base.delta(state["base_state"], self.base.events[index]) if index < self.count else 0.

    def work(self, index):
        return self.base.work(self.base.events[index]) if index < self.count else 0.

    def energy(self, state):
        return self.base.energy(state["base_state"])

    def apply(self, state, index):
        assert self.rate(state, index) > 0, ("disabled event", index)
        gate = state["gate"]
        base = state["base_state"]
        helper = -1
        if index < self.count:
            event = self.base.events[index]
            if event[0] in (1, 2):
                helper = base["material"][event[2]]
            record = self.base.apply(base, index)
        else:
            _, before, after = self.events[index]
            state["gate"] = after
            state["gate_event_count"] += 1
            fuel = sum(base["fuel"])
            record = dict(time=base["time"], event_index=index, kind="gate_flip", i=-1, j=-1,
                          before=before, after=after, delta_H=0., heat_to_bath=0., external_work=0.,
                          chemical_work=0., fuel_before=fuel, fuel_after=fuel,
                          fuel_standard_energy_change=0., material_energy_change=0.,
                          instance_created=0, instance_removed=0, parent_instance=0,
                          catalyst_instance=0, subsystem="environment")
        record.update(gate_before=gate, gate_after=state["gate"], helper_species=helper)
        return record

    def summary(self, state, include_state=False):
        out = self.base.summary(state["base_state"], include_state)
        out["material_reservoir_transport_events"] = out["events"]
        out["events"] += state["gate_event_count"]
        out.update(gate=state["gate"], gate_events=state["gate_event_count"],
                   passive_multipliers=self.multipliers(state))
        return out


def make_reference(run):
    return JointReference(Params(**run["params"]), Flags(**run["flags"]),
                          DriveParams(**run["drive"]), ResponseParams(**run["response"]),
                          run["response_mode"], run["transport_mode"],
                          run["production_site_mobility"])


def verify_trajectory(run):
    ref = make_reference(run)
    close(ref.base.mu, run["effective_mu"])
    state = copy.deepcopy(run["initial_state"])
    base = state["base_state"]
    ref.base.validate(base)
    assert state["gate"] in (-1, 1)
    initial_h = ref.energy(state)
    initial_q, initial_w = base["cumulative_heat"], base["external_work"]
    initial_count, initial_gate_count = base["event_count"], state["gate_event_count"]
    initial_fuel = sum(base["fuel"])
    direct_h, sample_index = initial_h, 0
    samples, errors = run["samples"], {}
    counts, reservoir_balance, internal_fuel = {}, 0, 0

    def samples_until(until):
        nonlocal sample_index
        while sample_index < len(samples) and samples[sample_index]["time"] <= until:
            observed = samples[sample_index]
            expected = ref.summary(state, True)
            expected["time"] = observed["time"]
            compare_record(expected, observed, errors)
            sample_index += 1

    for n, saved in enumerate(run["ledger"]):
        assert base["time"] <= saved["time"] <= run["horizon"]
        samples_until(saved["time"])
        base["time"] = saved["time"]
        index = int(saved["event_index"])
        assert 0 <= index < len(ref.events)
        expected = ref.apply(state, index)
        compare_record(expected, saved, errors)
        counts[saved["kind"]] = counts.get(saved["kind"], 0) + 1
        if saved["kind"] == "reservoir_exchange":
            reservoir_balance += saved["after"] - saved["before"]
        elif saved["kind"] in ("copy_R", "produce_L"):
            internal_fuel += -1 if saved["before"] == 0 else 1
        assert sum(base["fuel"]) == initial_fuel + reservoir_balance + internal_fuel
        close(base["external_work"] - initial_w, ref.base.mu * reservoir_balance)
        direct_h += expected["delta_H"]
        close(initial_h + base["external_work"] - initial_w - direct_h - base["cumulative_heat"] + initial_q, 0.)
        if n % 97 == 0 or n == len(run["ledger"]) - 1:
            close(ref.energy(state), direct_h)
            ref.base.validate(base)
    samples_until(run["horizon"])
    assert sample_index == len(samples)
    base["time"] = run["horizon"]
    compare_state(base, run["final_state"]["base_state"])
    assert state["gate"] == run["final_state"]["gate"]
    assert state["gate_event_count"] == run["final_state"]["gate_event_count"]
    compare_record(ref.summary(state), run["final_summary"], errors)
    assert base["event_count"] - initial_count + state["gate_event_count"] - initial_gate_count == len(run["ledger"])
    close(ref.energy(state), direct_h)
    return dict(seed=run["seed"], flags=run["flags"], response_mode=run["response_mode"],
                response=run["response"], events_replayed=len(run["ledger"]), event_counts=counts,
                samples_checked=len(samples), net_reservoir_fuel_exchange=reservoir_balance,
                internal_fuel_change=internal_fuel, initial_H=initial_h, final_H=ref.energy(state),
                chemical_work=base["external_work"] - initial_w,
                heat_to_bath=base["cumulative_heat"] - initial_q, maximum_record_errors=errors)


def check_physics(seed=770123):
    rng = random.Random(seed)
    p = Params(side=4, initial_templates=4, initial_fuel=12)
    rates_checked = pairs = cases = 0
    max_rate = max_delta = max_ldb = 0.
    initial_states = []
    for mode in ("opposed", "aligned"):
        for produce in (False, True):
            model = DifferentiatedModel(p, Flags(produce_L=produce), response_mode=mode)
            initial_states.append(model.state_to_dict(model.initial(71)))
    assert all(x == initial_states[0] for x in initial_states)
    for configuration in range(2):
        state = copy.deepcopy(initial_states[0])
        b = state["base_state"]
        b["material"] = [rng.randrange(4) for _ in range(p.side ** 2)]
        b["fuel"] = [rng.randrange(2) for _ in b["material"]]
        b["instance_tags"] = [i + 1 if a else 0 for i, a in enumerate(b["material"])]
        b["next_tag"] = len(b["material"]) + 1
        # Deliberately include exact-zero, fractional, and natural helper masks.
        mask = [0. if i % 3 == 0 else .25 if i % 3 == 1 else 1. for i in range(len(b["material"]))]
        for mode in ("opposed", "aligned"):
            for transport in ("local", "mixed_L"):
                for mu in (None, 0.):
                    for produce in (False, True):
                        for gate in (-1, 1):
                            state["gate"] = gate
                            flags, drive = Flags(produce_L=produce), DriveParams(mu=mu)
                            response = ResponseParams(asymmetry=.8, flip_rate=.1)
                            ref = JointReference(p, flags, drive, response, mode, transport, mask)
                            model = DifferentiatedModel(p, flags, drive, response, mode, transport, mask)
                            actual = model.state_from_dict(state)
                            events, observed = model.rates(actual)
                            assert list(model.reverse_index) == ref.reverse_index
                            assert [tuple(asdict(e).values()) for e in events[:ref.count]] == ref.base.events
                            assert [(e.kind, e.before, e.after) for e in events[ref.count:]] == ref.events[ref.count:]
                            close(model.energy(actual), ref.energy(state))
                            enabled = []
                            for index in range(len(ref.events)):
                                expected = ref.rate(state, index)
                                close(expected, observed[index], 2e-12)
                                max_rate = max(max_rate, abs(expected - float(observed[index])))
                                rates_checked += 1
                                if expected > 0:
                                    enabled.append(index)
                            rng.shuffle(enabled)
                            selected = set(enabled[:40] + [ref.count + (gate == 1)])
                            for index in selected:
                                before_h = ref.energy(state)
                                delta, work = ref.delta(state, index), ref.work(index)
                                close(model.delta_energy(actual, index), delta)
                                max_delta = max(max_delta, abs(model.delta_energy(actual, index) - delta))
                                after = copy.deepcopy(state)
                                expected_record = ref.apply(after, index)
                                converted = model.copy_state(actual)
                                observed_record = model.apply(converted, index)
                                compare_record(expected_record, observed_record, {})
                                compare_state(after["base_state"], model.state_to_dict(converted)["base_state"])
                                assert after["gate"] == converted.gate
                                close(ref.energy(after) - before_h, delta)
                                reverse = ref.reverse_index[index]
                                forward, backward = ref.rate(state, index), ref.rate(after, reverse)
                                assert backward > 0
                                error = abs(math.log(forward / backward) - (work - delta) / p.theta)
                                assert error < 2e-11
                                max_ldb = max(max_ldb, error)
                                pairs += 1
                            cases += 1
    # Mutating the caller's mask and changing neutral tags must not alter rates.
    supplied = [1.] * 16
    model = DifferentiatedModel(p, production_site_mobility=supplied)
    original = model.rates(model.initial(9))[1]
    supplied[0] = 0.
    altered = model.initial(9)
    altered.base_state.instance_tags *= 37
    altered.base_state.next_tag *= 37
    assert (model.rates(altered)[1] == original).all()
    assert model.production_site_mobility == tuple([1.] * 16)
    return dict(configurations=2, cases=cases, scalar_rates_checked=rates_checked,
                reverse_pairs_checked=pairs, initial_state_parity=True,
                zero_mu_included=True, zero_and_fractional_helper_masks_included=True,
                neutral_tags_nonfunctional=True, caller_mask_mutation_isolated=True,
                maximum_rate_error=max_rate, maximum_delta_error=max_delta,
                maximum_local_detailed_balance_log_error=max_ldb)


def check_clock_smokes():
    p = Params(side=4, initial_templates=4, initial_fuel=12)
    runs = []
    for mode in ("opposed", "aligned"):
        for produce in (False, True):
            runs.append(simulate_differentiated(p, 442, 8., Flags(produce_L=produce),
                        response=ResponseParams(.8, .7), response_mode=mode, sample_dt=1.))
    paths = [[(e["time"], e["before"], e["after"]) for e in r["ledger"] if e["kind"] == "gate_flip"] for r in runs]
    assert all(path == paths[0] for path in paths)
    assert all(r["initial_state"] == runs[0]["initial_state"] for r in runs)
    audits = [verify_trajectory(run) for run in runs]
    # A realized gate survives a supplied-state restart, including a nondefault sign.
    model = DifferentiatedModel(p, response=ResponseParams(.5, .1, initial_gate=1))
    supplied = model.initial(9)
    supplied.gate = -1
    continuation = simulate_differentiated(p, 19, 0., response=ResponseParams(.5, .1, initial_gate=1), initial_state=supplied)
    assert continuation["initial_state"]["gate"] == continuation["final_state"]["gate"] == -1
    return dict(scope="Four tiny physical clock/replay smokes; not study outcomes or additional ensemble evidence",
                matched_gate_histories=True, gate_flips_per_arm=len(paths[0]),
                initial_state_parity=True, continuation_gate_preserved=True,
                trajectories=audits)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = dict(verifier_sha256=digest(__file__), source_sha256={
        "differentiated_model.py": digest(ROOT / "differentiated_model.py"),
        "verify_driven.py": digest(ROOT.parent / "driven-spatial" / "verify_driven.py"),
        "verify_interface.py": digest(ROOT.parent / "interface-maintenance" / "verify_interface.py")},
        physics=check_physics(), clock_smokes=check_clock_smokes(),
        supplied_trajectories=[dict(file=str(path), sha256=digest(path), **verify_trajectory(load(path))) for path in args.paths],
        scope="Conditional thermodynamic/rate/replay checks only. No proof of diversity benefit, selection, collective individuality, or lower resource cost.")
    Path(args.output).write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"output": args.output, "physics": result["physics"],
                      "matched_gate_histories": result["clock_smokes"]["matched_gate_histories"],
                      "supplied_trajectories": len(result["supplied_trajectories"])}, indent=2))


if __name__ == "__main__":
    main()
