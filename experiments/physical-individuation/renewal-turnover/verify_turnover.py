"""Independent audit of a bounded protective-material turnover study.

This file owns no simulation source. It imports the already independent scalar
reference, rather than the production replay/readout implementations.
"""
from __future__ import annotations

import copy
import argparse
from dataclasses import asdict
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path
import random
import statistics
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
BASE = ROOT.parent / "interface-maintenance"
sys.path.insert(0, str(BASE))
REFERENCE_SHA256 = "ccbaac8114d5f4ad48cf1bd655de7c186ddca0f4e99197322ed3fc4e43aef0b6"
MODEL_SHA256 = "08c5d9948e84d62ddc77a0265e5e91e805b2badf47bc516542a36fa918201927"
PROTOCOL_SHA256 = {
    "pilot": "c69ff33639bc804d50632222c9540353674d24fc003a36383ae757abf38351cf",
    "confirmation": "1d4d6e2b11b01ca7481b54907f15e22017116c3e232ff74d0a8910b8c25b5259",
}

from verify_interface import Reference, Params, Flags, Model, close, compare_record, compare_state, verify_simulation


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    path = Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as stream:
            return json.load(stream)
    return json.loads(path.read_text())


def estimate(values):
    mean = statistics.mean(values)
    se = statistics.stdev(values) / math.sqrt(len(values))
    return dict(n=len(values), mean=mean, monte_carlo_standard_error=se,
                approximate_95_percent_normal_interval=[mean-1.96*se, mean+1.96*se])


def exposures(run):
    """Independent time integral from event deltas, including zero-R starts."""
    initial, final = run["initial_state"], run["final_state"]
    counts = [initial["material"].count(k) for k in (1, 2, 3)]
    area = [0., 0., 0.]
    previous = 0.
    first_zero = 0. if counts[0]+counts[1] == 0 else None
    events = dict(paid_R_births=0, reverse_R_copies=0, paid_L_births=0,
                  reverse_L_production=0, passive_R_births=0, passive_R_deaths=0,
                  passive_L_births=0, passive_L_deaths=0)
    for event in run["ledger"]:
        t = event["time"]
        assert previous <= t <= run["horizon"]
        for j in range(3):
            area[j] += counts[j]*(t-previous)
        previous = t
        if event["kind"] in ("material_swap", "fuel_swap"):
            continue
        before, after = event["before"], event["after"]
        for j, species in enumerate((1, 2, 3)):
            counts[j] += int(after == species)-int(before == species)
        if event["kind"] == "copy_R":
            key = "paid_R_births" if after in (1, 2) else "reverse_R_copies"
        elif event["kind"] == "produce_L":
            key = "paid_L_births" if after == 3 else "reverse_L_production"
        elif before in (1, 2) or after in (1, 2):
            key = "passive_R_births" if after in (1, 2) else "passive_R_deaths"
        else:
            key = "passive_L_births" if after == 3 else "passive_L_deaths"
        events[key] += 1
        if first_zero is None and counts[0]+counts[1] == 0:
            first_zero = t
    for j in range(3):
        area[j] += counts[j]*(run["horizon"]-previous)
    assert counts == [final["material"].count(k) for k in (1, 2, 3)]
    fuel = sum(initial["fuel"])-sum(final["fuel"])
    assert fuel == (events["paid_R_births"]+events["paid_L_births"]
                    -events["reverse_R_copies"]-events["reverse_L_production"])
    return dict(integral_R=area[0]+area[1], integral_R0=area[0],
                integral_R1=area[1], integral_L=area[2],
                final_R=counts[0]+counts[1], final_L=counts[2],
                first_R_extinction=first_zero,
                first_episode_duration=run["horizon"] if first_zero is None else first_zero,
                net_fuel_consumed=fuel, net_heat_to_bath=final["cumulative_heat"], **events)


def compare_readouts(run, row):
    calculated = exposures(run)
    zero_initial_R = not any(a in (1, 2) for a in run["initial_state"]["material"])
    assert row["initial_absence_episode_readout_corrected"] == zero_initial_R
    if zero_initial_R:
        assert "frozen_helper_first_episode_duration" in row
    for field, value in calculated.items():
        assert field in row, ("missing readout", field)
        if value is None:
            assert row[field] is None
        else:
            close(value, row[field])
    return calculated


def check_actual_parameters(params, configurations=3, seed=1987301):
    """Challenge actual study parameters independently of saved trajectories."""
    rng = random.Random(seed)
    p = Params(**params)
    maxima = dict(H=0., delta_H=0., rate=0., LDB_log=0.)
    rate_count = pair_count = 0
    arms = [Flags(produce_L=a, protect_R=b) for a, b in
            ((True, True), (True, False), (False, True), (False, False))]
    initial = [Model(p, flags).state_to_dict(Model(p, flags).initial(171)) for flags in arms]
    assert all(state == initial[0] for state in initial)
    for _ in range(configurations):
        state = copy.deepcopy(initial[0])
        state["material"] = [rng.randrange(4) for _ in state["material"]]
        state["fuel"] = [rng.randrange(2) for _ in state["fuel"]]
        state["instance_tags"] = [j+1 if a else 0 for j, a in enumerate(state["material"])]
        state["next_tag"] = len(state["material"])+1
        for flags in arms:
            ref, model = Reference(p, flags), Model(p, flags)
            converted = model.state_from_dict(state)
            observed_events, rates = model.rates(converted)
            assert [tuple(asdict(event).values()) for event in observed_events] == ref.events
            maxima["H"] = max(maxima["H"], abs(model.energy(converted)-ref.energy(state)))
            close(model.energy(converted), ref.energy(state))
            enabled = []
            for index, event in enumerate(ref.events):
                expected = ref.rate(state, event)
                maxima["rate"] = max(maxima["rate"], abs(expected-float(rates[index])))
                close(expected, rates[index], 2e-12)
                rate_count += 1
                if expected > 0:
                    enabled.append(index)
            rng.shuffle(enabled)
            for index in enabled[:80]:
                event = ref.events[index]
                delta = ref.delta(state, event)
                actual = model.delta_energy(converted, observed_events[index])
                maxima["delta_H"] = max(maxima["delta_H"], abs(delta-actual))
                close(delta, actual)
                after = copy.deepcopy(state)
                ref.apply(after, index)
                reverse = ref.events[int(model.reverse_index[index])]
                forward_rate, reverse_rate = ref.rate(state, event), ref.rate(after, reverse)
                assert reverse_rate > 0
                discrepancy = abs(math.log(forward_rate/reverse_rate)+delta/p.theta)
                maxima["LDB_log"] = max(maxima["LDB_log"], discrepancy)
                assert discrepancy < 2e-11
                close(ref.energy(after)-ref.energy(state), delta)
                pair_count += 1
    return dict(passive_L=p.passive_L, state_flag_cases=configurations*len(arms),
                individual_propensities_compared=rate_count, reverse_pairs=pair_count,
                maximum_errors=maxima, factorial_initial_states_identical=True)


WEIGHTS = {
    "production_effect_when_protected": {(True, True): 1, (False, True): -1},
    "production_effect_without_protection": {(True, False): 1, (False, False): -1},
    "protection_effect_with_production": {(True, True): 1, (True, False): -1},
    "protection_effect_after_production_stopped": {(False, True): 1, (False, False): -1},
    "production_protection_interaction": {(True, True): 1, (True, False): -1,
                                          (False, True): -1, (False, False): 1},
}
FULL = dict(produce_L=True, protect_R=True, copy_R=True)


def compare_estimate(expected, observed):
    assert expected["n"] == observed["n"]
    for field in ("mean", "monte_carlo_standard_error"):
        close(expected[field], observed[field])
    for a, b in zip(expected["approximate_95_percent_normal_interval"],
                    observed["approximate_95_percent_normal_interval"]):
        close(a, b)


def turnover(run, original_founders):
    """Compute clearance by removal-time maxima, independently of the runner."""
    initial = run["initial_state"]
    tags = set(t for t in initial["instance_tags"] if t)
    removed = {event["instance_removed"]: event["time"] for event in run["ledger"]
               if event["instance_removed"]}
    groups = {"original_founder": set(original_founders) & tags,
              "preparation_R": {t for a, t in zip(initial["material"], initial["instance_tags"])
                                if a in (1, 2)},
              "preparation_L": {t for a, t in zip(initial["material"], initial["instance_tags"])
                                if a == 3}}
    out = {}
    for name, group in groups.items():
        clearance = (max((removed[t] for t in group), default=0.)
                     if all(t in removed for t in group) else None)
        count = sum(event["kind"] == "copy_R" and event["after"] in (1, 2)
                    and clearance is not None and event["time"] >= clearance
                    for event in run["ledger"])
        out[name+"_clearance_time_post_branch"] = clearance
        out["paid_R_births_after_"+name+"_clearance"] = count
    return out


def check_selection(protocol):
    if protocol["stage"] != "confirmation":
        return None
    pilot_path = ROOT / protocol["pilot_results_file"]
    assert digest(pilot_path) == protocol["pilot_results_sha256"]
    pilot = load(pilot_path)
    pilot_protocol_path = ROOT / "pilot-protocol.json"
    assert digest(pilot_protocol_path) == protocol["pilot_protocol_sha256"]
    pilot_protocol = load(pilot_protocol_path)
    assert pilot["protocol_sha256"] == digest(pilot_protocol_path)
    assert protocol["seeds"] == pilot_protocol["confirmation_seeds"]
    assert not set(protocol["seeds"]) & set(pilot_protocol["seeds"])
    assert protocol["base_params"] == pilot_protocol["base_params"]
    assert protocol["source_hashes"] == pilot_protocol["source_hashes"]
    for field in ("preparation_horizon", "continuation_horizon", "arms", "post_seed_offset",
                  "selection_rule", "primary_contrasts", "joint_support_rule"):
        assert protocol[field] == pilot_protocol[field]
    scores = []
    for level in pilot_protocol["passive_L_levels"]:
        by_seed = {}
        for row in pilot["runs"]:
            if row["passive_L"] == level:
                arm = row["flags"]["produce_L"], row["flags"]["protect_R"]
                by_seed.setdefault(row["origin_seed"], {})[arm] = row["integral_R"]
        primary = []
        for name in pilot_protocol["primary_contrasts"]:
            values = [sum(weight*arms[arm] for arm, weight in WEIGHTS[name].items())
                      for _, arms in sorted(by_seed.items())]
            primary.append(estimate(values))
        zs = [stat["mean"]/stat["monte_carlo_standard_error"]
              if stat["monte_carlo_standard_error"] else
              (1e300 if stat["mean"] > 0 else -1e300 if stat["mean"] < 0 else 0.)
              for stat in primary]
        scores.append(dict(passive_L=level, selection_score=min(zs), primary_z_scores=zs))
    selected = max(scores, key=lambda row: (row["selection_score"], -row["passive_L"]))["passive_L"]
    assert protocol["passive_L_levels"] == [selected]
    assert protocol["selected_passive_L"] == selected
    for expected, observed in zip(scores, protocol["selection_scores"]):
        assert expected["passive_L"] == observed["passive_L"]
        close(expected["selection_score"], observed["selection_score"])
        for a, b in zip(expected["primary_z_scores"], observed["primary_z_scores"]):
            close(a, b)
    return dict(checks_passed=True, selected_passive_L=selected,
                scores_recomputed=scores, fresh_confirmation_seeds=True)


def verify_study(path):
    result = load(path)
    protocol_path = ROOT / result["protocol_file"]
    protocol = load(protocol_path)
    assert digest(protocol_path) == result["protocol_sha256"]
    assert digest(protocol_path) == PROTOCOL_SHA256[result["stage"]]
    assert protocol["source_hashes"] == result["source_hashes"]
    for name, expected in protocol["source_hashes"].items():
        assert digest(ROOT/name) == expected, ("changed source", name)
    assert digest(ROOT/protocol["base_parameter_protocol"]) == protocol["base_parameter_protocol_sha256"]
    assert load(ROOT/protocol["base_parameter_protocol"])["params"] == protocol["base_params"]
    for field in ("stage", "seeds", "passive_L_levels", "preparation_horizon"):
        assert result[field] == protocol[field]
    assert result["horizon"] == protocol["continuation_horizon"]
    assert protocol["preparation_horizon"] == 20. and protocol["continuation_horizon"] == 80.
    assert set(protocol["primary_contrasts"]) == {"production_effect_when_protected", "production_protection_interaction"}
    assert protocol["primary_field"] == "integral_R"
    levels, seeds = protocol["passive_L_levels"], protocol["seeds"]
    assert len(set(levels)) == len(levels) and len(set(seeds)) == len(seeds)
    assert result["replicates_per_arm_per_level"] == len(seeds)
    arms = {(f["produce_L"], f["protect_R"]) for f in protocol["arms"]}
    assert arms == {(True, True), (True, False), (False, True), (False, False)}
    assert all(f["copy_R"] for f in protocol["arms"])
    prep_index = {(row["passive_L"], row["origin_seed"]): row for row in result["preparations"]}
    expected_preps = {(level, seed) for level in levels for seed in seeds}
    assert set(prep_index) == expected_preps
    assert len(prep_index) == len(result["preparations"])
    continuation_index = {(row["passive_L"], row["origin_seed"], row["flags"]["produce_L"], row["flags"]["protect_R"]): row
                          for row in result["runs"]}
    assert set(continuation_index) == {(level, seed, *arm) for level, seed in expected_preps for arm in arms}
    assert len(continuation_index) == len(result["runs"])
    replays, verified_values, initial_by_seed = [], {}, {}
    no_r_preparations = no_l_preparations = 0
    for level, seed in sorted(expected_preps):
        prep_row = prep_index[(level, seed)]
        prep_path = ROOT/prep_row["trajectory_file"]
        assert digest(prep_path) == prep_row["trajectory_sha256"]
        prep = load(prep_path)
        params = {**protocol["base_params"], "passive_L": level}
        assert prep["params"] == params and prep["flags"] == FULL
        assert prep["seed"] == prep_row["seed"] == seed
        assert prep["horizon"] == 20. and prep["phase"] == "preparation"
        assert prep["stage"] == result["stage"] and prep["passive_L"] == level
        assert prep["source_hashes"] == protocol["source_hashes"]
        assert prep["protocol_sha256"] == digest(protocol_path)
        initial = prep["initial_state"]
        assert initial["material"].count(1) == params["initial_templates"]
        assert initial["material"].count(2) == initial["material"].count(3) == 0
        assert sum(initial["fuel"]) == len(initial["fuel"])
        if seed in initial_by_seed:
            assert initial == initial_by_seed[seed], "rate levels change the microscopic founders/resources"
        initial_by_seed[seed] = initial
        replay = verify_simulation(prep)
        replay.update(record=str(prep_path), stage=result["stage"], phase="preparation", passive_L=level)
        replays.append(replay)
        prep_values = compare_readouts(prep, prep_row)
        no_r_preparations += int(prep_values["final_R"] == 0)
        no_l_preparations += int(prep_values["final_L"] == 0)
        reset = copy.deepcopy(prep["final_state"])
        for key in ("time", "cumulative_heat", "external_work", "event_count"):
            reset[key] = 0
        founders = {tag for species, tag in zip(initial["material"], initial["instance_tags"]) if species in (1, 2)}
        for arm in sorted(arms):
            row = continuation_index[(level, seed, *arm)]
            branch_path = ROOT/row["trajectory_file"]
            assert digest(branch_path) == row["trajectory_sha256"]
            run = load(branch_path)
            flags = dict(produce_L=arm[0], protect_R=arm[1], copy_R=True)
            assert run["params"] == params and run["flags"] == row["flags"] == flags
            assert run["initial_state"] == reset, "physical state differs after branching"
            assert run["seed"] == row["seed"] == row["post_seed"] == seed+protocol["post_seed_offset"]
            assert run["origin_seed"] == seed and run["passive_L"] == level
            assert run["stage"] == result["stage"] and run["phase"] == "continuation"
            assert run["preparation_time"] == 20. and run["horizon"] == 80.
            assert run["source_hashes"] == protocol["source_hashes"]
            assert run["protocol_sha256"] == digest(protocol_path)
            assert run["preparation_file"] == row["preparation_file"] == prep_row["trajectory_file"]
            assert run["preparation_sha256"] == row["preparation_sha256"] == digest(prep_path)
            replay = verify_simulation(run)
            replay.update(record=str(branch_path), stage=result["stage"], phase="continuation", passive_L=level)
            replays.append(replay)
            measured = compare_readouts(run, row)
            measured.update(turnover(run, founders))
            for field, expected in measured.items():
                if expected is None:
                    assert row[field] is None
                else:
                    close(expected, row[field])
            assert row["prepared_R"] == prep_values["final_R"]
            assert row["prepared_L"] == prep_values["final_L"]
            assert row["prepared_fuel"] == sum(reset["fuel"])
            close(row["preparation_net_fuel_consumed"], prep_values["net_fuel_consumed"])
            close(row["preparation_heat_to_bath"], prep_values["net_heat_to_bath"])
            for total, field in (("total_net_fuel_consumed", "net_fuel_consumed"),
                                 ("total_heat_to_bath", "net_heat_to_bath"), ("total_integral_R", "integral_R")):
                expected = prep_values[field]+measured[field]
                close(row[total], expected)
                measured[total] = expected
            ref = Reference(Params(**params), Flags(**flags))
            close(ref.energy(initial)-ref.energy(run["final_state"]), measured["total_heat_to_bath"])
            assert sum(initial["fuel"])-sum(run["final_state"]["fuel"]) == measured["total_net_fuel_consumed"]
            verified_values[(level, seed, *arm)] = measured
    for group in result["groups"]:
        arm = group["flags"]["produce_L"], group["flags"]["protect_R"]
        for field, observed in group["statistics"].items():
            compare_estimate(estimate([verified_values[(group["passive_L"], seed, *arm)][field]
                                       for seed in seeds]), observed)
    verified_contrasts = []
    for contrast in result["paired_contrasts"]:
        weights = WEIGHTS[contrast["contrast"]]
        assert {(w["produce_L"], w["protect_R"]): w["weight"] for w in contrast["definition"]} == weights
        values = [sum(weight*verified_values[(contrast["passive_L"], seed, *arm)][contrast["field"]]
                      for arm, weight in weights.items()) for seed in seeds]
        expected = estimate(values)
        compare_estimate(expected, contrast)
        if contrast["field"] == "integral_R":
            verified_contrasts.append(dict(passive_L=contrast["passive_L"], contrast=contrast["contrast"], **expected))
    joint = None
    if result["stage"] == "confirmation":
        primary = [c for c in verified_contrasts if c["contrast"] in protocol["primary_contrasts"]]
        assert len(primary) == 2
        joint = all(c["approximate_95_percent_normal_interval"][0] > 0 for c in primary)
        assert result["declared_joint_claim_supported"] == joint
        assert len(result["primary_confirmation_results"]) == 2
        for expected in primary:
            observed = next(c for c in result["primary_confirmation_results"] if c["contrast"] == expected["contrast"])
            compare_estimate(expected, observed)
    return dict(checks_passed=True, stage=result["stage"], result_file=str(path), result_sha256=digest(path),
                protocol_sha256=digest(protocol_path), preparation_trajectories=len(prep_index),
                continuation_trajectories=len(continuation_index),
                trajectory_records_replayed=len(replays), total_events_replayed=sum(r["events_replayed"] for r in replays),
                samples_checked=sum(r["samples_checked"] for r in replays),
                no_R_preparations_included=no_r_preparations, no_L_preparations_included=no_l_preparations,
                all_declared_seed_level_arm_cells_present=True, all_branch_initial_physical_states_identical=True,
                preparation_and_continuation_costs_counted=True, readouts_and_all_contrasts_recomputed=True,
                primary_R_time_contrasts=verified_contrasts, declared_joint_claim_supported=joint,
                selection_audit=check_selection(protocol), replays=replays)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", nargs="+", type=Path)
    args = parser.parse_args()
    verifier_hash = digest(Path(__file__))
    assert digest(BASE/"verify_interface.py") == REFERENCE_SHA256
    assert digest(BASE/"interface_model.py") == MODEL_SHA256
    studies = [verify_study(path.resolve()) for path in args.results]
    records = [record["record"] for study in studies for record in study["replays"]]
    assert len(records) == len(set(records)), "duplicate physical records"
    actual_params = {}
    for path in args.results:
        result = load(path)
        protocol = load(ROOT/result["protocol_file"])
        for level in protocol["passive_L_levels"]:
            actual_params[level] = {**protocol["base_params"], "passive_L": level}
    physics = [check_actual_parameters(params) for _, params in sorted(actual_params.items())]
    assert digest(Path(__file__)) == verifier_hash
    assert digest(BASE/"verify_interface.py") == REFERENCE_SHA256
    assert digest(BASE/"interface_model.py") == MODEL_SHA256
    report = dict(checks_passed=True, verifier_sha256=verifier_hash,
                  independent_reference_sha256=digest(BASE/"verify_interface.py"),
                  model_sha256=digest(BASE/"interface_model.py"),
                  trajectory_records_replayed=len(records),
                  preparation_trajectories=sum(s["preparation_trajectories"] for s in studies),
                  continuation_trajectories=sum(s["continuation_trajectories"] for s in studies),
                  total_events_replayed=sum(s["total_events_replayed"] for s in studies),
                  model_checks=physics, studies=studies,
                  scope="A symmetric kinetic-turnover intervention in supplied reversible local chemistry; no reward, group readout, imposed aging clock or replenishment. Preparation plus continuation cost and physical state parity are verified. Any positive contrast is conditional material maintenance, not a membrane, functional memory, Darwinian individuality or universal biodiversity theorem.")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
