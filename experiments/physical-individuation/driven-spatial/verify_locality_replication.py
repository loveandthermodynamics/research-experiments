"""Independent full replay and paired analysis of the fixed locality replication.

Imports frozen scalar physics/readouts; never imports replication simulation or
analysis code. Recomputes intervals from the fresh 256 pairs only.
"""
from __future__ import annotations

import copy
import json
import math
from pathlib import Path
import statistics
import sys

sys.dont_write_bytecode = True
import numpy as np
from scipy.stats import t

from verify_driven import (ROOT, BASE, digest, load, verify_trajectory, readouts,
                           compare_values, close, FULL, REFERENCE_SHA256,
                           BASE_MODEL_SHA256, DRIVEN_MODEL_SHA256)

FROZEN_VERIFIER_SHA256 = "366a0a19b4a0f6c3ca8c89880f223eb5a5a6fa6be3c1c4df8ef59a6454a2238e"
PROTOCOL_SHA256 = "9054f64f239a29e881585a2da5050e8e92cfb4f7ff75026331ad8976e3184c94"


def percentile(sorted_values, probability):
    position = (len(sorted_values)-1)*probability
    low = math.floor(position)
    high = min(low+1, len(sorted_values)-1)
    return sorted_values[low]+(position-low)*(sorted_values[high]-sorted_values[low])


def main():
    own_hash = digest(Path(__file__))
    assert digest(ROOT/"verify_driven.py") == FROZEN_VERIFIER_SHA256
    assert digest(BASE/"verify_interface.py") == REFERENCE_SHA256
    assert digest(BASE/"interface_model.py") == BASE_MODEL_SHA256
    assert digest(ROOT/"driven_model.py") == DRIVEN_MODEL_SHA256
    result_path = ROOT/"locality-replication-results.json"
    result = load(result_path)
    protocol_path = ROOT/result["protocol_file"]
    assert digest(protocol_path) == result["protocol_sha256"] == PROTOCOL_SHA256
    protocol = load(protocol_path)
    assert protocol["stage"] == result["stage"] == "locality_replication"
    assert result["source_hashes"] == protocol["source_hashes"]
    for name, expected in protocol["source_hashes"].items():
        assert digest(ROOT/name) == expected
    for field in ("prior_confirmation_protocol", "prior_confirmation_results"):
        assert digest(ROOT/protocol[field]) == protocol[field+"_sha256"]
    previous_protocol = load(ROOT/protocol["prior_confirmation_protocol"])
    previous = load(ROOT/protocol["prior_confirmation_results"])
    assert protocol["params"] == {**previous_protocol["base_params"], "copy_rate": previous_protocol["selected_copy_rate"]}
    assert protocol["params"]["copy_rate"] == .028
    assert protocol["drive"] == previous_protocol["drive"]
    seeds = list(range(1040001, 1040257))
    assert protocol["seeds"] == result["seeds"] == seeds
    assert not set(seeds) & set(previous_protocol["seeds"])
    assert not set(seeds) & set(load(ROOT/"pilot-protocol.json")["seeds"])
    assert protocol["flags"] == FULL
    assert protocol["transport_modes"] == ["local", "mixed_L"]
    assert protocol["preparation_mode"] == "local"
    assert protocol["preparation_horizon"] == result["preparation_horizon"] == 20.
    assert protocol["continuation_horizon"] == result["horizon"] == 80.
    assert protocol["replicates_per_arm"] == 256
    assert protocol["bootstrap_seed"] == protocol["bootstrap_replicates"] == 10000
    preps = {r["origin_seed"]: r for r in result["preparations"]}
    branches = {(r["origin_seed"], r["transport_mode"]): r for r in result["runs"]}
    pairs = {r["origin_seed"]: r for r in result["paired_differences"]}
    assert len(preps) == len(result["preparations"]) == 256 and set(preps) == set(seeds)
    assert len(branches) == len(result["runs"]) == 512
    assert set(branches) == {(s, m) for s in seeds for m in protocol["transport_modes"]}
    assert len(pairs) == len(result["paired_differences"]) == 256 and set(pairs) == set(seeds)
    replays, observed, differences, formation = [], {}, [], []
    no_r = no_l = 0
    for seed in seeds:
        row = preps[seed]
        path = ROOT/row["trajectory_file"]
        assert digest(path) == row["trajectory_sha256"]
        prep = load(path)
        assert prep["params"] == protocol["params"] and prep["drive"] == protocol["drive"]
        assert prep["flags"] == FULL and prep["transport_mode"] == "local"
        assert prep["seed"] == row["seed"] == seed and prep["horizon"] == 20.
        assert prep["stage"] == "locality_replication" and prep["phase"] == "preparation"
        assert prep["source_hashes"] == protocol["source_hashes"]
        assert prep["protocol_sha256"] == PROTOCOL_SHA256
        initial = prep["initial_state"]
        assert all(initial[key] == 0 for key in ("time", "cumulative_heat", "external_work", "event_count"))
        assert initial["material"].count(1) == 8 and initial["material"].count(2) == initial["material"].count(3) == 0
        assert sum(initial["fuel"]) == 58
        replay = verify_trajectory(prep)
        replay.update(record=str(path), phase="preparation", origin_seed=seed)
        replays.append(replay)
        prep_measured, first_l = readouts(prep)
        compare_values(prep_measured, row)
        no_r += int(prep_measured["final_R"] == 0)
        no_l += int(prep_measured["final_L"] == 0)
        formation.append(dict(origin_seed=seed, initial_L=0, first_L_birth=first_l,
                              final_R=prep_measured["final_R"], final_L=prep_measured["final_L"]))
        close(prep_measured["initial_H"], 760.)
        reset = copy.deepcopy(prep["final_state"])
        for key in ("time", "cumulative_heat", "external_work", "event_count"):
            reset[key] = 0
        founders = {tag for a, tag in zip(initial["material"], initial["instance_tags"]) if a in (1, 2)}
        for mode in protocol["transport_modes"]:
            branch_row = branches[(seed, mode)]
            branch_path = ROOT/branch_row["trajectory_file"]
            assert digest(branch_path) == branch_row["trajectory_sha256"]
            run = load(branch_path)
            assert run["params"] == protocol["params"] and run["drive"] == protocol["drive"]
            assert run["flags"] == branch_row["flags"] == FULL and run["transport_mode"] == mode
            assert run["initial_state"] == reset
            assert run["seed"] == branch_row["seed"] == seed+protocol["post_seed_offset"]
            assert run["origin_seed"] == seed and run["horizon"] == 80.
            assert run["stage"] == "locality_replication" and run["phase"] == "continuation"
            assert run["source_hashes"] == protocol["source_hashes"]
            assert run["protocol_sha256"] == PROTOCOL_SHA256
            assert run["preparation_file"] == branch_row["preparation_file"] == row["trajectory_file"]
            assert run["preparation_sha256"] == branch_row["preparation_sha256"] == digest(path)
            replay = verify_trajectory(run)
            replay.update(record=str(branch_path), phase="continuation", origin_seed=seed)
            replays.append(replay)
            values, _ = readouts(run, founders)
            compare_values(values, branch_row)
            assert branch_row["prepared_R"] == prep_measured["final_R"]
            assert branch_row["prepared_L"] == prep_measured["final_L"]
            assert branch_row["prepared_fuel"] == prep_measured["final_fuel"]
            close(branch_row["supplied_initial_H"], prep_measured["initial_H"])
            close(branch_row["preparation_chemical_work"], prep_measured["chemical_work"])
            close(branch_row["preparation_heat_to_bath"], prep_measured["heat_to_bath"])
            for total, field in (("total_chemical_work", "chemical_work"), ("total_heat_to_bath", "heat_to_bath")):
                values[total] = prep_measured[field]+values[field]
                close(values[total], branch_row[total])
            close(760.+values["total_chemical_work"]-values["final_H"]-values["total_heat_to_bath"], 0.)
            observed[(seed, mode)] = values
        delta = observed[(seed, "local")]["integral_R"]-observed[(seed, "mixed_L")]["integral_R"]
        compare_values(dict(local_integral_R=observed[(seed, "local")]["integral_R"],
                            mixed_L_integral_R=observed[(seed, "mixed_L")]["integral_R"], difference=delta), pairs[seed])
        differences.append(delta)
    assert len({r["record"] for r in replays}) == len(replays) == 768
    for group in result["groups"]:
        for field, mean in group["means"].items():
            close(statistics.mean(observed[(seed, group["transport_mode"])][field] for seed in seeds), mean)
    n = len(differences)
    mean, se = statistics.mean(differences), statistics.stdev(differences)/math.sqrt(n)
    critical = float(t.ppf(.975, n-1))
    interval = [mean-critical*se, mean+critical*se]
    rng = np.random.default_rng(10000)
    indices = rng.integers(0, n, size=(10000, n))
    bootstrap_means = sorted(math.fsum(differences[int(j)] for j in row)/n for row in indices)
    bootstrap_interval = [percentile(bootstrap_means, q) for q in (.025, .975)]
    calculated = dict(n=n, mean=mean, monte_carlo_standard_error=se, degrees_freedom=n-1,
                      t_critical=critical, t_95_interval=interval,
                      primary_support_criterion_met=interval[0] > 0,
                      bootstrap_replicates=10000, bootstrap_seed=10000,
                      percentile_bootstrap_95_interval=bootstrap_interval)
    compare_values(calculated, result["replication_estimate"])
    old_locality = next(c for c in previous["primary_confirmation_results"] if c["contrast"] == "locality_effect_full")
    historical = result["previous_confirmation"]
    assert historical["normal_result_as_saved"] == old_locality
    assert historical["original_joint_criterion_supported_as_saved"] == previous["declared_joint_claim_supported"]
    assert historical["source_sha256"] == protocol["prior_confirmation_results_sha256"]
    old_critical = float(t.ppf(.975, old_locality["n"]-1))
    old_interval = [old_locality["mean"]-old_critical*old_locality["monte_carlo_standard_error"],
                    old_locality["mean"]+old_critical*old_locality["monte_carlo_standard_error"]]
    compare_values(dict(t_sensitivity_degrees_freedom=63, t_critical=old_critical,
                        t_95_interval_sensitivity=old_interval), historical)
    assert result["pooled_estimate"] == "Not computed or authorized."
    assert digest(Path(__file__)) == own_hash
    assert digest(ROOT/"verify_driven.py") == FROZEN_VERIFIER_SHA256
    assert digest(BASE/"verify_interface.py") == REFERENCE_SHA256
    assert digest(BASE/"interface_model.py") == BASE_MODEL_SHA256
    assert digest(protocol_path) == PROTOCOL_SHA256
    for name, expected in protocol["source_hashes"].items():
        assert digest(ROOT/name) == expected
    print(json.dumps(dict(checks_passed=True, verifier_sha256=own_hash,
                          frozen_physics_verifier_sha256=FROZEN_VERIFIER_SHA256,
                          result_sha256=digest(result_path), protocol_sha256=PROTOCOL_SHA256,
                          preparation_trajectories=256, continuation_trajectories=512,
                          trajectory_records_replayed=len(replays), total_events_replayed=sum(r["events_replayed"] for r in replays),
                          samples_checked=sum(r["samples_checked"] for r in replays),
                          all_fixed_fresh_seeds_included=True, all_pair_initial_physical_states_identical=True,
                          no_R_preparations_retained=no_r, no_L_preparations_retained=no_l,
                          formation_histories=formation, preparation_costs_and_initial_resources_verified=True,
                          all_primary_secondary_readouts_and_group_means_recomputed=True,
                          replication_estimate=calculated, original_confirmation_t_sensitivity=old_interval,
                          previous_decision_unchanged=True, no_pooled_estimate=True, replays=replays,
                          scope="A fixed independent transport comparison at the same supplied drive and parameters. This is a separate 256-pair estimate, not a pooled rescue or a proof of individuality. Full costs, failures and zero-L local preparations are retained. Heat is not total entropy production; supply conditions do not enforce equal realized chemical work."), indent=2))


if __name__ == "__main__":
    main()
