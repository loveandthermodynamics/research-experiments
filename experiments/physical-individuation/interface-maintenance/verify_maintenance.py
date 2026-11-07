"""Audit fixed-time preparation and continuation readouts independently.

Run after the unchanged verify_interface.py has written the new maintenance
verification report. This adds preparation/selection/readout checks to that
report. It never calls the production replay or readout implementations.
"""
from __future__ import annotations

import copy
import gzip
from hashlib import sha256
import json
from pathlib import Path
import statistics
import sys

sys.dont_write_bytecode = True
from verify_interface import Reference, Params, Flags, close, compare_record, compare_state

ROOT = Path(__file__).resolve().parent


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def load_gzip(path):
    with gzip.open(path, "rt") as stream:
        return json.load(stream)


def exposures(run):
    """Integrate population counts with an event-delta implementation."""
    m = run["initial_state"]["material"]
    counts = [m.count(k) for k in (1, 2, 3)]
    area = [0., 0., 0.]
    previous = 0.
    first_zero = 0. if counts[0] + counts[1] == 0 else None
    paid_r = paid_l = reverse_r = reverse_l = deaths_r = deaths_l = 0
    for event in run["ledger"]:
        t = event["time"]
        assert previous <= t <= run["horizon"]
        for j in range(3):
            area[j] += (t - previous) * counts[j]
        previous = t
        if event["kind"] in ("material_swap", "fuel_swap"):
            continue
        before, after = event["before"], event["after"]
        for j, species in enumerate((1, 2, 3)):
            counts[j] += int(after == species) - int(before == species)
        if event["kind"] == "copy_R":
            paid_r += int(after in (1, 2))
            reverse_r += int(before in (1, 2))
        elif event["kind"] == "produce_L":
            paid_l += int(after == 3)
            reverse_l += int(before == 3)
        elif event["kind"] == "passive":
            deaths_r += int(before in (1, 2))
            deaths_l += int(before == 3)
        if first_zero is None and counts[0] + counts[1] == 0:
            first_zero = t
    for j in range(3):
        area[j] += (run["horizon"] - previous) * counts[j]
    final = run["final_state"]
    assert counts == [final["material"].count(k) for k in (1, 2, 3)]
    return dict(integral_R=area[0]+area[1], integral_R0=area[0],
                integral_R1=area[1], integral_L=area[2],
                final_R=counts[0]+counts[1], final_L=counts[2],
                paid_R_births=paid_r, paid_L_births=paid_l,
                reverse_R_copies=reverse_r, reverse_L_production=reverse_l,
                passive_R_deaths=deaths_r, passive_L_deaths=deaths_l,
                first_episode_duration=run["horizon"] if first_zero is None else first_zero,
                net_fuel_consumed=sum(run["initial_state"]["fuel"])-sum(final["fuel"]),
                net_heat_to_bath=final["cumulative_heat"])


def main():
    result = json.loads((ROOT / "maintenance-results.json").read_text())
    protocol = json.loads((ROOT / result["protocol_file"]).read_text())
    assert digest(ROOT / result["protocol_file"]) == result["protocol_sha256"]
    assert digest(ROOT / protocol["source_results"]) == protocol["source_results_sha256"]
    assert result["origin_results_sha256"] == protocol["source_results_sha256"]
    assert result["source_hashes"] == protocol["source_hashes"]
    for name, expected in result["source_hashes"].items():
        assert digest(ROOT / name) == expected
    original_result = json.loads((ROOT / protocol["source_results"]).read_text())
    original_rows = {r["seed"]: r for r in original_result["runs"]
                     if r["flags"] == dict(produce_L=True, protect_R=True, copy_R=True)}
    expected_seeds = sorted(original_result["seeds"])
    assert sorted(original_rows) == expected_seeds
    assert result["origin_seeds"] == expected_seeds
    assert [x["origin_seed"] for x in protocol["preparations"]] == expected_seeds
    prepared_rows = {x["origin_seed"]: x for x in result["preparations"]}
    assert len(prepared_rows) == len(result["preparations"]) == 64
    assert sorted(prepared_rows) == expected_seeds
    assert result["branch_time"] == protocol["branch_time"] == 20.
    assert result["horizon"] == protocol["continuation_horizon"] == 80.
    assert result["replicates_per_arm"] == protocol["replicates_per_arm"] == 64
    expected_arms = {(True, True), (True, False), (False, True), (False, False)}
    by_origin = {seed: {} for seed in expected_seeds}
    for row in result["runs"]:
        seed = row["origin_seed"]
        arm = row["flags"]["produce_L"], row["flags"]["protect_R"]
        assert seed in by_origin and arm not in by_origin[seed]
        assert row["flags"]["copy_R"] is True
        by_origin[seed][arm] = row
    assert len(result["runs"]) == 256
    assert all(set(rows) == expected_arms for rows in by_origin.values())
    prefix_events = 0
    checked_rows = []
    fields_checked = set()
    preparation_sizes = []
    for item in protocol["preparations"]:
        seed = item["origin_seed"]
        metadata = prepared_rows[seed]
        assert item["source_trajectory"] == original_rows[seed]["trajectory_file"]
        source_path = ROOT / item["source_trajectory"]
        assert digest(source_path) == item["source_trajectory_sha256"]
        original = load_gzip(source_path)
        assert original["params"] == protocol["params"] == result["params"]
        assert original["seed"] == seed
        assert original["initial_state"] == metadata["original_initial_state"]
        state = copy.deepcopy(original["initial_state"])
        ref = Reference(Params(**original["params"]), Flags(**original["flags"]))
        prefix = [e for e in original["ledger"] if e["time"] <= 20.]
        assert prefix == original["ledger"][:len(prefix)]
        for saved in prefix:
            state["time"] = saved["time"]
            compare_record(ref.apply(state, saved["event_index"]), saved)
        state["time"] = 20.
        compare_state(state, metadata["prepared_state_before_bookkeeping_reset"])
        assert len(prefix) == metadata["preparation_prefix_events"]
        prefix_events += len(prefix)
        pre_run = dict(initial_state=original["initial_state"], final_state=state,
                       horizon=20., ledger=prefix)
        for field, value in exposures(pre_run).items():
            close(value, metadata["preparation_measurements"][field])
        sizes = dict(prepared_R=state["material"].count(1)+state["material"].count(2),
                     prepared_L=state["material"].count(3), prepared_fuel=sum(state["fuel"]))
        assert all(metadata[k] == v for k, v in sizes.items())
        close(metadata["prepared_H"], ref.energy(state))
        close(metadata["preparation_heat_to_bath"], state["cumulative_heat"])
        close(metadata["preparation_external_work"], state["external_work"])
        preparation_sizes.append(sizes)
        for key in ("time", "cumulative_heat", "external_work", "event_count"):
            state[key] = 0
        for row in by_origin[seed].values():
            path = ROOT / row["trajectory_file"]
            assert digest(path) == row["trajectory_sha256"]
            run = load_gzip(path)
            assert run["origin_seed"] == seed and run["preparation"] == metadata
            assert run["seed"] == row["post_seed"] == seed + protocol["post_seed_offset"]
            assert run["initial_state"] == state
            assert run["params"] == protocol["params"] and run["flags"] == row["flags"]
            assert run["horizon"] == 80.
            assert run["source_hashes"] == result["source_hashes"]
            assert run["maintenance_protocol_sha256"] == result["protocol_sha256"]
            assert all(row[k] == v for k, v in sizes.items())
            measured = exposures(run)
            for field, value in measured.items():
                close(value, row[field])
                fields_checked.add(field)
            checked_rows.append((seed, (row["flags"]["produce_L"], row["flags"]["protect_R"]), measured))
    index = {(seed, arm): values for seed, arm, values in checked_rows}
    for contrast in result["paired_contrasts"]:
        values = [sum(weight["weight"] * index[(seed, (weight["produce_L"], weight["protect_R"]))][contrast["field"]]
                      for weight in contrast["definition"]) for seed in expected_seeds]
        mean = statistics.mean(values)
        se = statistics.stdev(values) / len(values)**.5
        assert contrast["n"] == 64
        close(mean, contrast["mean"])
        close(se, contrast["monte_carlo_standard_error"])
        for actual, expected in zip(contrast["approximate_95_percent_normal_interval"],
                                    (mean - 1.96*se, mean + 1.96*se)):
            close(actual, expected)
    report_path = ROOT / "maintenance-verification.json"
    report = json.loads(report_path.read_text())
    assert report["checks_passed"] and report["trajectory_records_replayed"] == 256
    report["fixed_time_preparation_audit"] = dict(
        checks_passed=True, verifier_sha256=digest(Path(__file__)),
        preparation_records_replayed=64, preparation_prefix_events_replayed=prefix_events,
        continuations_with_identical_prepared_physical_state=256,
        all_original_full_arm_seeds_included=True, state_dependent_filtering=False,
        branch_time=20., continuation_horizon=80., physical_state_reset=False,
        bookkeeping_reset_only=["time", "cumulative_heat", "external_work", "event_count"],
        initial_R_range=[min(x["prepared_R"] for x in preparation_sizes), max(x["prepared_R"] for x in preparation_sizes)],
        initial_L_range=[min(x["prepared_L"] for x in preparation_sizes), max(x["prepared_L"] for x in preparation_sizes)],
        readouts_recomputed_independently=sorted(fields_checked),
        paired_contrasts_recomputed=len(result["paired_contrasts"]),
        all_source_protocol_trajectory_hashes_match=True,
        scope="Post-formation ongoing production versus legacy protection; no claim of closed-loop benefit, membrane, aging or collective individuality follows from verification alone.")
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["fixed_time_preparation_audit"], indent=2))


if __name__ == "__main__":
    main()
