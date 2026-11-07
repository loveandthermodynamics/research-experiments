"""Independent full active-phase study replay and seed-level decision audit.

Imports only the frozen independent scalar physical/readout reference. Neither
production analysis nor the study runner is used to recompute measurements.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
import statistics
import sys

sys.dont_write_bytecode = True
import numpy as np
from scipy.stats import t
from verify_active_phase import (ROOT, Reference, compare, close, digest, load,
                                 observe_independent, verify_trajectory)

ARMS = ("full", "no_drive", "no_affinity")
METRICS = ("B_fraction", "largest_B_fraction", "mass_weighted_cluster_size",
           "excess_neighbor_B_fraction", "E_on_B_excess_fraction", "E_enrichment",
           "spectral_wavelength_low_k", "low_k_power_fraction", "clusters_ge_9")


def average(values):
    return None if any(value is None for value in values) else statistics.mean(values)


def ratio(numerator, denominator):
    return None if numerator is None or denominator is None or denominator == 0 else numerator/denominator


def difference(first, second):
    return None if first is None or second is None else first-second


def estimate(values):
    if any(value is None or not math.isfinite(value) for value in values):
        return dict(estimable=False, n=len(values), values=values,
                    reason="undefined required values retained; no seed filtering")
    n, mean = len(values), statistics.mean(values)
    sem = statistics.stdev(values)/math.sqrt(n)
    margin = float(t.ppf(.975, n-1))*sem
    indices = np.random.default_rng(384701).integers(n, size=(10000, n))
    bootstrap = np.asarray(values)[indices].mean(axis=1)
    return dict(estimable=True, n=n, values=values, mean=mean, sem=sem,
                t95=[mean-margin, mean+margin], paired_seed_bootstrap95=np.percentile(bootstrap, [2.5, 97.5]).tolist())


def gate(name, values, relation, threshold):
    out = estimate(values)
    out.update(name=name, relation=relation, threshold=threshold)
    if not out["estimable"]:
        out["met"] = False
    else:
        out["met"] = out["t95"][0] > threshold if relation == "lower_gt" else out["t95"][1] <= threshold
    return out


def initial_state(seed, params):
    n = params["side"]**2
    material, enzyme, tags = [0]*n, [0]*n, [0]*n
    selected = np.random.default_rng(np.random.SeedSequence([seed, 74501])).choice(n, params["material_count"], replace=False)
    for tag, site in enumerate(selected, 1):
        material[int(site)], tags[int(site)] = 1, tag
    selected = np.random.default_rng(np.random.SeedSequence([seed, 74502])).choice(n, params["enzyme_count"], replace=False)
    for site in selected:
        enzyme[int(site)] = 1
    return dict(material=material, enzyme=enzyme, material_tags=tags, time=0., cumulative_work=0.,
                cumulative_heat=0., net_fuel_to_waste=0, event_count=0)


def measurements(run, path):
    rows = []
    for sample in run["samples"]:
        row = observe_independent(sample["material"], sample["enzyme"])
        row["large_B_mass_fraction"] = row["B_mass_in_clusters_ge_9"]/row["n_B"] if row["n_B"] else 0.
        row["time"] = sample["time"]
        rows.append(row)
    windows = {}
    for name, bounds in (("early", (100., 200.)), ("late", (300., 400.))):
        selected = [row for row in rows if bounds[0] <= row["time"] <= bounds[1]]
        assert len(selected) == 21
        windows[name] = {metric: average([row[metric] for row in selected]) for metric in (*METRICS, "large_B_mass_fraction")}
    snapshot = next(sample for sample in run["samples"] if sample["time"] == 200.)
    cohort = {tag for a, tag in zip(snapshot["material"], snapshot["material_tags"]) if a == 2}
    # Direct tag -> complete conversion history, independent of live membership.
    histories = {tag: [] for tag in cohort}
    conversions, late = Counter(), Counter()
    for event in run["ledger"]:
        if event["kind"] not in ("passive", "catalyzed"):
            continue
        before, after = event["before"]["material"][0], event["after"]["material"][0]
        label = f"{event['kind']}_{before}_to_{after}"
        if event["time"] > 200.:
            conversions[label] += 1
            tag = event["before"]["material_tags"][0]
            if tag in histories:
                histories[tag].append((event["time"], before, after))
        if event["time"] > 300.:
            late[label] += 1
    converted = {tag for tag, events in histories.items() if any(before == 2 for _, before, after in events)}
    final_b = {tag for a, tag in zip(run["final_state"]["material"], run["final_state"]["material_tags"]) if a == 2}
    turnover = dict(B_carrier_cohort_at200=len(cohort), cohort_tags_converted_to_A_after200=len(converted),
                    cohort_tags_B_at400_anywhere=len(cohort & final_b),
                    cohort_tags_never_converted_to_A_after200=len(cohort-converted),
                    conversions_after200=dict(conversions), conversions_after300=dict(late))
    return dict(side=run["params"]["side"], arm=run["arm"], seed=run["seed"], path=str(path.relative_to(ROOT)),
                record_sha256=digest(path), windows=windows,
                growth=ratio(windows["late"]["spectral_wavelength_low_k"], windows["early"]["spectral_wavelength_low_k"]),
                turnover=turnover, timeseries=rows, physical_final=run["final_summary"])


def verify_analysis(protocol, execution, analysis, independent_rows):
    assert analysis["protocol_sha256"] == digest(ROOT / "protocol.json")
    assert analysis["analysis_sha256"] == digest(ROOT / "analyze_study.py")
    assert analysis["unique_records"] == 48 and analysis["independent_seed_rows"] == 8
    observed = {(row["side"], row["arm"], row["seed"]): row for row in analysis["rows"]}
    index = {}
    for row in independent_rows:
        key = (row["side"], row["arm"], row["seed"])
        assert key not in index
        index[key] = row
        compare(row, observed[key])
    seeds = protocol["seeds"]
    get = lambda side, arm, seed: index[(side, arm, seed)]
    gates = []
    for side in (16, 32):
        for metric, threshold in (("B_fraction", .05), ("large_B_mass_fraction", .5), ("excess_neighbor_B_fraction", .15)):
            gates.append(gate(f"L{side}_late_{metric}", [get(side, "full", seed)["windows"]["late"][metric] for seed in seeds], "lower_gt", threshold))
        gates.append(gate(f"L{side}_full_wavelength_growth", [get(side, "full", seed)["growth"] for seed in seeds], "upper_le", 1.25))
        gates.append(gate(f"L{side}_no_drive_minus_full_growth", [difference(get(side, "no_drive", seed)["growth"], get(side, "full", seed)["growth"]) for seed in seeds], "lower_gt", 0.))
    gates.append(gate("full_late_wavelength_L32_over_L16", [ratio(get(32, "full", seed)["windows"]["late"]["spectral_wavelength_low_k"], get(16, "full", seed)["windows"]["late"]["spectral_wavelength_low_k"]) for seed in seeds], "upper_le", 1.5))
    gates.append(gate("full_late_largest_B_fraction_L16_minus_L32", [difference(get(16, "full", seed)["windows"]["late"]["largest_B_fraction"], get(32, "full", seed)["windows"]["late"]["largest_B_fraction"]) for seed in seeds], "lower_gt", 0.))
    compare(gates, analysis["primary_gates"])
    assert analysis["joint_claim_met"] == all(g["met"] for g in gates)
    summaries = {}
    for side in (16, 32):
        for arm in ARMS:
            selected = [get(side, arm, seed) for seed in seeds]
            summaries[f"L{side}-{arm}"] = dict(
                late={metric: estimate([row["windows"]["late"][metric] for row in selected]) for metric in (*METRICS, "large_B_mass_fraction")},
                growth=estimate([row["growth"] for row in selected]),
                chemical_work_per_area=estimate([row["physical_final"]["chemical_work"]/side**2 for row in selected]),
                heat_to_bath_per_area=estimate([row["physical_final"]["heat_to_bath"]/side**2 for row in selected]),
                turnover=[row["turnover"] for row in selected])
    compare(summaries, analysis["summary"])
    return gates, summaries


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--accepted-only", action="store_true", help="Explicitly omit full rejected-proposal regeneration; accepted physical ledgers still replay completely")
    args = parser.parse_args()
    protocol = load(ROOT / "protocol.json")
    protocol_sha = digest(ROOT / "protocol.json")
    assert protocol["seeds"] == list(range(720100, 720108))
    assert protocol["sides"] == [16, 32] and protocol["arms"] == list(ARMS)
    assert protocol["horizon"] == 400. and protocol["sample_dt"] == 5.
    assert protocol["early_window"] == [100., 200.] and protocol["late_window"] == [300., 400.]
    assert protocol["no_outcome_parameter_selection"] and protocol["no_rescue_or_optional_extension"]
    for name, expected in protocol["source_sha256"].items():
        assert digest(ROOT / name) == expected, name
    assert digest(ROOT / "preflight-verification.json") == protocol["independent_preflight_sha256"]
    assert load(ROOT / "preflight-verification.json")["checks_passed"]
    execution = load(ROOT / "execution.json")
    assert execution["protocol_sha256"] == protocol_sha
    assert execution["passed_execution"] and execution["all_preparations_retained"] and execution["unique_records"] == 48
    expected_cells = {(side, arm, seed) for side in (16, 32) for arm in ARMS for seed in protocol["seeds"]}
    seen, files, initial_by_pair = set(), set(), {}
    audits, rows = [], []
    for number, record in enumerate(execution["records"]):
        side, arm, seed = record["side"], record["arm"], record["seed"]
        key = (side, arm, seed)
        assert key in expected_cells and key not in seen
        seen.add(key)
        path = ROOT / record["path"]
        assert path not in files and digest(path) == record["sha256"]
        assert path.stat().st_size == record["bytes"]
        files.add(path)
        run = load(path)
        assert run["arm"] == arm and run["seed"] == seed and run["horizon"] == 400.
        assert run["protocol_sha256"] == protocol_sha
        assert run["params"] == protocol["parameters"][f"{side}-{arm}"]
        assert run["params"]["material_count"] == 77*side**2//256
        assert run["params"]["enzyme_count"] == 8*side**2//256
        assert run["initial_state"] == initial_state(seed, run["params"])
        pair = (side, seed)
        if pair in initial_by_pair:
            assert run["initial_state"] == initial_by_pair[pair]
        else:
            initial_by_pair[pair] = run["initial_state"]
        assert [sample["time"] for sample in run["samples"]] == [float(x) for x in range(0, 401, 5)]
        audit = verify_trajectory(run, full_proposals=not args.accepted_only)
        assert record["events"] == audit["accepted_events_replayed"]
        assert record["proposals"] == audit["proposals"] and record["samples"] == audit["samples_checked"]
        compare(record["final"], run["final_summary"])
        close(record["first_law_error"], run["first_law_error"])
        rows.append(measurements(run, path))
        audits.append(dict(path=record["path"], sha256=record["sha256"], side=side, arm=arm, **audit))
        progress = dict(completed_records=len(audits), expected_records=48, protocol_sha256=protocol_sha,
                        exact_all_proposals=not args.accepted_only,
                        events_replayed=sum(r["accepted_events_replayed"] for r in audits), records=audits)
        (ROOT / "independent-verification-progress.json").write_text(json.dumps(progress, indent=2)+"\n")
        print(f"Independently replayed {number+1}/48: L{side} {arm} seed{seed}, {audit['accepted_events_replayed']} events / {audit['proposals']} proposals", flush=True)
    assert seen == expected_cells and len(files) == 48
    assert execution["events"] == sum(r["accepted_events_replayed"] for r in audits)
    assert execution["samples"] == sum(r["samples_checked"] for r in audits)
    analysis = load(ROOT / "analysis.json")
    gates, summaries = verify_analysis(protocol, execution, analysis, rows)
    for name, expected in protocol["source_sha256"].items():
        assert digest(ROOT / name) == expected, name
    report = dict(checks_passed=True, verifier_sha256=digest(__file__),
        independent_physical_reference_sha256=digest(ROOT / "verify_active_phase.py"),
        protocol_sha256=protocol_sha, execution_sha256=digest(ROOT / "execution.json"),
        analysis_sha256=digest(ROOT / "analysis.json"), all48cells_retained=True,
        matched_initial_states_all16seed_size_pairs=True, zero_B_initial_all_records=True,
        exact_density_preservation=True, unique_records=48, independent_seed_rows=8,
        accepted_events=sum(r["accepted_events_replayed"] for r in audits),
        proposals=sum(r["proposals"] for r in audits), samples=sum(r["samples_checked"] for r in audits),
        every_accepted_ledger_replayed=True, every_proposal_regenerated=not args.accepted_only,
        morphology_windows_turnover_and_statistics_match=True,
        joint_claim_met=analysis["joint_claim_met"], primary_gates=gates, summary=summaries, records=audits,
        scope="One bounded48-record study only. Excludes software smokes. Null metrics/seeds retained; Student gates and paired bootstrap use8 seed rows, not81 correlated time snapshots. Size tolerances and low-k/cluster cutoffs are supplied observation rules, never physical feedback. Open supplied chemistry and conserved enzyme are not copied roles, autonomy, recursive individuality, aging, biodiversity maximization or ethical legitimacy.")
    Path(args.output).write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps({key: report[key] for key in ("checks_passed", "unique_records", "accepted_events", "proposals", "samples", "every_proposal_regenerated", "joint_claim_met")}, indent=2))


if __name__ == "__main__":
    main()
