"""Independent fixed-site causal assay audit, without running simulations.

Replays complete saved histories through the frozen scalar reference, rebuilds
selection/masks from the prepared physical endpoint, and calls the separately
audited independent cohort reconstruction. No experiment-runner imports.
"""
from __future__ import annotations

import argparse
from collections import Counter
import copy
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import numpy as np
from verify_differentiated import (ROOT, make_reference, verify_trajectory, load,
                                  digest, compare_state)
from verify_driven import compare_values, close
from verify_study import estimate
from verify_support_readouts import reconstruct_support

MODES, LESIONS = ("opposed", "aligned"), ("sham", "near", "remote")
PRIMARY = ("opposed_sham_minus_near", "opposed_remote_minus_near")
SECONDARY = ("remote_minus_near_response_interaction", "sham_opposed_minus_aligned")
OUTCOME = "outcome_integral_local_R"
FIELDS = (
    OUTCOME, "integral_global_R", "integral_local_L", "integral_global_cohort_R",
    "integral_local_cohort_R", "integral_global_cohort_R1",
    "integral_local_R_after_original_cohort_clearance",
    "integral_cohort_R_after_original_cohort_clearance", "initial_cohort_size",
    "final_cohort_size", "final_cohort_R1", "final_local_R", "final_global_R",
    "cohort_forward_L", "cohort_reverse_L", "cohort_net_L", "near_site_forward_L",
    "near_site_reverse_L", "cohort_copies", "cohort_copies_after_original_clearance",
    "cohort_outward_crossings", "cohort_inward_crossings", "chemical_work",
    "heat_to_bath", "total_preparation_and_continuation_work",
    "total_preparation_and_continuation_heat")


def physical_accounting(run):
    ref = make_reference(run)
    initial, final = run["initial_state"], run["final_state"]
    b0, b1 = initial["base_state"], final["base_state"]
    h0, h1 = ref.energy(initial), ref.energy(final)
    work = math.fsum(event["external_work"] for event in run["ledger"])
    heat = math.fsum(event["heat_to_bath"] for event in run["ledger"])
    first_law = abs(h1-h0-work+heat)
    assert first_law < 3e-7
    close(b1["external_work"]-b0["external_work"], work, 3e-7)
    close(b1["cumulative_heat"]-b0["cumulative_heat"], heat, 3e-7)
    counts = Counter(b0["material"])
    gates = sum(event["kind"] == "gate_flip" for event in run["ledger"])
    return dict(initial_stored_energy_H=h0, final_stored_energy_H=h1,
                initial_material_units=ref.base.n, initial_fuel=sum(b0["fuel"]), final_fuel=sum(b1["fuel"]),
                initial_species_counts={name: counts[species] for name, species in (("D", 0), ("R0", 1), ("R1", 2), ("L", 3))},
                initial_gate=initial["gate"], final_gate=final["gate"], chemical_work=work,
                heat_to_bath=heat, first_law_error=first_law,
                material_reservoir_transport_events=len(run["ledger"])-gates,
                gate_events=gates, ledger_events=len(run["ledger"]))


def expected_gate_path(run):
    gate, seed = run["initial_state"]["gate"], run["seed"]
    rng = np.random.default_rng(np.random.SeedSequence([seed, 62412]))
    clock, result = 0., []
    while True:
        clock += float(rng.exponential(1. / run["response"]["flip_rate"]))
        if clock > run["horizon"]:
            return result
        result.append((clock, gate, -gate))
        gate = -gate


def check_gate_clock(run, stationary_initial=False):
    path = [(e["time"], e["before"], e["after"]) for e in run["ledger"] if e["kind"] == "gate_flip"]
    assert path == expected_gate_path(run)
    if stationary_initial:
        rng = np.random.default_rng(np.random.SeedSequence([run["seed"], 62411]))
        assert run["initial_state"]["gate"] == int(rng.choice((-1, 1)))
    return path


def independent_plan(preparation, min_distance=3.):
    """Select squares using independent scalar physical propensities only."""
    ref, state = make_reference(preparation), preparation["final_state"]
    side, n = preparation["params"]["side"], ref.base.n
    material, mask = state["base_state"]["material"], preparation["production_site_mobility"]
    forward_terms, reverse_terms = [[] for _ in range(n)], [[] for _ in range(n)]
    for index, event in enumerate(ref.base.events):
        if event[0] == 2:
            (forward_terms if event[3] == 0 else reverse_terms)[event[2]].append(ref.rate(state, index))
    forward = [math.fsum(terms) for terms in forward_terms]
    reverse = [math.fsum(terms) for terms in reverse_terms]
    squares = []
    for row in range(side):
        for col in range(side):
            sites = sorted({((row+dr) % side)*side + (col+dc) % side for dr in (0, 1) for dc in (0, 1)})
            counts = Counter(material[i] for i in sites)
            f, r = math.fsum(forward[i] for i in sites), math.fsum(reverse[i] for i in sites)
            squares.append(dict(top_left=[row, col], sites=sites, center=[row+.5, col+.5],
                                R=counts[1]+counts[2], R0=counts[1], R1=counts[2], L=counts[3], D=counts[0],
                                raw_paid_L_forward_propensity=f, raw_paid_L_reverse_propensity=r,
                                raw_paid_L_total_propensity=math.fsum((f, r)), raw_paid_L_net_current=f-r))
    near = min(squares, key=lambda square: (-square["R"], -square["L"], square["top_left"]))

    def distance(first, second):
        return float(sum(min(abs(a-b), side-abs(a-b)) for a, b in zip(first["center"], second["center"])))

    eligible = [square for square in squares if not set(near["sites"]) & set(square["sites"])
                and distance(near, square) >= min_distance]
    hn = near["raw_paid_L_total_propensity"]
    remote = min(eligible, key=lambda square: (abs(square["raw_paid_L_total_propensity"]-hn),
                                              -distance(near, square), square["top_left"])) if eligible else None
    hr = remote["raw_paid_L_total_propensity"] if remote is not None else 0.
    reasons = []
    if near["R"] == 0:
        reasons.append("near_square_has_no_templates")
    if remote is None:
        reasons.append("no_disjoint_remote_square_at_declared_distance")
    if hn == 0:
        reasons.append("zero_near_initial_paid_L_propensity")
    if hr == 0:
        reasons.append("zero_remote_initial_paid_L_propensity")
    matched = not reasons
    qn, qr = (min(1., hr/hn), min(1., hn/hr)) if matched else (0., 0.)
    masks = {"sham": mask[:], "near": mask[:], "remote": mask[:]}
    for i in near["sites"]:
        masks["near"][i] *= 1.-qn
    if remote is not None:
        for i in remote["sites"]:
            masks["remote"][i] *= 1.-qr

    def removed(square, q):
        f = 0. if square is None else q*square["raw_paid_L_forward_propensity"]
        r = 0. if square is None else q*square["raw_paid_L_reverse_propensity"]
        return dict(forward=f, reverse=r, total=f+r, net_current=f-r)

    rn, rr = removed(near, qn), removed(remote, qr)
    near_targets = {i for j in near["sites"] for i in ref.base.neighbors[j]}
    remote_targets = {i for j in remote["sites"] for i in ref.base.neighbors[j]} if remote else set()
    metadata = dict(preparation_time=preparation["horizon"], gate=state["gate"],
        initial_stored_energy_H=ref.energy(state), initial_material_units=n,
        initial_fuel=sum(state["base_state"]["fuel"]), params=preparation["params"],
        flags=preparation["flags"], response=preparation["response"], response_mode=preparation["response_mode"],
        transport_mode=preparation["transport_mode"], near=near, remote=remote,
        minimum_periodic_center_distance=min_distance, periodic_center_distance_metric="Manhattan (L1)",
        actual_periodic_center_distance=distance(near, remote) if remote else None,
        eligible_remote_squares=len(eligible), selection_hazard_difference=abs(hr-hn) if remote else None,
        matched=matched, unmatched_reasons=reasons, near_removed_mobility_fraction=qn,
        remote_removed_mobility_fraction=qr, near_removed_initial_propensity=rn,
        remote_removed_initial_propensity=rr, removed_total_matching_error=rn["total"]-rr["total"],
        possible_L_target_neighborhood_overlap=sorted(near_targets & remote_targets),
        site_paid_L_forward_propensities=forward, site_paid_L_reverse_propensities=reverse,
        state_energy_jump=0., intervention_apparatus_work_modeled=False,
        helper_source_sha256=digest(ROOT / "prepare_lesion.py"))
    assert abs(rn["total"]-rr["total"]) < 1e-10
    return dict(baseline_production_site_mobility=masks["sham"], near_production_site_mobility=masks["near"],
                remote_production_site_mobility=masks["remote"], metadata=metadata), masks


def check_initial_hazard_removal(preparation, plan, masks):
    ref = make_reference(preparation)
    state = preparation["final_state"]
    baseline = [ref.rate(state, i) for i in range(len(ref.events))]
    removed = {}
    for lesion in LESIONS:
        altered = copy.deepcopy(preparation)
        altered["production_site_mobility"] = masks[lesion]
        changed = make_reference(altered)
        f_terms, r_terms = [], []
        for i, event in enumerate(ref.events):
            if i < ref.count and event[0] == 2:
                (f_terms if event[3] == 0 else r_terms).append(baseline[i]-changed.rate(state, i))
            else:
                close(baseline[i], changed.rate(state, i), 1e-14)
        f, r = math.fsum(f_terms), math.fsum(r_terms)
        removed[lesion] = dict(forward=f, reverse=r, total=f+r, net_current=f-r)
        if lesion == "sham":
            close(f+r, 0., 1e-14)
        else:
            compare_values(removed[lesion], plan["metadata"][lesion+"_removed_initial_propensity"])
        close(ref.energy(state), changed.energy(state))
    close(removed["near"]["total"], removed["remote"]["total"], 1e-10)
    return removed


def common_run_checks(run, protocol, horizon, mode, seed):
    compare_values(protocol["params"], run["params"])
    assert run["flags"] == protocol["flags"]
    compare_values(protocol["drive"], run["drive"])
    assert run["response"] == protocol["response"]
    assert run["response_mode"] == mode and run["transport_mode"] == "local"
    assert run["horizon"] == horizon and run["seed"] == seed
    initial = run["initial_state"]
    assert initial["gate"] in (-1, 1) and initial["gate_event_count"] == 0
    b = initial["base_state"]
    assert b["time"] == b["cumulative_heat"] == b["external_work"] == b["event_count"] == 0
    assert [s["time"] for s in run["samples"]] == list(np.r_[np.arange(0., horizon, protocol["sample_dt"]), horizon])


def paired_statistics(rows, protocol, summary):
    groups, by_seed = [], {}
    for row in rows:
        arms = by_seed.setdefault(row["origin_seed"], {})
        key = (row["response_mode"], row["lesion"])
        assert key not in arms
        arms[key] = row
    assert sorted(by_seed) == protocol["seeds"]
    for mode in MODES:
        for lesion in LESIONS:
            arm = [by_seed[seed][(mode, lesion)] for seed in protocol["seeds"]]
            groups.append(dict(response_mode=mode, lesion=lesion,
                statistics={field: estimate([row[field] for row in arm]) for field in FIELDS},
                original_cohort_present=sum(row["original_cohort_present"] for row in arm),
                original_cohort_cleared_after_presence=sum(row["original_cohort_cleared_after_presence"] for row in arm),
                final_zero_global_R=sum(row["final_global_R"] == 0 for row in arm)))
    compare_values(groups, summary["groups"])
    differences = {name: [] for name in PRIMARY+SECONDARY}
    seed_rows = []
    for seed in protocol["seeds"]:
        arms = by_seed[seed]
        y = lambda mode, lesion: arms[(mode, lesion)][OUTCOME]
        sham_near = y("opposed", "sham")-y("opposed", "near")
        remote_near = y("opposed", "remote")-y("opposed", "near")
        values = (sham_near, remote_near,
                  remote_near-y("aligned", "remote")+y("aligned", "near"),
                  y("opposed", "sham")-y("aligned", "sham"))
        seed_rows.append({"origin_seed": seed, **dict(zip(PRIMARY+SECONDARY, values))})
        for name, value in zip(PRIMARY+SECONDARY, values):
            differences[name].append(value)
    compare_values(seed_rows, summary["paired_seed_differences"])
    contrasts = [dict(contrast=name, field=OUTCOME, primary=name in PRIMARY, **estimate(values))
                 for name, values in differences.items()]
    compare_values(contrasts, summary["paired_contrasts"])
    primary = contrasts[:2]
    compare_values(primary, summary["primary_results"])
    joint = all(r["t_95_interval"][0] > 0 for r in primary)
    assert joint == summary["declared_joint_claim_supported"]
    spec = protocol["primary_bootstrap_sensitivity"]
    assert spec["replicates"] == 10000 and spec["seed"] == 1080002
    data = np.column_stack([differences[name] for name in PRIMARY])
    indices = np.random.default_rng(spec["seed"]).integers(0, len(data), size=(10000, len(data)))
    quantiles = np.percentile(data[indices].mean(axis=1), [2.5, 97.5], axis=0)
    bootstrap = [dict(contrast=name, replicates=10000, seed=1080002, percentile_95_interval=quantiles[:, i].tolist(),
                      sensitivity_only=True) for i, name in enumerate(PRIMARY)]
    compare_values(bootstrap, summary["primary_bootstrap_sensitivity"])
    return primary, contrasts[2:], bootstrap


def verify_study(path):
    summary = load(path)
    protocol_path = ROOT / summary["protocol_file"]
    protocol = load(protocol_path)
    assert digest(protocol_path) == summary["protocol_sha256"]
    assert summary["source_hashes"] == protocol["source_hashes"]
    for file, expected in protocol["source_hashes"].items():
        assert digest(ROOT / file) == expected, file
    assert digest(ROOT / protocol["selected_parameter_protocol_file"]) == protocol["selected_parameter_protocol_sha256"]
    assert summary["selected_parameter_protocol_sha256"] == protocol["selected_parameter_protocol_sha256"]
    assert protocol["params"]["mutation"] == [[.95, .05], [.05, .95]]
    assert protocol["response"] == dict(asymmetry=.8, flip_rate=.025, initial_gate=None)
    assert protocol["seeds"] == list(range(1070001, 1070065))
    assert protocol["preparation_horizon"] == 100. and protocol["continuation_horizon"] == 80.
    assert protocol["post_seed_offset"] == 3000000 and protocol["sample_dt"] == 2.
    assert protocol["preparation_response_mode"] == "opposed" and protocol["transport_mode"] == "local"
    assert protocol["response_modes"] == list(MODES) and protocol["lesions"] == list(LESIONS)
    assert protocol["flags"] == dict(produce_L=True, protect_R=True, copy_R=True)
    assert protocol["primary_contrasts"] == list(PRIMARY) and protocol["named_secondary_contrasts"] == list(SECONDARY)
    assert protocol["primary_outcome"] == OUTCOME
    assert protocol["expected_preparations"] == 64 and protocol["expected_continuations"] == 384
    assert len(summary["preparations"]) == 64 and len(summary["runs"]) == 384
    assert summary["all_preparations_included"] is True
    prepared, files, preparation_audits, branch_audits = {}, set(), [], []
    for row in summary["preparations"]:
        seed = row["origin_seed"]
        assert seed in protocol["seeds"] and seed not in prepared
        file = ROOT / row["trajectory_file"]
        assert file not in files and digest(file) == row["trajectory_sha256"]
        files.add(file)
        run = load(file)
        common_run_checks(run, protocol, 100., "opposed", seed)
        assert run["production_site_mobility"] == [1.] * protocol["params"]["side"]**2
        initial = run["initial_state"]["base_state"]
        assert set(initial["material"]) <= {0, 1}
        assert initial["material"].count(1) == protocol["params"]["initial_templates"]
        assert sum(initial["fuel"]) == protocol["params"]["initial_fuel"]
        check_gate_clock(run, stationary_initial=True)
        replay = verify_trajectory(run)
        physics = physical_accounting(run)
        compare_values(physics, row["physical_accounting"])
        plan, masks = independent_plan(run, protocol["selection"]["minimum_periodic_manhattan_center_distance"])
        compare_values(plan, row["intervention_plan"])
        hazards = check_initial_hazard_removal(run, plan, masks)
        expected_initial = copy.deepcopy(run["final_state"])
        expected_initial["base_state"].update(time=0., cumulative_heat=0., external_work=0., event_count=0)
        expected_initial["gate_event_count"] = 0
        assert row["branch_initial_state"] == expected_initial
        zero_r = not any(s in (1, 2) for s in run["final_state"]["base_state"]["material"])
        assert zero_r == row["final_zero_R"]
        prepared[seed] = dict(row=row, physics=physics, initial=expected_initial, plan=plan, masks=masks,
                              paths={}, gate_path=None)
        preparation_audits.append(dict(trajectory_file=row["trajectory_file"], trajectory_sha256=row["trajectory_sha256"],
            removed_initial_hazards=hazards, matched=plan["metadata"]["matched"], **replay))
    assert sorted(prepared) == protocol["seeds"]
    independent_rows, seen = [], set()
    for number, row in enumerate(summary["runs"]):
        seed, mode, lesion = row["origin_seed"], row["response_mode"], row["lesion"]
        key = (seed, mode, lesion)
        assert seed in prepared and mode in MODES and lesion in LESIONS and key not in seen
        seen.add(key)
        prefix = prepared[seed]
        assert row["post_seed"] == seed+protocol["post_seed_offset"]
        file = ROOT / row["trajectory_file"]
        assert file not in files and digest(file) == row["trajectory_sha256"]
        files.add(file)
        assert row["preparation_trajectory_file"] == prefix["row"]["trajectory_file"]
        assert row["preparation_trajectory_sha256"] == prefix["row"]["trajectory_sha256"]
        run = load(file)
        common_run_checks(run, protocol, 80., mode, row["post_seed"])
        assert run["initial_state"] == prefix["initial"]
        assert run["production_site_mobility"] == prefix["masks"][lesion]
        assert run["near_sites"] == prefix["plan"]["metadata"]["near"]["sites"]
        compare_values(prefix["plan"]["metadata"], run["lesion_metadata"])
        compare_values(dict(origin_seed=seed, lesion=lesion, preparation_trajectory_file=prefix["row"]["trajectory_file"],
            preparation_trajectory_sha256=prefix["row"]["trajectory_sha256"], protocol_file=protocol_path.name,
            protocol_sha256=digest(protocol_path), observer_annotations_only=True), run["causal_support_context"])
        gate_path = check_gate_clock(run)
        if prefix["gate_path"] is None:
            prefix["gate_path"] = gate_path
        else:
            assert prefix["gate_path"] == gate_path
        if not prefix["plan"]["metadata"]["matched"]:
            identity = (run["ledger"], run["final_state"], run["samples"])
            if mode in prefix["paths"]:
                assert prefix["paths"][mode] == identity, "No-op lesions must have identical paths within a response mode"
            else:
                prefix["paths"][mode] = identity
        replay = verify_trajectory(run)
        physics = physical_accounting(run)
        compare_values(physics, row)
        reconstruction = reconstruct_support(run, run["near_sites"])
        measurements = reconstruction["measurements"]
        compare_values(measurements, row)
        assert row["matched"] == prefix["plan"]["metadata"]["matched"]
        pw, pq = prefix["physics"]["chemical_work"], prefix["physics"]["heat_to_bath"]
        total_work, total_heat = pw+physics["chemical_work"], pq+physics["heat_to_bath"]
        original_h = prefix["physics"]["initial_stored_energy_H"]
        error = abs(physics["final_stored_energy_H"]-original_h-total_work+total_heat)
        assert error < 5e-7
        totals = dict(preparation_chemical_work=pw, preparation_heat_to_bath=pq,
            original_initial_stored_energy_H=original_h, total_preparation_and_continuation_work=total_work,
            total_preparation_and_continuation_heat=total_heat, total_preparation_and_continuation_first_law_error=error)
        compare_values(totals, row)
        independent_rows.append({**row, **measurements, **physics, **totals})
        branch_audits.append(dict(trajectory_file=row["trajectory_file"], trajectory_sha256=row["trajectory_sha256"],
            origin_seed=seed, lesion=lesion, cohort_diagnostics=reconstruction["diagnostics"], **replay))
        if (number+1) % 32 == 0:
            print(f"Independently replayed 64 preparations and {number+1}/384 causal branches", flush=True)
    assert seen == {(seed, mode, lesion) for seed in protocol["seeds"] for mode in MODES for lesion in LESIONS}
    primary, secondary, bootstrap = paired_statistics(independent_rows, protocol, summary)
    matching = dict(matched=sum(p["plan"]["metadata"]["matched"] for p in prepared.values()),
        unmatched=sum(not p["plan"]["metadata"]["matched"] for p in prepared.values()),
        final_zero_R_preparations=sum(p["row"]["final_zero_R"] for p in prepared.values()),
        maximum_removed_total_matching_error=max(abs(p["plan"]["metadata"]["removed_total_matching_error"]) for p in prepared.values()))
    compare_values(matching, summary["matching_summary"])
    unique_events = sum(r["events_replayed"] for r in preparation_audits+branch_audits)
    totals = dict(unique_preparation_records=64, continuation_records=384, unique_ledger_events=unique_events,
        preparation_work_counted_once=sum(p["physics"]["chemical_work"] for p in prepared.values()),
        continuation_work_all_counterfactual_arms=sum(r["chemical_work"] for r in independent_rows),
        maximum_trajectory_first_law_error=max([p["physics"]["first_law_error"] for p in prepared.values()]+[r["first_law_error"] for r in independent_rows]),
        maximum_preparation_plus_branch_first_law_error=max(r["total_preparation_and_continuation_first_law_error"] for r in independent_rows))
    compare_values(totals, summary["physical_totals"])
    assert len(files) == 448
    for file, expected in protocol["source_hashes"].items():
        assert digest(ROOT / file) == expected, file
    return dict(checks_passed=True, summary_file=path.name, summary_sha256=digest(path),
        protocol_file=protocol_path.name, protocol_sha256=digest(protocol_path),
        source_hashes=protocol["source_hashes"], preparations=64, continuations=384,
        unique_physical_records=448, unique_physical_events=unique_events,
        preparation_events=sum(r["events_replayed"] for r in preparation_audits),
        continuation_events=sum(r["events_replayed"] for r in branch_audits),
        preparation_samples=sum(r["samples_checked"] for r in preparation_audits),
        continuation_samples=sum(r["samples_checked"] for r in branch_audits),
        all_seeds_and_six_cells_retained=True, prefix_physical_state_and_actual_gate_preserved=True,
        six_shared_post_gate_paths_independently_regenerated=True,
        static_helper_masks_and_initial_total_hazard_matching_verified=True,
        no_op_paths_identical_within_response_modes=True,
        matching_summary=matching, primary_results=primary, secondary_results=secondary,
        primary_bootstrap_sensitivity=bootstrap, declared_joint_claim_supported=summary["declared_joint_claim_supported"],
        physical_totals=totals, preparation_audits=preparation_audits, continuation_audits=branch_audits,
        scope="Fixed external helper-site lesions on opposed-prepared physical states. Both directions of the paid-L channel are changed. Initial total propensity matching is not matching net current, future flux or apparatus cost. Region/cohort readouts do not establish ownership, self-produced boundaries, selection or collective individuality. Prefix histories count once in data totals; alternative branches are not simultaneous physical copies.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("summary", nargs="?", default="causal-support-results.json")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = verify_study(ROOT / args.summary)
    result.update(verifier_sha256=digest(__file__),
        independent_physics_verifier_sha256=digest(ROOT / "verify_differentiated.py"),
        independent_cohort_verifier_sha256=digest(ROOT / "verify_support_readouts.py"),
        independent_statistics_verifier_sha256=digest(ROOT / "verify_study.py"))
    Path(args.output).write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps({key: result[key] for key in ("checks_passed", "unique_physical_records", "unique_physical_events",
        "preparation_samples", "continuation_samples", "matching_summary", "primary_results", "declared_joint_claim_supported")}, indent=2))


if __name__ == "__main__":
    main()
