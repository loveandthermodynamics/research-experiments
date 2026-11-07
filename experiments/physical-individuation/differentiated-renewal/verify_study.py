"""Independent saved-study audit; lifetime-graph ancestry, not live-set reader.

This imports the independent scalar physics replay, not the experiment runner
or its measurement implementation. All primary differences are reconstructed
from raw ledgers, with a separate statistics implementation.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
import statistics
import sys

sys.dont_write_bytecode = True
import numpy as np
from scipy.stats import t
from verify_differentiated import make_reference, verify_trajectory, digest, load, ROOT
from verify_driven import compare_values, close

PRIMARY = ("response_effect_full_renewal", "response_by_production_interaction")
FIELDS = ("integral_R", "integral_R1", "integral_L", "both_types_present_time",
          "integral_founder_descendant_R1", "integral_focal_descendants", "integral_focal_R1",
          "focal_paid_mutant_appeared_by_deadline", "focal_established", "focal_survives_to_horizon",
          "focal_R1_uninterrupted_to_horizon", "final_focal_descendants", "final_focal_R1",
          "paid_R_births_after_founder_instance_clearance", "chemical_work", "heat_to_bath",
          "final_R", "final_R1", "final_L", "first_episode_duration")


def graph_readouts(run, deadline=20., threshold=4):
    """First build the complete parent graph and instance lifetimes, then integrate."""
    initial, final = run["initial_state"]["base_state"], run["final_state"]["base_state"]
    horizon = run["horizon"]
    nodes = {}
    for tag, species in zip(initial["instance_tags"], initial["material"]):
        if tag:
            nodes[tag] = dict(species=species, born=0., died=horizon, removed=False,
                              parent=0, channel="initial", helper=0)
    founders = {tag for tag, node in nodes.items() if node["species"] in (1, 2)}
    candidates = []
    event_counts = Counter()
    gate_times = Counter()
    gate, previous_gate_time = run["initial_state"]["gate"], 0.
    for event in run["ledger"]:
        event_counts[event["kind"]] += 1
        if event["kind"] == "gate_flip":
            gate_times[str(gate)] += event["time"] - previous_gate_time
            gate, previous_gate_time = event["after"], event["time"]
        if event["instance_removed"]:
            tag = event["instance_removed"]
            assert tag in nodes and not nodes[tag]["removed"]
            nodes[tag].update(died=event["time"], removed=True)
        if event["instance_created"]:
            tag = event["instance_created"]
            assert tag not in nodes
            parent = event["parent_instance"] if event["kind"] == "copy_R" else 0
            nodes[tag] = dict(species=event["after"], born=event["time"], died=horizon,
                              removed=False, parent=parent, channel=event["kind"],
                              helper=event["helper_species"])
            if event["kind"] == "copy_R":
                assert parent in nodes and nodes[parent]["species"] == event["helper_species"]
                assert nodes[parent]["born"] <= event["time"] <= nodes[parent]["died"]
                if event["helper_species"] == 1 and event["after"] == 2 and event["time"] <= deadline:
                    candidates.append(tag)
    gate_times[str(gate)] += horizon - previous_gate_time
    focal = candidates[0] if candidates else None
    founder_descendants, focal_descendants = set(founders), set()
    # Tags are strictly increasing on creation; graph traversal independently
    # checks all genealogy instead of updating the production reader's live set.
    for tag in sorted(nodes):
        parent = nodes[tag]["parent"]
        if parent in founder_descendants:
            founder_descendants.add(tag)
        if tag == focal or parent in focal_descendants:
            focal_descendants.add(tag)
    lengths = {tag: node["died"] - node["born"] for tag, node in nodes.items()}
    total = lambda tags: sum(lengths[tag] for tag in tags)
    by_species = {s: {tag for tag, node in nodes.items() if node["species"] == s} for s in (1, 2, 3)}
    final_tags = {tag for tag, node in nodes.items() if not node["removed"]}
    expected_final = {int(tag): int(species) for tag, species in zip(final["instance_tags"], final["material"]) if tag}
    assert expected_final == {tag: nodes[tag]["species"] for tag in final_tags}
    final_counts = Counter(nodes[tag]["species"] for tag in final_tags)
    by_time = defaultdict(list)
    for tag, node in nodes.items():
        if node["channel"] != "initial":
            by_time[node["born"]].append((1, tag))
        if node["removed"]:
            by_time[node["died"]].append((-1, tag))
    live = {tag for tag, node in nodes.items() if node["channel"] == "initial"}
    previous = both = 0.
    first_inactive = None if founders else 0.
    focal_loss = focal_extinct = establish = None
    restores = maximum = 0
    focal_time = nodes[focal]["born"] if focal is not None else None
    for time in sorted(by_time):
        current_species = {nodes[tag]["species"] for tag in live}
        both += (time - previous) * int(1 in current_species and 2 in current_species)
        previous = time
        # Exact equal-time multiple physical jumps have probability zero; preserve
        # their ledger order if synthetic fixtures intentionally contain them.
        for sign, tag in by_time[time]:
            old_r1 = bool(live & focal_descendants & by_species[2])
            if sign > 0:
                live.add(tag)
            else:
                live.remove(tag)
            focal_live = live & focal_descendants
            current_r1 = bool(focal_live & by_species[2])
            if focal is not None and time >= focal_time:
                if not current_r1 and focal_loss is None:
                    focal_loss = time
                if current_r1 and not old_r1 and time > focal_time:
                    restores += 1
                maximum = max(maximum, len(focal_live))
                if len(focal_live) >= threshold and establish is None:
                    establish = time
                if not focal_live and focal_extinct is None:
                    focal_extinct = time
            if first_inactive is None and not any(nodes[t]["species"] in (1, 2) for t in live):
                first_inactive = time
    final_species = {nodes[tag]["species"] for tag in live}
    both += (horizon - previous) * int(1 in final_species and 2 in final_species)
    assert live == final_tags
    cleared = all(nodes[tag]["removed"] for tag in founders)
    clearance = max((nodes[tag]["died"] for tag in founders), default=0.) if cleared else None
    born_r1 = [tag for tag, node in nodes.items() if node["species"] == 2]
    first_r1 = min(born_r1, key=lambda tag: (nodes[tag]["born"], tag)) if born_r1 else None
    births = Counter()
    for tag, node in nodes.items():
        channel, species, helper = node["channel"], node["species"], node["helper"]
        if channel == "copy_R":
            births["paid_R"] += 1
            births[f"paid_R{helper-1}_to_R{species-1}"] += 1
            if tag in focal_descendants:
                births["focal_R0_to_R1" if helper == 1 and species == 2 else "focal_other_copy"] += 1
            elif helper == 1 and species == 2:
                births["outside_focal_R0_to_R1"] += 1
        elif channel == "produce_L":
            births["paid_L"] += 1
        elif channel == "passive":
            births["thermal_R" if species in (1, 2) else "thermal_L"] += 1
    work = math.fsum(e["external_work"] for e in run["ledger"])
    heat = math.fsum(e["heat_to_bath"] for e in run["ledger"])
    ref = make_reference(run)
    first_law = abs(ref.energy(run["final_state"]) - ref.energy(run["initial_state"]) - work + heat)
    return dict(integral_R=total(by_species[1] | by_species[2]), integral_R0=total(by_species[1]),
        integral_R1=total(by_species[2]), integral_L=total(by_species[3]), both_types_present_time=both,
        integral_founder_descendant_R1=total(founder_descendants & by_species[2]),
        integral_focal_descendants=total(focal_descendants), integral_focal_R1=total(focal_descendants & by_species[2]),
        focal_paid_mutant_appeared_by_deadline=focal is not None,
        first_R1_appearance_time=nodes[first_r1]["born"] if first_r1 is not None else None,
        first_R1_appearance_origin=nodes[first_r1]["channel"] if first_r1 is not None else None,
        focal_tag=focal, focal_appearance_time=focal_time,
        focal_parent_has_initial_founder_ancestry=(nodes[focal]["parent"] in founder_descendants) if focal is not None else None,
        focal_established=establish is not None, establishment_population="copying descendants over both current alleles",
        focal_establishment_time=establish, focal_extinction_time=focal_extinct,
        focal_survives_to_horizon=bool(final_tags & focal_descendants), maximum_focal_descendants=maximum,
        focal_R1_first_loss_time=focal_loss, focal_R1_uninterrupted_to_horizon=focal is not None and focal_loss is None,
        focal_R1_restorations_after_loss=restores, final_focal_descendants=len(final_tags & focal_descendants),
        final_focal_R1=len(final_tags & focal_descendants & by_species[2]),
        founder_instance_clearance_time=clearance,
        paid_R_births_after_founder_instance_clearance=sum(node["channel"] == "copy_R" and clearance is not None and node["born"] >= clearance for node in nodes.values()),
        final_R=final_counts[1] + final_counts[2], final_R1=final_counts[2], final_L=final_counts[3],
        first_R_inactivity=first_inactive, first_episode_duration=horizon if first_inactive is None else first_inactive,
        chemical_work=work, heat_to_bath=heat, first_law_error=first_law, gate_events=event_counts["gate_flip"],
        gate_time=dict(gate_times), birth_counts=dict(births), event_counts=dict(event_counts))


def estimate(values):
    values = [float(x) for x in values]
    n, mean = len(values), statistics.mean(values)
    se = statistics.stdev(values) / math.sqrt(n) if n > 1 else 0.
    half = float(t.ppf(.975, n - 1)) * se if n > 1 else 0.
    return dict(n=n, mean=mean, monte_carlo_standard_error=se,
                degrees_freedom=n-1, t_95_interval=[mean-half, mean+half])


def vectors(rows, field):
    paired = defaultdict(dict)
    for row in rows:
        key = (row["response_mode"], row["produce_L"])
        assert key not in paired[row["origin_seed"]]
        paired[row["origin_seed"]][key] = row[field]
    one, two = [], []
    for _, arms in sorted(paired.items()):
        assert set(arms) == {(m, p) for m in ("opposed", "aligned") for p in (False, True)}
        response = arms[("opposed", True)] - arms[("aligned", True)]
        one.append(response)
        two.append(response - arms[("opposed", False)] + arms[("aligned", False)])
    return one, two


def check_statistics(summary, rows):
    contrasts, groups = [], []
    for epsilon, kappa in sorted({(r["asymmetry"], r["flip_rate"]) for r in rows}):
        selected = [r for r in rows if (r["asymmetry"], r["flip_rate"]) == (epsilon, kappa)]
        for mode in ("opposed", "aligned"):
            for produce in (True, False):
                arm = [r for r in selected if (r["response_mode"], r["produce_L"]) == (mode, produce)]
                groups.append(dict(asymmetry=epsilon, flip_rate=kappa, response_mode=mode, produce_L=produce,
                    statistics={field: estimate([r[field] for r in arm]) for field in FIELDS},
                    focal_appearances=sum(r["focal_paid_mutant_appeared_by_deadline"] for r in arm),
                    focal_establishments_over_both_alleles=sum(r["focal_established"] for r in arm),
                    focal_survivors=sum(r["focal_survives_to_horizon"] for r in arm),
                    final_zero_R=sum(r["final_R"] == 0 for r in arm)))
        for field in FIELDS:
            for name, values in zip(PRIMARY, vectors(selected, field)):
                contrasts.append(dict(asymmetry=epsilon, flip_rate=kappa, contrast=name, field=field, **estimate(values)))
    compare_values(groups, summary["groups"])
    compare_values(contrasts, summary["paired_contrasts"])
    primary = [r for r in contrasts if r["field"] == "integral_R"]
    compare_values(primary, summary["primary_results"])
    return primary


def check_selection(pilot, protocol):
    scores = []
    primary = pilot["primary_results"]
    for epsilon, kappa in sorted({(r["asymmetry"], r["flip_rate"]) for r in primary}):
        cells = [r for r in primary if (r["asymmetry"], r["flip_rate"]) == (epsilon, kappa)]
        assert len(cells) == 2
        z = [r["mean"] / r["monte_carlo_standard_error"] if r["monte_carlo_standard_error"] else
             0. if r["mean"] == 0 else math.copysign(1e12, r["mean"]) for r in cells]
        scores.append(dict(asymmetry=epsilon, flip_rate=kappa, score=min(z), component_scores=z))
    compare_values(scores, protocol["selection_scores"])
    winner = min(scores, key=lambda r: (-r["score"], r["asymmetry"], r["flip_rate"]))
    assert protocol["settings"] == [dict(asymmetry=winner["asymmetry"], flip_rate=winner["flip_rate"], initial_gate=None)]
    return winner


def verify_summary(path):
    summary = load(path)
    protocol_path = ROOT / summary["protocol_file"]
    protocol = load(protocol_path)
    assert digest(protocol_path) == summary["protocol_sha256"]
    assert summary["source_hashes"] == protocol["source_hashes"]
    for name, expected in protocol["source_hashes"].items():
        assert digest(ROOT / name) == expected, name
    assert protocol["params"]["mutation"] == [[.95, .05], [.05, .95]]
    assert protocol["params"]["initial_R1_fraction"] == 0.
    assert protocol["horizon"] == 100. and protocol["sample_dt"] == 2.
    assert protocol["focal_appearance_deadline"] == 20. and protocol["focal_establishment_size"] == 4
    assert protocol["primary_contrasts"] == list(PRIMARY)
    assert protocol["response_modes"] == ["opposed", "aligned"] and protocol["production_flags"] == [True, False]
    expected_cells = {(s["asymmetry"], s["flip_rate"], seed, mode, production)
                      for s in protocol["settings"] for seed in protocol["seeds"]
                      for mode in protocol["response_modes"] for production in protocol["production_flags"]}
    seen, files, matched = set(), set(), {}
    audits, independent_rows = [], []
    for number, row in enumerate(summary["runs"]):
        key = (row["asymmetry"], row["flip_rate"], row["origin_seed"], row["response_mode"], row["produce_L"])
        assert key in expected_cells and key not in seen
        seen.add(key)
        trajectory_path = ROOT / row["trajectory_file"]
        assert trajectory_path not in files
        files.add(trajectory_path)
        assert digest(trajectory_path) == row["trajectory_sha256"]
        run = load(trajectory_path)
        compare_values(protocol["params"], run["params"])
        compare_values(protocol["drive"], run["drive"])
        assert run["seed"] == row["origin_seed"] and run["horizon"] == protocol["horizon"]
        assert run["response"] == dict(asymmetry=row["asymmetry"], flip_rate=row["flip_rate"], initial_gate=None)
        assert run["response_mode"] == row["response_mode"] and run["transport_mode"] == "local"
        assert run["flags"] == dict(produce_L=row["produce_L"], protect_R=True, copy_R=True)
        assert run["production_site_mobility"] == [1.] * protocol["params"]["side"] ** 2
        initial = run["initial_state"]
        assert initial["base_state"]["time"] == 0. and initial["base_state"]["event_count"] == 0
        assert initial["base_state"]["cumulative_heat"] == initial["base_state"]["external_work"] == 0.
        assert initial["gate_event_count"] == 0
        assert set(initial["base_state"]["material"]) <= {0, 1}
        assert initial["base_state"]["material"].count(1) == protocol["params"]["initial_templates"]
        assert sum(initial["base_state"]["fuel"]) == protocol["params"]["initial_fuel"]
        gate_path = [(e["time"], e["before"], e["after"]) for e in run["ledger"] if e["kind"] == "gate_flip"]
        seeded_gate = int(np.random.default_rng(np.random.SeedSequence([run["seed"], 62411])).choice((-1, 1)))
        assert initial["gate"] == seeded_gate
        clock_rng = np.random.default_rng(np.random.SeedSequence([run["seed"], 62412]))
        expected_gate_path = []
        next_time = float(clock_rng.exponential(1. / row["flip_rate"]))
        while next_time <= run["horizon"]:
            expected_gate_path.append((next_time, seeded_gate, -seeded_gate))
            seeded_gate = -seeded_gate
            next_time += float(clock_rng.exponential(1. / row["flip_rate"]))
        assert gate_path == expected_gate_path
        pair_key = key[:3]
        if pair_key in matched:
            assert matched[pair_key] == (initial, gate_path)
        else:
            matched[pair_key] = (initial, gate_path)
        audit = verify_trajectory(run)
        measured = graph_readouts(run, protocol["focal_appearance_deadline"], protocol["focal_establishment_size"])
        compare_values(measured, row)
        independent_rows.append({**row, **measured})
        audits.append(dict(trajectory_file=row["trajectory_file"], trajectory_sha256=row["trajectory_sha256"], **audit))
        if (number + 1) % 32 == 0:
            print(f"{path.name}: independently replayed {number+1}/{len(summary['runs'])} trajectories", flush=True)
    assert seen == expected_cells
    assert summary["all_preparations_included"] is True
    primary = check_statistics(summary, independent_rows)
    selected = None
    if summary["stage"] == "confirmation":
        pilot = load(ROOT / "pilot-results.json")
        original = load(ROOT / "pilot-protocol.json")
        assert digest(ROOT / "pilot-results.json") == protocol["pilot_results_sha256"]
        assert digest(ROOT / "pilot-protocol.json") == protocol["pilot_protocol_sha256"]
        assert protocol["seeds"] == original["confirmation_seeds"] == list(range(1060001, 1060129))
        assert not set(protocol["seeds"]) & set(original["seeds"])
        selected = check_selection(pilot, protocol)
        spec = protocol["primary_bootstrap_sensitivity"]
        assert spec["replicates"] == 10000 and spec["seed"] == 1080001
        data = np.column_stack(vectors(independent_rows, "integral_R"))
        sampled = np.random.default_rng(spec["seed"]).integers(0, len(data), size=(spec["replicates"], len(data)))
        boot = np.percentile(np.mean(data[sampled], axis=1), (2.5, 97.5), axis=0)
        bootstrap = [dict(contrast=name, replicates=spec["replicates"], seed=spec["seed"],
                          percentile_95_interval=boot[:, i].tolist()) for i, name in enumerate(PRIMARY)]
        compare_values(bootstrap, summary["primary_bootstrap_sensitivity"])
        assert summary["declared_joint_claim_supported"] == all(r["t_95_interval"][0] > 0 for r in primary)
    else:
        assert protocol["seeds"] == list(range(1050001, 1050009))
        assert protocol["settings"] == [dict(asymmetry=e, flip_rate=k, initial_gate=None) for e in (.5, .8) for k in (.025, .1)]
        assert summary["declared_joint_claim_supported"] is None
    return dict(stage=summary["stage"], summary_file=path.name, summary_sha256=digest(path),
                protocol_file=protocol_path.name, protocol_sha256=digest(protocol_path),
                physical_records=len(audits), physical_events=sum(r["events_replayed"] for r in audits),
                samples_checked=sum(r["samples_checked"] for r in audits),
                complete_unfiltered_cells=True, matched_actual_gate_paths=True, matched_initial_states=True,
                matched_seed_setting_pairs=len(matched), fresh_confirmation=summary["stage"] == "confirmation",
                primary_results=primary, selected_pilot_setting=selected,
                declared_joint_claim_supported=summary["declared_joint_claim_supported"],
                final_zero_R=sum(r["final_R"] == 0 for r in independent_rows),
                no_focal_mutant_by_deadline=sum(not r["focal_paid_mutant_appeared_by_deadline"] for r in independent_rows),
                focal_establishments_over_both_alleles=sum(r["focal_established"] for r in independent_rows),
                focal_survivors=sum(r["focal_survives_to_horizon"] for r in independent_rows),
                focal_parents_without_initial_founder_ancestry=sum(r["focal_parent_has_initial_founder_ancestry"] is False for r in independent_rows),
                independent_lifetime_graph_readouts_match=True, trajectories=audits)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("summaries", nargs="+")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    results = [verify_summary(ROOT / name) for name in args.summaries]
    assert len({r["trajectory_file"] for stage in results for r in stage["trajectories"]}) == sum(stage["physical_records"] for stage in results)
    report = dict(checks_passed=True, verifier_sha256=digest(__file__),
        physics_verifier_sha256=digest(ROOT / "verify_differentiated.py"),
        records=sum(stage["physical_records"] for stage in results),
        events=sum(stage["physical_events"] for stage in results),
        samples=sum(stage["samples_checked"] for stage in results), stages=results,
        scope="Independent scalar full-trajectory and lifetime-graph audit. Primary records only; excludes synthetic cases and tiny preflight smokes. This verifies the finite supplied-chemistry experiment, not selection, efficacy of a particular mutant, collective individuality, global biodiversity or reduced work.")
    Path(args.output).write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({key: report[key] for key in ("checks_passed", "records", "events", "samples")}, indent=2))


if __name__ == "__main__":
    main()
