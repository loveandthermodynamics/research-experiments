"""Predeclared driven contact-maintenance pilot and fresh confirmation."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from hashlib import sha256
import gzip
import json
from pathlib import Path
import time

import numpy as np

from driven_model import DriveParams, DrivenModel, Flags, Params, PRIOR, simulate_driven
from interface_readouts import estimate, spatial_readouts

ROOT = Path(__file__).resolve().parent
ARMS = tuple(Flags(produce_L=p, protect_R=b) for p, b in
             ((True, True), (True, False), (False, True), (False, False)))
MODES = ("local", "mixed_L")
SOURCES = {"driven_model.py": ROOT / "driven_model.py", "run_driven.py": ROOT / "run_driven.py",
           "model-spec.md": ROOT / "model-spec.md",
           "../interface-maintenance/interface_model.py": PRIOR / "interface_model.py",
           "../interface-maintenance/interface_readouts.py": PRIOR / "interface_readouts.py"}
PRIMARY = ("production_effect_when_protected", "production_protection_interaction", "locality_effect_full")


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def hashes():
    return {name: digest(path) for name, path in SOURCES.items()}


def write_once(path, value):
    encoded = json.dumps(value, indent=2) + "\n"
    if path.exists():
        assert path.read_text() == encoded, f"refusing to replace {path.name}"
    else:
        path.write_text(encoded)


def save_run(path, run):
    assert not path.exists()
    with gzip.open(path, "wt", encoding="utf-8") as stream:
        json.dump(run, stream, separators=(",", ":"))


def freeze_pilot():
    params = asdict(Params(side=8, theta=0.7, epsilon_R=8, epsilon_L=8, E_F=12,
                          J_LL=1, J_RL=1, J_RR=0, copy_rate=0.025, production_rate=0.04,
                          passive_R=0.12, passive_L=0.08, material_mobility=0.1, fuel_mobility=0.2,
                          protection_strength=1, transport_protection_strength=0.25,
                          initial_templates=8, initial_fuel=58))
    protocol = {"stage": "pilot", "source_hashes": hashes(), "base_params": params,
                "copy_rate_levels": [0.020, 0.025, 0.028], "drive": asdict(DriveParams()),
                "seeds": list(range(1010001, 1010009)), "confirmation_seeds": list(range(1020001, 1020065)),
                "replicates_per_cell": 8, "preparation_horizon": 20.0, "continuation_horizon": 80.0,
                "post_seed_offset": 4_000_000, "sample_dt": 2.0,
                "arms": [asdict(a) for a in ARMS], "transport_modes": list(MODES),
                "preparation_rule": "Each seed starts eight dispersed R0, zero L and58F. Local full P1B1 evolves for20; every resulting state enters all8branches without selection or replenishment.",
                "branch_rule": "Identical material/fuel/tags/next_tag for every mode and flag. Reset only time, cumulative_heat, external_work,event_count. Transport graph changes only after local preparation.",
                "reservoir_rule": "Uniform reversible W<->F exchange at every site, up=.9,down=.1,mu=E_F+theta*ln9; count reservoir work and bath heat on each jump.",
                "mixing_control": "Mixed_L replaces L-containing NNswaps by all-pairs conservative exchanges with proposal mobility4M/(N-1), same H/Metropolis and endpoint-inclusive symmetric barrier. Matched proposal clocks, not accepted flux or R displacement. No mixed formation comparison.",
                "primary_field": "integral_R", "primary_contrasts": list(PRIMARY),
                "selection_rule": "Choose one copy rate maximizing the minimum of the three paired primary mean/SE scores; ties choose lower copy rate. Retain all levels and signs; no extra tuning.",
                "confirmation_rule": "64 fresh seeds at selected rate. Support joint claim only if ALL THREE paired primary approximate95percent normal intervals have lower endpoints strictly positive.",
                "secondary_readouts": ["integral_L", "paid_L_births", "passive_L_births", "sampled_integral_RL_contacts",
                                       "sampled_mean_RL_enrichment", "final_spatial", "paid_R_births_after_preparation_R_clearance",
                                       "paid_L_births_after_preparation_L_clearance", "chemical_work", "heat_to_bath"],
                "scope": "Driven finite-horizon contact maintenance. Initial material/fuel/H and reservoir/transport apparatus are supplied. Uniform drive is not a group reward; contact geometry and instance continuity do not establish enclosing boundaries, aging, functional memory, participant independence or individuality. Heat is not total entropy production; fuel/work ratios alone do not establish matched-function efficiency."}
    write_once(ROOT / "pilot-protocol.json", protocol)
    return protocol


def measurements(run, original_founders=None):
    initial, final = run["initial_state"], run["final_state"]
    m = initial["material"]
    counts = {"R": sum(a in (1, 2) for a in m), "L": m.count(3)}
    integral = {"R": 0.0, "L": 0.0}
    event_counts = {name: 0 for name in ("paid_R_births", "reverse_R_copies", "paid_L_births", "reverse_L_production",
                                        "passive_R_births", "passive_R_deaths", "passive_L_births", "passive_L_deaths",
                                        "reservoir_fuel_up", "reservoir_fuel_down", "mixed_L_swaps")}
    initial_r_tags = {t for a, t in zip(m, initial["instance_tags"]) if a in (1, 2)}
    initial_l_tags = {t for a, t in zip(m, initial["instance_tags"]) if a == 3}
    tracked = {"original_founder": set(original_founders or initial_r_tags), "preparation_R": initial_r_tags,
               "preparation_L": initial_l_tags}
    present = {t for t in initial["instance_tags"] if t}
    live = {key: tags & present for key, tags in tracked.items()}
    cleared = {key: 0.0 if not tags else None for key, tags in live.items()}
    post_r = {key: 0 for key in tracked}
    post_l = {key: 0 for key in tracked}
    paid_r, paid_l = set(), set()
    first_inactivity = 0.0 if counts["R"] == 0 else None
    previous = work = heat = 0.0
    for event in run["ledger"]:
        dt = event["time"] - previous
        assert dt >= 0
        for key in counts:
            integral[key] += counts[key] * dt
        previous = event["time"]
        work += event["external_work"]
        heat += event["heat_to_bath"]
        kind, before, after = event["kind"], event["before"], event["after"]
        if kind == "reservoir_exchange":
            event_counts["reservoir_fuel_up" if after == 1 else "reservoir_fuel_down"] += 1
        elif kind == "mixed_L_swap":
            event_counts["mixed_L_swaps"] += 1
        elif kind in ("passive", "copy_R", "produce_L"):
            counts["R"] += int(after in (1, 2)) - int(before in (1, 2))
            counts["L"] += int(after == 3) - int(before == 3)
            if kind == "copy_R":
                event_counts["paid_R_births" if after in (1, 2) else "reverse_R_copies"] += 1
                if after in (1, 2):
                    paid_r.add(event["instance_created"])
            elif kind == "produce_L":
                event_counts["paid_L_births" if after == 3 else "reverse_L_production"] += 1
                if after == 3:
                    paid_l.add(event["instance_created"])
            else:
                is_r = before in (1, 2) or after in (1, 2)
                name = ("passive_R_" if is_r else "passive_L_") + ("births" if before == 0 else "deaths")
                event_counts[name] += 1
        for key in tracked:
            live[key].discard(event["instance_removed"])
            if cleared[key] is None and not live[key]:
                cleared[key] = event["time"]
            if cleared[key] is not None and kind == "copy_R" and after in (1, 2):
                post_r[key] += 1
            if cleared[key] is not None and kind == "produce_L" and after == 3:
                post_l[key] += 1
        if counts["R"] == 0 and first_inactivity is None:
            first_inactivity = event["time"]
    for key in counts:
        integral[key] += counts[key] * (run["horizon"] - previous)
    model = DrivenModel(Params(**run["params"]), Flags(**run["flags"]), DriveParams(**run["drive"]), run["transport_mode"])
    h0, h1 = model.energy(model.state_from_dict(initial)), model.energy(model.state_from_dict(final))
    error = abs(h1 - h0 - work + heat)
    assert error < 2e-7
    assert abs(final["external_work"] - initial["external_work"] - work) < 2e-7
    assert abs(final["cumulative_heat"] - initial["cumulative_heat"] - heat) < 2e-7
    fuel_delta = sum(final["fuel"]) - sum(initial["fuel"])
    internal_consumption = event_counts["paid_R_births"] + event_counts["paid_L_births"] - event_counts["reverse_R_copies"] - event_counts["reverse_L_production"]
    assert fuel_delta == event_counts["reservoir_fuel_up"] - event_counts["reservoir_fuel_down"] - internal_consumption
    fm, tags = np.asarray(final["material"]), np.asarray(final["instance_tags"])
    spatial = [spatial_readouts(s["material"], run["params"]["side"]) for s in run["samples"]]
    times = np.array([s["time"] for s in run["samples"]])
    rl = np.array([s["RL_edges"] for s in spatial])
    expected_rl = np.array([4 * s["R"] * s["L"] / (len(m) - 1) for s in spatial])
    enrichment = [s["RL_enrichment_over_random_placement"] for s in spatial if s["RL_enrichment_over_random_placement"] is not None]
    low_power = []
    side = run["params"]["side"]
    for sample in run["samples"]:
        occupied = np.isin(sample["material"], (1, 2)).reshape(side, side)
        nr = int(occupied.sum())
        if 0 < nr < len(m):
            transform = np.fft.fft2(occupied)
            power = np.mean([abs(transform[1, 0]) ** 2, abs(transform[-1, 0]) ** 2,
                             abs(transform[0, 1]) ** 2, abs(transform[0, -1]) ** 2])
            random_fixed_count_power = nr * (len(m) - nr) / (len(m) - 1)
            low_power.append(float(power / random_fixed_count_power))
    return {"seed": run["seed"], "flags": run["flags"], "transport_mode": run["transport_mode"],
            **event_counts, "integral_R": integral["R"], "integral_L": integral["L"],
            "final_R": counts["R"], "final_L": counts["L"], "initial_R": sum(a in (1, 2) for a in m),
            "initial_L": m.count(3), "initial_fuel": sum(initial["fuel"]), "final_fuel": sum(final["fuel"]),
            "initial_H": h0, "final_H": h1, "chemical_work": work, "heat_to_bath": heat,
            "first_law_error": error, "events": len(run["ledger"]),
            "first_R_inactivity": first_inactivity, "first_episode_duration": run["horizon"] if first_inactivity is None else first_inactivity,
            "paid_R_births_after_original_founder_clearance": post_r["original_founder"],
            "paid_R_births_after_preparation_R_clearance": post_r["preparation_R"],
            "paid_L_births_after_preparation_L_clearance": post_l["preparation_L"],
            "final_paid_R_instances": sum(int(t) in paid_r for t in tags[np.isin(fm, (1, 2))]),
            "final_paid_L_instances": sum(int(t) in paid_l for t in tags[fm == 3]),
            "founder_R_instances_remaining": sum(int(t) in tracked["original_founder"] for t in tags[np.isin(fm, (1, 2))]),
            "sampled_integral_RL_contacts": float(np.trapezoid(rl, times)),
            "sampled_integral_RL_excess_over_fixed_counts": float(np.trapezoid(rl - expected_rl, times)),
            "sampled_mean_RL_enrichment": float(np.mean(enrichment)) if enrichment else None,
            "sampled_mean_low_wave_Rpower_normalized": float(np.mean(low_power)) if low_power else None,
            "final_spatial": spatial_readouts(final["material"], run["params"]["side"])}


FIELDS = ("integral_R", "integral_L", "paid_L_births", "passive_L_births", "passive_L_deaths",
          "paid_R_births", "paid_R_births_after_original_founder_clearance",
          "paid_R_births_after_preparation_R_clearance", "paid_L_births_after_preparation_L_clearance",
          "sampled_integral_RL_contacts", "sampled_integral_RL_excess_over_fixed_counts", "chemical_work", "heat_to_bath", "final_R", "final_L",
          "first_episode_duration", "total_chemical_work", "total_heat_to_bath")


def summarize(rows):
    groups, contrasts = [], []
    for rate in sorted({r["copy_rate"] for r in rows}):
        selected = [r for r in rows if r["copy_rate"] == rate]
        for mode in MODES:
            for flags in ARMS:
                arm = [r for r in selected if r["transport_mode"] == mode and r["flags"] == asdict(flags)]
                morphology = {}
                for field in ("sampled_mean_RL_enrichment", "sampled_mean_low_wave_Rpower_normalized"):
                    valid = [r[field] for r in arm if r[field] is not None]
                    morphology[field] = estimate(valid) if valid else {"n": 0, "mean": None}
                groups.append({"copy_rate": rate, "transport_mode": mode, "flags": asdict(flags),
                               "statistics": {field: estimate([r[field] for r in arm]) for field in FIELDS},
                               "morphology_statistics_conditioned_on_defined_samples": morphology})
        by_seed = {}
        for row in selected:
            by_seed.setdefault(row["origin_seed"], {})[(row["transport_mode"], row["flags"]["produce_L"], row["flags"]["protect_R"])] = row
        weights = {
            "production_effect_when_protected": {("local", True, True): 1, ("local", False, True): -1},
            "production_protection_interaction": {("local", True, True): 1, ("local", True, False): -1,
                                                  ("local", False, True): -1, ("local", False, False): 1},
            "locality_effect_full": {("local", True, True): 1, ("mixed_L", True, True): -1}}
        for name, ws in weights.items():
            for field in FIELDS:
                values = [sum(w * cells[cell][field] for cell, w in ws.items()) for _, cells in sorted(by_seed.items())]
                contrasts.append({"copy_rate": rate, "contrast": name, "field": field, **estimate(values)})
    return groups, contrasts


def execute(stage):
    pp = ROOT / f"{stage}-protocol.json"
    protocol = json.loads(pp.read_text()); ph = digest(pp); source = hashes()
    assert source == protocol["source_hashes"]
    if stage == "confirmation":
        assert digest(ROOT / "pilot-results.json") == protocol["pilot_results_sha256"]
        assert digest(ROOT / "pilot-protocol.json") == protocol["pilot_protocol_sha256"]
    folder = ROOT / stage; folder.mkdir(exist_ok=True)
    preparations, rows = [], []; started = time.monotonic()
    drive = DriveParams(**protocol["drive"])
    for rate in protocol["copy_rate_levels"]:
        params = Params(**{**protocol["base_params"], "copy_rate": rate})
        label = f"k{int(round(rate * 1000)):03d}"
        for seed in protocol["seeds"]:
            prep = simulate_driven(params, seed, protocol["preparation_horizon"], ARMS[0], drive, "local", sample_dt=protocol["sample_dt"])
            prep.update(stage=stage, phase="preparation", source_hashes=source, protocol_sha256=ph, copy_rate=rate)
            path = folder / f"{label}-s{seed}-preparation.json.gz"; save_run(path, prep)
            prow = measurements(prep); prow.update(origin_seed=seed, copy_rate=rate, trajectory_file=str(path.relative_to(ROOT)), trajectory_sha256=digest(path))
            preparations.append(prow)
            model = DrivenModel(params, drive=drive)
            post = model.state_from_dict(prep["final_state"])
            post.time = post.cumulative_heat = post.external_work = 0.0; post.event_count = 0
            expected = model.state_to_dict(post)
            founders = {t for a, t in zip(prep["initial_state"]["material"], prep["initial_state"]["instance_tags"]) if a in (1, 2)}
            for mode in MODES:
                for flags in ARMS:
                    run = simulate_driven(params, seed + protocol["post_seed_offset"], protocol["continuation_horizon"], flags, drive, mode,
                                          initial_state=post, sample_dt=protocol["sample_dt"])
                    assert run["initial_state"] == expected
                    run.update(stage=stage, phase="continuation", origin_seed=seed, copy_rate=rate, source_hashes=source,
                               protocol_sha256=ph, preparation_file=str(path.relative_to(ROOT)), preparation_sha256=digest(path))
                    dest = folder / f"{label}-s{seed}-{mode}-p{int(flags.produce_L)}b{int(flags.protect_R)}.json.gz"; save_run(dest, run)
                    row = measurements(run, founders)
                    row.update(origin_seed=seed, copy_rate=rate, preparation_file=str(path.relative_to(ROOT)), preparation_sha256=digest(path),
                               prepared_R=prow["final_R"], prepared_L=prow["final_L"], prepared_fuel=prow["final_fuel"],
                               supplied_initial_H=prow["initial_H"], preparation_chemical_work=prow["chemical_work"],
                               preparation_heat_to_bath=prow["heat_to_bath"], total_chemical_work=prow["chemical_work"] + row["chemical_work"],
                               total_heat_to_bath=prow["heat_to_bath"] + row["heat_to_bath"],
                               trajectory_file=str(dest.relative_to(ROOT)), trajectory_sha256=digest(dest))
                    rows.append(row)
            print(json.dumps({"stage": stage, "copy_rate": rate, "origin_seed": seed, "prepared_R": prow["final_R"],
                              "prepared_L": prow["final_L"], "prep_work": round(prow["chemical_work"], 3),
                              "elapsed_seconds": round(time.monotonic() - started, 2)}), flush=True)
    assert len(preparations) == len(protocol["copy_rate_levels"]) * len(protocol["seeds"])
    assert len(rows) == 8 * len(preparations)
    groups, contrasts = summarize(rows)
    assert source == hashes() and digest(pp) == ph
    result = {"stage": stage, "source_hashes": source, "protocol_sha256": ph, "protocol_file": pp.name,
              "copy_rate_levels": protocol["copy_rate_levels"], "seeds": protocol["seeds"],
              "preparation_horizon": protocol["preparation_horizon"], "horizon": protocol["continuation_horizon"],
              "all_preparations_included": True, "preparations": preparations, "runs": rows, "groups": groups,
              "paired_contrasts": contrasts, "elapsed_seconds": time.monotonic() - started, "scope": protocol["scope"]}
    if stage == "confirmation":
        primary = [c for c in contrasts if c["field"] == "integral_R" and c["contrast"] in PRIMARY]
        result["primary_confirmation_results"] = primary
        result["declared_joint_claim_supported"] = len(primary) == 3 and all(c["approximate_95_percent_normal_interval"][0] > 0 for c in primary)
    write_once(ROOT / f"{stage}-results.json", result)
    return result


def freeze_confirmation():
    ppath = ROOT / "pilot-protocol.json"; p = json.loads(ppath.read_text())
    rpath = ROOT / "pilot-results.json"; r = json.loads(rpath.read_text())
    assert hashes() == r["source_hashes"] == p["source_hashes"] and r["protocol_sha256"] == digest(ppath)
    scores = []
    for rate in p["copy_rate_levels"]:
        effects = {c["contrast"]: c for c in r["paired_contrasts"] if c["copy_rate"] == rate and c["field"] == "integral_R"}
        zs = [effects[n]["mean"] / effects[n]["monte_carlo_standard_error"] if effects[n]["monte_carlo_standard_error"] else 0.0 for n in PRIMARY]
        scores.append({"copy_rate": rate, "primary_z_scores": zs, "selection_score": min(zs), "primary_pilot_effects": [effects[n] for n in PRIMARY]})
    selected = max(scores, key=lambda x: (x["selection_score"], -x["copy_rate"]))["copy_rate"]
    protocol = {**p, "stage": "confirmation", "copy_rate_levels": [selected], "seeds": p["confirmation_seeds"],
                "replicates_per_cell": 64, "selection_scores": scores, "selected_copy_rate": selected,
                "pilot_results_sha256": digest(rpath), "pilot_protocol_sha256": digest(ppath)}
    write_once(ROOT / "confirmation-protocol.json", protocol)
    return protocol


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("action", choices=("freeze-pilot", "pilot", "freeze-confirmation", "confirmation"))
    args = parser.parse_args()
    if args.action == "freeze-pilot":
        freeze_pilot(); print(json.dumps({"protocol_sha256": digest(ROOT / "pilot-protocol.json"), "source_hashes": hashes()}))
    elif args.action == "freeze-confirmation":
        p = freeze_confirmation(); print(json.dumps({"protocol_sha256": digest(ROOT / "confirmation-protocol.json"), "selected_copy_rate": p["selected_copy_rate"], "selection_scores": p["selection_scores"]}, indent=2))
    else:
        execute(args.action)
