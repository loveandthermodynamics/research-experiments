"""Bounded pilot and fresh confirmation of physical protective-material renewal.

Imports frozen sibling sources, writes only this new study directory, retains
all preparations, and never uses a readout to alter a trajectory's rates.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parent
PRIOR = ROOT.parent / "interface-maintenance"
sys.path.insert(0, str(PRIOR))
from interface_model import Flags, Model, Params, simulate
from interface_readouts import estimate
from run_interface import measurements, save_run

ARMS = tuple(Flags(produce_L=p, protect_R=b) for p, b in
             ((True, True), (True, False), (False, True), (False, False)))
CONTRASTS = {
    "production_effect_when_protected": {(True, True): 1, (False, True): -1},
    "production_effect_without_protection": {(True, False): 1, (False, False): -1},
    "protection_effect_with_production": {(True, True): 1, (True, False): -1},
    "protection_effect_after_production_stopped": {(False, True): 1, (False, False): -1},
    "production_protection_interaction": {(True, True): 1, (True, False): -1,
                                          (False, True): -1, (False, False): 1}}
SOURCES = {"run_turnover.py": ROOT / "run_turnover.py",
           "../interface-maintenance/interface_model.py": PRIOR / "interface_model.py",
           "../interface-maintenance/interface_readouts.py": PRIOR / "interface_readouts.py",
           "../interface-maintenance/run_interface.py": PRIOR / "run_interface.py"}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def source_hashes() -> dict[str, str]:
    return {name: digest(path) for name, path in SOURCES.items()}


def write_once(path: Path, value: dict) -> None:
    encoded = json.dumps(value, indent=2) + "\n"
    if path.exists():
        assert path.read_text() == encoded, f"refusing to replace frozen {path.name}"
    else:
        path.write_text(encoded)


def freeze_pilot() -> dict:
    prior_protocol = PRIOR / "ensemble-protocol.json"
    params = json.loads(prior_protocol.read_text())["params"]
    protocol = {
        "stage": "pilot", "source_hashes": source_hashes(),
        "base_parameter_protocol": "../interface-maintenance/ensemble-protocol.json",
        "base_parameter_protocol_sha256": digest(prior_protocol), "base_params": params,
        "passive_L_levels": [0.08, 0.16, 0.32],
        "turnover_intervention": "Only passive_L changes across levels. It is the symmetric D<->L mobility during BOTH preparation and continuation and is equal across all factorial arms at a level.",
        "seeds": list(range(920001, 920009)), "replicates_per_arm_per_level": 8,
        "preparation_horizon": 20.0, "continuation_horizon": 80.0,
        "post_seed_offset": 3_000_000, "sample_dt": 2.0, "arms": [asdict(f) for f in ARMS],
        "selection_rule": "Choose exactly one level maximizing the minimum of paired R-time mean/SE for production_effect_when_protected and production_protection_interaction; ties choose lower passive_L. Selection occurs regardless of signs. No additional tuning.",
        "confirmation_seeds": list(range(930001, 930065)),
        "joint_support_rule": "Fresh confirmation supports the joint claim only if BOTH primary paired R-time contrasts have approximate 95 percent normal intervals strictly above zero.",
        "primary_contrasts": ["production_effect_when_protected", "production_protection_interaction"],
        "primary_field": "integral_R",
        "secondary_fields": ["integral_L", "paid_L_births", "passive_L_births", "passive_L_deaths",
                             "paid_R_births_after_original_founder_clearance",
                             "paid_R_births_after_preparation_R_clearance", "net_fuel_consumed", "net_heat_to_bath"],
        "preparation_rule": "Every seed starts dispersed with no L and runs full P1B1 for fixed time20; all resulting states, including R=0 or L=0, enter all four branches.",
        "branch_rule": "Identical material/fuel/instance_tags/next_tag across four arms. Reset only time, cumulative_heat, external_work and event_count. Copying stays enabled.",
        "scope": "Cross-level comparisons combine formation and continuation turnover changes. Primary causal contrasts are paired within a level. Physical turnover is not aging; instance ancestry is not functional information. No reward, replenishment, enclosing membrane or new individuality is supplied or inferred."
    }
    write_once(ROOT / "pilot-protocol.json", protocol)
    return protocol


def measured(run: dict) -> dict:
    """Frozen exact integrals; fix only initial-absence ancillary episode readout."""
    row = measurements(run)
    if not any(a in (1, 2) for a in run["initial_state"]["material"]):
        row["frozen_helper_first_episode_duration"] = row["first_episode_duration"]
        row["first_R_extinction"] = 0.0
        row["first_episode_duration"] = 0.0
        row["initial_absence_episode_readout_corrected"] = True
    else:
        row["initial_absence_episode_readout_corrected"] = False
    return row


def turnover_readouts(run: dict, original_founders: set[int]) -> dict:
    initial = run["initial_state"]
    prepared_r = {int(t) for a, t in zip(initial["material"], initial["instance_tags"]) if a in (1, 2)}
    prepared_l = {int(t) for a, t in zip(initial["material"], initial["instance_tags"]) if a == 3}
    groups = {"original_founder": set(original_founders), "preparation_R": prepared_r,
              "preparation_L": prepared_l}
    initial_tags = {int(t) for t in initial["instance_tags"] if t}
    remaining = {name: group & initial_tags for name, group in groups.items()}
    clearance = {name: 0.0 if not live else None for name, live in remaining.items()}
    copies_after = {name: 0 for name in groups}
    for event in run["ledger"]:
        removed = int(event["instance_removed"])
        for name in groups:
            remaining[name].discard(removed)
            if clearance[name] is None and not remaining[name]:
                clearance[name] = event["time"]
            if event["kind"] == "copy_R" and event["after"] in (1, 2) and clearance[name] is not None:
                copies_after[name] += 1
    return {"original_founder_clearance_time_post_branch": clearance["original_founder"],
            "preparation_R_clearance_time_post_branch": clearance["preparation_R"],
            "preparation_L_clearance_time_post_branch": clearance["preparation_L"],
            "paid_R_births_after_original_founder_clearance": copies_after["original_founder"],
            "paid_R_births_after_preparation_R_clearance": copies_after["preparation_R"],
            "paid_R_births_after_preparation_L_clearance": copies_after["preparation_L"]}


def summarize(rows: list[dict], fields: list[str]) -> tuple[list[dict], list[dict]]:
    groups, contrasts = [], []
    for level in sorted({row["passive_L"] for row in rows}):
        selected = [row for row in rows if row["passive_L"] == level]
        for flags in ARMS:
            arm = [row for row in selected if row["flags"] == asdict(flags)]
            groups.append({"passive_L": level, "flags": asdict(flags),
                           "statistics": {field: estimate([row[field] for row in arm]) for field in fields}})
        by_seed: dict[int, dict[tuple[bool, bool], dict]] = {}
        for row in selected:
            by_seed.setdefault(row["origin_seed"], {})[(row["flags"]["produce_L"], row["flags"]["protect_R"])] = row
        for name, weights in CONTRASTS.items():
            for field in fields:
                differences = [sum(weight * arms[arm][field] for arm, weight in weights.items())
                               for _, arms in sorted(by_seed.items())]
                contrasts.append({"passive_L": level, "contrast": name, "field": field,
                                  **estimate(differences),
                                  "definition": [{"produce_L": arm[0], "protect_R": arm[1], "weight": weight}
                                                 for arm, weight in weights.items()]})
    return groups, contrasts


def execute(stage: str) -> dict:
    protocol_path = ROOT / f"{stage}-protocol.json"
    assert protocol_path.exists(), "freeze protocol before simulation"
    protocol_hash = digest(protocol_path)
    protocol = json.loads(protocol_path.read_text())
    before = source_hashes()
    assert before == protocol["source_hashes"]
    assert digest(ROOT / protocol["base_parameter_protocol"]) == protocol["base_parameter_protocol_sha256"]
    if stage == "confirmation":
        assert digest(ROOT / protocol["pilot_results_file"]) == protocol["pilot_results_sha256"]
        assert digest(ROOT / "pilot-protocol.json") == protocol["pilot_protocol_sha256"]
    output = ROOT / stage
    output.mkdir(exist_ok=True)
    rows, preparations = [], []
    started = time.monotonic()
    for level in protocol["passive_L_levels"]:
        params = Params(**{**protocol["base_params"], "passive_L": level})
        label = f"L{int(round(level * 1000)):03d}"
        for seed in protocol["seeds"]:
            prep_path = output / f"{label}-s{seed}-preparation.json.gz"
            assert not prep_path.exists()
            prep = simulate(params, seed, protocol["preparation_horizon"], ARMS[0], sample_dt=protocol["sample_dt"])
            prep.update(stage=stage, phase="preparation", passive_L=level, source_hashes=before,
                        protocol_sha256=protocol_hash)
            save_run(prep_path, prep)
            prep_row = measured(prep)
            prep_row.update(passive_L=level, origin_seed=seed, phase="preparation",
                            trajectory_file=str(prep_path.relative_to(ROOT)), trajectory_sha256=digest(prep_path))
            preparations.append(prep_row)
            original_founders = {int(t) for a, t in zip(prep["initial_state"]["material"], prep["initial_state"]["instance_tags"]) if a in (1, 2)}
            model = Model(params)
            state = model.state_from_dict(prep["final_state"])
            post = model.copy_state(state)
            post.time = post.cumulative_heat = post.external_work = 0.0
            post.event_count = 0
            post_seed = seed + protocol["post_seed_offset"]
            expected_initial = model.state_to_dict(post)
            for flags in ARMS:
                path = output / f"{label}-s{seed}-p{int(flags.produce_L)}b{int(flags.protect_R)}.json.gz"
                assert not path.exists()
                run = simulate(params, post_seed, protocol["continuation_horizon"], flags,
                               sample_dt=protocol["sample_dt"], initial_state=post)
                assert run["initial_state"] == expected_initial
                run.update(stage=stage, phase="continuation", origin_seed=seed, passive_L=level,
                           preparation_file=str(prep_path.relative_to(ROOT)), preparation_sha256=digest(prep_path),
                           preparation_time=protocol["preparation_horizon"], source_hashes=before,
                           protocol_sha256=protocol_hash)
                save_run(path, run)
                row = measured(run)
                row.update(turnover_readouts(run, original_founders))
                row.update(passive_L=level, origin_seed=seed, phase="continuation", post_seed=post_seed,
                           preparation_file=str(prep_path.relative_to(ROOT)), preparation_sha256=digest(prep_path),
                           prepared_R=prep_row["final_R"], prepared_L=prep_row["final_L"],
                           prepared_fuel=sum(prep["final_state"]["fuel"]),
                           preparation_net_fuel_consumed=prep_row["net_fuel_consumed"],
                           preparation_heat_to_bath=prep_row["net_heat_to_bath"],
                           total_net_fuel_consumed=prep_row["net_fuel_consumed"] + row["net_fuel_consumed"],
                           total_heat_to_bath=prep_row["net_heat_to_bath"] + row["net_heat_to_bath"],
                           total_integral_R=prep_row["integral_R"] + row["integral_R"],
                           trajectory_file=str(path.relative_to(ROOT)), trajectory_sha256=digest(path))
                rows.append(row)
            print(json.dumps({"stage": stage, "passive_L": level, "origin_seed": seed,
                              "prepared_R": prep_row["final_R"], "prepared_L": prep_row["final_L"],
                              "prepared_fuel": sum(prep["final_state"]["fuel"]),
                              "elapsed_seconds": round(time.monotonic() - started, 2)}), flush=True)
    fields = ["integral_R", "integral_L", "final_R", "final_L", "paid_R_births", "paid_L_births",
              "passive_L_births", "passive_L_deaths", "net_fuel_consumed", "net_heat_to_bath",
              "first_episode_duration", "paid_R_births_after_original_founder_clearance",
              "paid_R_births_after_preparation_R_clearance", "paid_R_births_after_preparation_L_clearance",
              "total_net_fuel_consumed", "total_heat_to_bath", "total_integral_R"]
    groups, contrasts = summarize(rows, fields)
    assert len(preparations) == len(protocol["passive_L_levels"]) * len(protocol["seeds"])
    assert len(rows) == 4 * len(preparations)
    assert source_hashes() == before and digest(protocol_path) == protocol_hash
    result = {"stage": stage, "source_hashes": before, "protocol_file": protocol_path.name,
              "protocol_sha256": protocol_hash, "passive_L_levels": protocol["passive_L_levels"],
              "seeds": protocol["seeds"], "preparation_horizon": protocol["preparation_horizon"],
              "horizon": protocol["continuation_horizon"], "all_preparations_included": True,
              "replicates_per_arm_per_level": len(protocol["seeds"]), "preparations": preparations,
              "groups": groups, "paired_contrasts": contrasts, "runs": rows,
              "elapsed_seconds": time.monotonic() - started, "scope": protocol["scope"]}
    if stage == "confirmation":
        primary = [c for c in contrasts if c["field"] == "integral_R" and c["contrast"] in protocol["primary_contrasts"]]
        result["declared_joint_claim_supported"] = len(primary) == 2 and all(c["approximate_95_percent_normal_interval"][0] > 0 for c in primary)
        result["primary_confirmation_results"] = primary
    write_once(ROOT / f"{stage}-results.json", result)
    return result


def freeze_confirmation() -> dict:
    pilot_path = ROOT / "pilot-results.json"
    pilot_protocol = json.loads((ROOT / "pilot-protocol.json").read_text())
    pilot = json.loads(pilot_path.read_text())
    assert source_hashes() == pilot["source_hashes"] == pilot_protocol["source_hashes"]
    assert pilot["protocol_sha256"] == digest(ROOT / "pilot-protocol.json")
    scores = []
    for level in pilot_protocol["passive_L_levels"]:
        effects = {c["contrast"]: c for c in pilot["paired_contrasts"]
                   if c["passive_L"] == level and c["field"] == "integral_R"}
        zs = []
        for name in pilot_protocol["primary_contrasts"]:
            mean, se = effects[name]["mean"], effects[name]["monte_carlo_standard_error"]
            zs.append(mean / se if se else (1e300 if mean > 0 else -1e300 if mean < 0 else 0.0))
        scores.append({"passive_L": level, "primary_z_scores": zs, "selection_score": min(zs),
                       "primary_pilot_effects": [effects[name] for name in pilot_protocol["primary_contrasts"]]})
    selected = max(scores, key=lambda row: (row["selection_score"], -row["passive_L"]))["passive_L"]
    protocol = {**pilot_protocol, "stage": "confirmation", "passive_L_levels": [selected],
                "seeds": pilot_protocol["confirmation_seeds"], "replicates_per_arm_per_level": 64,
                "pilot_results_file": pilot_path.name, "pilot_results_sha256": digest(pilot_path),
                "pilot_protocol_sha256": digest(ROOT / "pilot-protocol.json"), "selection_scores": scores,
                "selected_passive_L": selected}
    write_once(ROOT / "confirmation-protocol.json", protocol)
    return protocol


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("freeze-pilot", "pilot", "freeze-confirmation", "confirmation"))
    args = parser.parse_args()
    if args.action == "freeze-pilot":
        protocol = freeze_pilot()
        print(json.dumps({"protocol": "pilot-protocol.json", "sha256": digest(ROOT / "pilot-protocol.json"),
                          "runner_sha256": digest(ROOT / "run_turnover.py")}))
    elif args.action == "freeze-confirmation":
        protocol = freeze_confirmation()
        print(json.dumps({"protocol": "confirmation-protocol.json", "sha256": digest(ROOT / "confirmation-protocol.json"),
                          "selected_passive_L": protocol["selected_passive_L"], "selection_scores": protocol["selection_scores"]}, indent=2))
    else:
        execute(args.action)
