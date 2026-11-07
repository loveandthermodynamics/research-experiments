"""One fixed, independent locality replication; no pooling or parameter tuning."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
import time

import numpy as np
from scipy.stats import t

from driven_model import DriveParams, DrivenModel, Flags, Params, simulate_driven
from run_driven import measurements, save_run

ROOT = Path(__file__).resolve().parent
FULL = Flags(produce_L=True, protect_R=True, copy_R=True)
PROTOCOL = ROOT / "locality-replication-protocol.json"
SOURCES = ("driven_model.py", "run_driven.py", "model-spec.md", "replicate_locality.py")


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def hashes():
    return {name: digest(ROOT / name) for name in SOURCES}


def write_once(path, value):
    encoded = json.dumps(value, indent=2) + "\n"
    if path.exists():
        assert path.read_text() == encoded, f"refusing to replace {path.name}"
    else:
        path.write_text(encoded)


def freeze_protocol():
    old_protocol_path = ROOT / "confirmation-protocol.json"
    old_results_path = ROOT / "confirmation-results.json"
    old = json.loads(old_protocol_path.read_text())
    params = {**old["base_params"], "copy_rate": old["selected_copy_rate"]}
    seeds = list(range(1040001, 1040257))
    assert not set(seeds) & set(old["seeds"])
    assert not set(seeds) & set(json.loads((ROOT / "pilot-protocol.json").read_text())["seeds"])
    protocol = {
        "stage": "locality_replication", "source_hashes": hashes(), "params": params, "drive": old["drive"],
        "prior_confirmation_protocol": old_protocol_path.name,
        "prior_confirmation_protocol_sha256": digest(old_protocol_path),
        "prior_confirmation_results": old_results_path.name,
        "prior_confirmation_results_sha256": digest(old_results_path),
        "seeds": seeds, "replicates_per_arm": 256, "flags": asdict(FULL),
        "transport_modes": ["local", "mixed_L"], "preparation_mode": "local",
        "preparation_horizon": 20.0, "continuation_horizon": 80.0,
        "post_seed_offset": 4_000_000, "sample_dt": 2.0,
        "preparation_rule": "Every seed begins eight dispersed R0, zero L and58F, evolves local full P1B1 to fixed20; retain every state, including failures.",
        "branch_rule": "Both full branches share identical material/fuel/instance_tags/next_tag; reset only time, cumulative_heat, external_work,event_count. Samepostseed, no material replenishment beyond the declared reservoir.",
        "primary_field": "integral_R", "primary_contrast": "local minus mixed_L, paired by origin seed",
        "primary_interval": "Two-sided Student t255 95percent mean interval; support only if lower endpoint is positive.",
        "bootstrap_replicates": 10_000, "bootstrap_seed": 10_000,
        "bootstrap_interval": "Percentile2.5/97.5 of paired-difference resampled means, NumPy default_rng(seed10000); sensitivity only.",
        "prior_result_rule": "Report prior confirmation separately with original normal criterion and a t63 sensitivity. Never pool prior and replication data or alter the previous decision.",
        "stopping_rule": "Exactly256freshpreparations and512branches; no interim primary analysis, tuning, additional sample extension or pooled rescue.",
        "scope": "One comparative transport claim under the same supplied drive and local preparation. Mixed graph changes accepted flux and R displacement despite matched proposal clocks. Neither result proves individuality, lower energetic cost, an enclosing membrane, aging or functional memory. Full preparation work, branch work and initial resources are retained; bath heat is not total entropy production."
    }
    write_once(PROTOCOL, protocol)
    return protocol


def previous_result(protocol):
    old = json.loads((ROOT / protocol["prior_confirmation_results"]).read_text())
    selected = next(row for row in old["primary_confirmation_results"] if row["contrast"] == "locality_effect_full")
    n = selected["n"]
    critical = float(t.ppf(0.975, n - 1))
    half = critical * selected["monte_carlo_standard_error"]
    return {"source_file": protocol["prior_confirmation_results"],
            "source_sha256": protocol["prior_confirmation_results_sha256"], "normal_result_as_saved": selected,
            "original_joint_criterion_supported_as_saved": old["declared_joint_claim_supported"],
            "t_sensitivity_degrees_freedom": n - 1, "t_critical": critical,
            "t_95_interval_sensitivity": [selected["mean"] - half, selected["mean"] + half],
            "interpretation": "Historical result and original decision retained unchanged; not pooled with replication."}


def execute():
    protocol = json.loads(PROTOCOL.read_text()); ph = digest(PROTOCOL); source = hashes()
    assert source == protocol["source_hashes"]
    for field in ("prior_confirmation_protocol", "prior_confirmation_results"):
        assert digest(ROOT / protocol[field]) == protocol[field + "_sha256"]
    folder = ROOT / "locality-replication"; folder.mkdir(exist_ok=True)
    params, drive = Params(**protocol["params"]), DriveParams(**protocol["drive"])
    rows, preparations, differences = [], [], []
    started = time.monotonic()
    for seed in protocol["seeds"]:
        prep = simulate_driven(params, seed, protocol["preparation_horizon"], FULL, drive, "local", sample_dt=protocol["sample_dt"])
        prep.update(stage="locality_replication", phase="preparation", source_hashes=source, protocol_sha256=ph)
        path = folder / f"s{seed}-preparation.json.gz"; save_run(path, prep)
        prow = measurements(prep)
        prow.update(origin_seed=seed, trajectory_file=str(path.relative_to(ROOT)), trajectory_sha256=digest(path))
        preparations.append(prow)
        model = DrivenModel(params, FULL, drive)
        post = model.state_from_dict(prep["final_state"])
        post.time = post.cumulative_heat = post.external_work = 0.0; post.event_count = 0
        expected = model.state_to_dict(post)
        founders = {tag for species, tag in zip(prep["initial_state"]["material"], prep["initial_state"]["instance_tags"]) if species in (1, 2)}
        pair = {}
        for mode in protocol["transport_modes"]:
            run = simulate_driven(params, seed + protocol["post_seed_offset"], protocol["continuation_horizon"], FULL, drive, mode,
                                  sample_dt=protocol["sample_dt"], initial_state=post)
            assert run["initial_state"] == expected
            run.update(stage="locality_replication", phase="continuation", origin_seed=seed,
                       preparation_file=str(path.relative_to(ROOT)), preparation_sha256=digest(path),
                       source_hashes=source, protocol_sha256=ph)
            dest = folder / f"s{seed}-{mode}.json.gz"; save_run(dest, run)
            row = measurements(run, founders)
            row.update(origin_seed=seed, preparation_file=str(path.relative_to(ROOT)), preparation_sha256=digest(path),
                       prepared_R=prow["final_R"], prepared_L=prow["final_L"], prepared_fuel=prow["final_fuel"],
                       supplied_initial_H=prow["initial_H"], preparation_chemical_work=prow["chemical_work"],
                       preparation_heat_to_bath=prow["heat_to_bath"], total_chemical_work=prow["chemical_work"] + row["chemical_work"],
                       total_heat_to_bath=prow["heat_to_bath"] + row["heat_to_bath"],
                       trajectory_file=str(dest.relative_to(ROOT)), trajectory_sha256=digest(dest))
            rows.append(row); pair[mode] = row["integral_R"]
        differences.append({"origin_seed": seed, "local_integral_R": pair["local"], "mixed_L_integral_R": pair["mixed_L"],
                            "difference": pair["local"] - pair["mixed_L"]})
        print(json.dumps({"origin_seed": seed, "prepared_R": prow["final_R"], "prepared_L": prow["final_L"],
                          "elapsed_seconds": round(time.monotonic() - started, 2)}), flush=True)
    assert len(preparations) == 256 and len(rows) == 512 and len(differences) == 256
    delta = np.array([row["difference"] for row in differences])
    n, mean = len(delta), float(delta.mean())
    se = float(delta.std(ddof=1) / np.sqrt(n)); critical = float(t.ppf(0.975, n - 1))
    interval = [mean - critical * se, mean + critical * se]
    rng = np.random.default_rng(protocol["bootstrap_seed"])
    bootstrap_means = delta[rng.integers(0, n, size=(protocol["bootstrap_replicates"], n))].mean(axis=1)
    bootstrap = np.quantile(bootstrap_means, [0.025, 0.975]).tolist()
    fields = ("integral_R", "integral_L", "chemical_work", "total_chemical_work", "heat_to_bath", "total_heat_to_bath",
              "sampled_integral_RL_contacts", "sampled_integral_RL_excess_over_fixed_counts", "paid_L_births",
              "paid_R_births_after_preparation_R_clearance", "paid_L_births_after_preparation_L_clearance", "final_R")
    groups = [{"transport_mode": mode, "means": {field: float(np.mean([r[field] for r in rows if r["transport_mode"] == mode])) for field in fields}}
              for mode in protocol["transport_modes"]]
    result = {"stage": "locality_replication", "source_hashes": source, "protocol_sha256": ph,
              "protocol_file": PROTOCOL.name, "seeds": protocol["seeds"], "all_preparations_included": True,
              "preparation_horizon": protocol["preparation_horizon"], "horizon": protocol["continuation_horizon"],
              "preparations": preparations, "runs": rows, "paired_differences": differences, "groups": groups,
              "replication_estimate": {"n": n, "mean": mean, "monte_carlo_standard_error": se,
                                       "degrees_freedom": n - 1, "t_critical": critical, "t_95_interval": interval,
                                       "primary_support_criterion_met": interval[0] > 0,
                                       "bootstrap_replicates": protocol["bootstrap_replicates"],
                                       "bootstrap_seed": protocol["bootstrap_seed"], "percentile_bootstrap_95_interval": bootstrap},
              "previous_confirmation": previous_result(protocol), "pooled_estimate": "Not computed or authorized.",
              "elapsed_seconds": time.monotonic() - started, "scope": protocol["scope"]}
    assert hashes() == source and digest(PROTOCOL) == ph
    for field in ("prior_confirmation_protocol", "prior_confirmation_results"):
        assert digest(ROOT / protocol[field]) == protocol[field + "_sha256"]
    write_once(ROOT / "locality-replication-results.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--freeze-protocol", action="store_true")
    args = parser.parse_args()
    if args.freeze_protocol:
        freeze_protocol(); print(json.dumps({"protocol_sha256": digest(PROTOCOL), "source_hashes": hashes()}))
    else:
        execute()
