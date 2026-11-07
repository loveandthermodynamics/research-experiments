"""Fixed-time continuation assay: ongoing production versus legacy protection.

Owns maintenance-protocol.json, maintenance-results.json and maintenance/*.json.gz.
The frozen base model, original runner and original trajectories are read-only.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import gzip
from hashlib import sha256
import json
from pathlib import Path
import time

import numpy as np

from interface_model import Flags, Model, Params, simulate
from interface_readouts import estimate
from run_interface import measurements, save_run

ROOT = Path(__file__).resolve().parent
PROTOCOL = ROOT / "maintenance-protocol.json"
ARMS = tuple(Flags(produce_L=p, protect_R=b) for p, b in
             ((True, True), (True, False), (False, True), (False, False)))
SOURCE_NAMES = ("interface_model.py", "interface_readouts.py", "run_interface.py", "run_maintenance.py")


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def source_hashes() -> dict[str, str]:
    return {name: digest(ROOT / name) for name in SOURCE_NAMES}


def freeze_protocol() -> dict:
    source = ROOT / "ensemble-results.json"
    ensemble = json.loads(source.read_text())
    selected = [row for row in ensemble["runs"] if row["flags"] == asdict(ARMS[0])]
    selected.sort(key=lambda row: row["seed"])
    assert len(selected) == 64
    assert [row["seed"] for row in selected] == sorted(ensemble["seeds"])
    preparations = [{"origin_seed": row["seed"], "source_trajectory": row["trajectory_file"],
                     "source_trajectory_sha256": digest(ROOT / row["trajectory_file"])} for row in selected]
    params = ensemble["params_by_case"]["selected"]
    protocol = {
        "stage": "maintenance", "source_results": source.name, "source_results_sha256": digest(source),
        "source_hashes": source_hashes(), "preparations": preparations,
        "replicates_per_arm": 64, "number_arms": 4, "branch_time": 20.0,
        "continuation_horizon": 80.0, "post_seed_offset": 3_000_000, "sample_dt": 2.0,
        "params": params, "arms": [asdict(flags) for flags in ARMS],
        "selection": "Every original full-production/protection baseline seed, including failed or low-L preparations; no state-dependent selection.",
        "preparation": "Replay the original ledger from its original microscopic initial state through every event at time <= 20. No extrapolated or selected successful state.",
        "bookkeeping_reset": ["time", "cumulative_heat", "external_work", "event_count"],
        "physical_reset": "None. Material, fuel, instance tags and next_tag remain identical across the four continuations.",
        "primary_fields": ["integral_R", "integral_L"],
        "primary_paired_contrasts": ["production_effect_when_protected", "production_effect_without_protection",
                                     "protection_effect_with_production", "protection_effect_after_production_stopped",
                                     "production_protection_interaction"],
        "analysis": "Exact event-time R/L integrals from run_interface.measurements; Monte Carlo SE of per-origin contrasts. Negative production effects retained.",
        "scope": "Conditional continuation of fixed-time generated preparations. Tests ongoing production versus inherited protective material; not autonomous origin, a maintained membrane, aging, functional memory, or matched-function efficiency."
    }
    encoded = json.dumps(protocol, indent=2) + "\n"
    if PROTOCOL.exists():
        assert PROTOCOL.read_text() == encoded, "refusing to replace a differing frozen maintenance protocol"
    else:
        PROTOCOL.write_text(encoded)
    return protocol


def reconstruct_preparation(item: dict, protocol: dict) -> tuple[object, dict]:
    path = ROOT / item["source_trajectory"]
    assert digest(path) == item["source_trajectory_sha256"]
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        original = json.load(stream)
    assert original["seed"] == item["origin_seed"]
    assert original["flags"] == asdict(ARMS[0])
    assert original["params"] == protocol["params"]
    assert original["horizon"] >= protocol["branch_time"]
    params = Params(**original["params"])
    model = Model(params, ARMS[0])
    state = model.state_from_dict(original["initial_state"])
    prefix = []
    for saved in original["ledger"]:
        if saved["time"] > protocol["branch_time"]:
            break
        state.time = saved["time"]
        replayed = model.apply(state, saved["event_index"])
        assert replayed == saved, "frozen original trajectory did not replay exactly"
        prefix.append(saved)
    state.time = protocol["branch_time"]
    prepared = model.state_to_dict(state)
    prefix_run = {"params": original["params"], "flags": original["flags"], "seed": original["seed"],
                  "horizon": protocol["branch_time"], "initial_state": original["initial_state"],
                  "final_state": prepared, "ledger": prefix}
    prefix_measurements = measurements(prefix_run)
    metadata = {**item, "preparation_time": protocol["branch_time"],
                "original_initial_state": original["initial_state"],
                "prepared_state_before_bookkeeping_reset": prepared,
                "preparation_prefix_events": len(prefix),
                "preparation_measurements": prefix_measurements,
                "prepared_R": int(np.isin(state.material, (1, 2)).sum()),
                "prepared_L": int((state.material == 3).sum()),
                "prepared_fuel": int(state.fuel.sum()), "prepared_H": model.energy(state),
                "preparation_heat_to_bath": state.cumulative_heat,
                "preparation_external_work": state.external_work}
    # Only observer bookkeeping is re-zeroed. No physical variable is reset.
    post = model.copy_state(state)
    post.time = 0.0
    post.cumulative_heat = 0.0
    post.external_work = 0.0
    post.event_count = 0
    assert np.array_equal(post.material, state.material) and np.array_equal(post.fuel, state.fuel)
    assert np.array_equal(post.instance_tags, state.instance_tags) and post.next_tag == state.next_tag
    return post, metadata


def paired_contrasts(rows: list[dict], fields: list[str]) -> list[dict]:
    by_origin: dict[int, dict[tuple[bool, bool], dict]] = {}
    for row in rows:
        by_origin.setdefault(row["origin_seed"], {})[(row["flags"]["produce_L"], row["flags"]["protect_R"])] = row
    contrasts = {
        "production_effect_when_protected": {(True, True): 1, (False, True): -1},
        "production_effect_without_protection": {(True, False): 1, (False, False): -1},
        "protection_effect_with_production": {(True, True): 1, (True, False): -1},
        "protection_effect_after_production_stopped": {(False, True): 1, (False, False): -1},
        "production_protection_interaction": {(True, True): 1, (True, False): -1,
                                              (False, True): -1, (False, False): 1}}
    result = []
    for name, weights in contrasts.items():
        for field in fields:
            values = [sum(weight * arms[arm][field] for arm, weight in weights.items())
                      for _, arms in sorted(by_origin.items())]
            result.append({"contrast": name, "field": field, **estimate(values),
                           "definition": [{"produce_L": arm[0], "protect_R": arm[1], "weight": weight}
                                          for arm, weight in weights.items()],
                           "uncertainty": "SE of per-origin differences; includes measured within-preparation covariance."})
    return result


def execute() -> dict:
    assert PROTOCOL.exists(), "freeze maintenance-protocol.json before execution"
    protocol_hash = digest(PROTOCOL)
    protocol = json.loads(PROTOCOL.read_text())
    before = source_hashes()
    assert before == protocol["source_hashes"], "source differs from frozen protocol"
    assert digest(ROOT / protocol["source_results"]) == protocol["source_results_sha256"]
    out = ROOT / "maintenance"
    out.mkdir(exist_ok=True)
    rows, preparation_rows = [], []
    start = time.monotonic()
    params = Params(**protocol["params"])
    for item in protocol["preparations"]:
        post, preparation = reconstruct_preparation(item, protocol)
        preparation_rows.append(preparation)
        initial_signature = None
        post_seed = item["origin_seed"] + protocol["post_seed_offset"]
        for flags in ARMS:
            label = f"origin-s{item['origin_seed']}-post-s{post_seed}-p{int(flags.produce_L)}b{int(flags.protect_R)}"
            path = out / (label + ".json.gz")
            assert not path.exists(), "refusing to replace an existing maintenance trajectory"
            run = simulate(params, post_seed, protocol["continuation_horizon"], flags,
                           sample_dt=protocol["sample_dt"], initial_state=post)
            signature = json.dumps(run["initial_state"], sort_keys=True)
            if initial_signature is None:
                initial_signature = signature
            assert signature == initial_signature, "physical continuation preparations differ across arms"
            run.update(origin_seed=item["origin_seed"], preparation=preparation,
                       source_hashes=before, maintenance_protocol_sha256=protocol_hash,
                       stage="maintenance")
            save_run(path, run)
            row = measurements(run)
            row.update(origin_seed=item["origin_seed"], post_seed=post_seed, case="prepared_at_fixed_t20",
                       trajectory_file=str(path.relative_to(ROOT)), trajectory_sha256=digest(path),
                       prepared_R=preparation["prepared_R"], prepared_L=preparation["prepared_L"],
                       prepared_fuel=preparation["prepared_fuel"])
            rows.append(row)
            print(json.dumps({"run": label, "R_int": round(row["integral_R"], 3),
                              "L_int": round(row["integral_L"], 3), "final_R": row["final_R"],
                              "fuel": row["net_fuel_consumed"], "events": row["events"],
                              "elapsed_s": round(time.monotonic() - start, 2)}), flush=True)
    assert len(rows) == protocol["replicates_per_arm"] * protocol["number_arms"]
    fields = ["integral_R", "integral_L", "final_R", "final_L", "paid_R_births", "paid_L_births",
              "passive_R_deaths", "passive_L_deaths", "net_fuel_consumed", "net_heat_to_bath",
              "first_episode_duration"]
    groups = [{"flags": asdict(flags), "statistics": {
        field: estimate([row[field] for row in rows if row["flags"] == asdict(flags)]) for field in fields}}
              for flags in ARMS]
    result = {"stage": "maintenance", "source_hashes": before, "protocol_file": PROTOCOL.name,
              "protocol_sha256": protocol_hash, "origin_results_sha256": protocol["source_results_sha256"],
              "params": protocol["params"], "replicates_per_arm": protocol["replicates_per_arm"],
              "branch_time": protocol["branch_time"], "horizon": protocol["continuation_horizon"],
              "origin_seeds": [item["origin_seed"] for item in protocol["preparations"]],
              "all_preparations_included": True, "groups": groups,
              "paired_contrasts": paired_contrasts(rows, fields), "preparations": preparation_rows,
              "runs": rows, "elapsed_seconds": time.monotonic() - start,
              "scope": protocol["scope"],
              "interpretation_limit": "Production effects after t20 compare ongoing synthesis with legacy L and its unchanged passive chemistry. No outcome is filtered, and a production-off advantage must be retained. Instance renewal is not aging; fuel/function ratios do not establish matched-function efficiency."}
    assert source_hashes() == before and digest(PROTOCOL) == protocol_hash, "frozen source or protocol changed"
    (ROOT / "maintenance-results.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze-protocol", action="store_true")
    arguments = parser.parse_args()
    if arguments.freeze_protocol:
        frozen = freeze_protocol()
        print(json.dumps({"protocol": str(PROTOCOL), "sha256": digest(PROTOCOL),
                          "preparations": len(frozen["preparations"]), "runs": 4 * len(frozen["preparations"])}))
    else:
        execute()
