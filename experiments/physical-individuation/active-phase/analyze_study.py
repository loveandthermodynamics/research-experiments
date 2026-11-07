"""Seed-level analysis of the prospectively frozen active-phase protocol."""
from __future__ import annotations

from collections import Counter
import gzip
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import t

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
from morphology import observe

METRICS = ("B_fraction", "largest_B_fraction", "mass_weighted_cluster_size",
           "excess_neighbor_B_fraction", "E_on_B_excess_fraction", "E_enrichment",
           "spectral_wavelength_low_k", "low_k_power_fraction", "clusters_ge_9")


def average(values):
    return None if any(v is None for v in values) else float(np.mean(values))


def ratio(a, b):
    return None if a is None or b is None or b == 0 else a / b


def difference(a, b):
    return None if a is None or b is None else a - b


def estimate(values):
    if any(v is None or not np.isfinite(v) for v in values):
        return {"estimable": False, "n": len(values), "values": values,
                "reason": "undefined required values retained; no seed filtering"}
    array = np.asarray(values, dtype=float)
    mean, sem = float(array.mean()), float(array.std(ddof=1) / np.sqrt(len(array)))
    margin = float(t.ppf(.975, len(array) - 1)) * sem
    # Same deterministic row resamples for every metric, preserving seed pairing.
    rng = np.random.default_rng(384701)
    rows = rng.integers(len(array), size=(10000, len(array)))
    bootstrap = np.mean(array[rows], axis=1)
    return {"estimable": True, "n": len(values), "values": values, "mean": mean,
            "sem": sem, "t95": [mean - margin, mean + margin],
            "paired_seed_bootstrap95": list(map(float, np.quantile(bootstrap, [.025, .975])))}


def gate(name, values, relation, threshold):
    result = estimate(values)
    result.update(name=name, relation=relation, threshold=threshold)
    if not result["estimable"]:
        result["met"] = False
    elif relation == "lower_gt":
        result["met"] = result["t95"][0] > threshold
    elif relation == "upper_le":
        result["met"] = result["t95"][1] <= threshold
    else:
        raise ValueError(relation)
    return result


def read_one(path):
    with gzip.open(path, "rt") as stream:
        physical = json.load(stream)
    side = physical["params"]["side"]
    rows = []
    for sample in physical["samples"]:
        observation = observe(np.asarray(sample["material"]).reshape(side, side),
                              np.asarray(sample["enzyme"]).reshape(side, side))
        observation["large_B_mass_fraction"] = (observation["B_mass_in_clusters_ge_9"]
                                                / observation["n_B"] if observation["n_B"] else 0.)
        observation["time"] = sample["time"]
        rows.append(observation)
    windows = {}
    for name, bounds in (("early", (100., 200.)), ("late", (300., 400.))):
        selected = [r for r in rows if bounds[0] <= r["time"] <= bounds[1]]
        assert len(selected) == 21
        windows[name] = {metric: average([r[metric] for r in selected])
                         for metric in (*METRICS, "large_B_mass_fraction")}
    growth = ratio(windows["late"]["spectral_wavelength_low_k"],
                   windows["early"]["spectral_wavelength_low_k"])
    cohort_sample = next(s for s in physical["samples"] if s["time"] == 200.)
    cohort = {tag for a, tag in zip(cohort_sample["material"], cohort_sample["material_tags"]) if a == 2}
    converted, conversions, late_conversions = set(), Counter(), Counter()
    for event in physical["ledger"]:
        if event["kind"] not in ("passive", "catalyzed"):
            continue
        before, after = event["before"]["material"][0], event["after"]["material"][0]
        label = f"{event['kind']}_{before}_to_{after}"
        if event["time"] > 200.:
            conversions[label] += 1
            if before == 2:
                tag = event["before"]["material_tags"][0]
                if tag in cohort:
                    converted.add(tag)
        if event["time"] > 300.:
            late_conversions[label] += 1
    final_b_tags = {tag for a, tag in zip(physical["final_state"]["material"],
                                         physical["final_state"]["material_tags"]) if a == 2}
    turnover = {"B_carrier_cohort_at200": len(cohort),
                "cohort_tags_converted_to_A_after200": len(converted),
                "cohort_tags_B_at400_anywhere": len(cohort & final_b_tags),
                "cohort_tags_never_converted_to_A_after200": len(cohort - converted),
                "conversions_after200": dict(conversions), "conversions_after300": dict(late_conversions),
                "scope": "Conserved carrier tracers and chemical-state cycling, without an organizational identity or heritable template interpretation."}
    return {"side": side, "arm": physical["arm"], "seed": physical["seed"],
            "path": str(path.relative_to(ROOT)), "record_sha256": sha256(path.read_bytes()).hexdigest(),
            "windows": windows, "growth": growth, "turnover": turnover,
            "timeseries": rows, "physical_final": physical["final_summary"]}


def main():
    protocol = json.loads((ROOT / "protocol.json").read_text())
    execution = json.loads((ROOT / "execution.json").read_text())
    rows = [read_one(ROOT / record["path"]) for record in execution["records"]]
    index = {(r["side"], r["arm"], r["seed"]): r for r in rows}
    seeds = protocol["seeds"]
    get = lambda side, arm, seed: index[(side, arm, seed)]
    gates = []
    for side in (16, 32):
        for metric, threshold in (("B_fraction", .05), ("large_B_mass_fraction", .5),
                                  ("excess_neighbor_B_fraction", .15)):
            gates.append(gate(f"L{side}_late_{metric}",
                              [get(side, "full", seed)["windows"]["late"][metric] for seed in seeds],
                              "lower_gt", threshold))
        gates.append(gate(f"L{side}_full_wavelength_growth",
                          [get(side, "full", seed)["growth"] for seed in seeds], "upper_le", 1.25))
        gates.append(gate(f"L{side}_no_drive_minus_full_growth",
                          [difference(get(side, "no_drive", seed)["growth"],
                                      get(side, "full", seed)["growth"]) for seed in seeds], "lower_gt", 0.))
    gates.append(gate("full_late_wavelength_L32_over_L16",
                      [ratio(get(32, "full", seed)["windows"]["late"]["spectral_wavelength_low_k"],
                             get(16, "full", seed)["windows"]["late"]["spectral_wavelength_low_k"])
                       for seed in seeds], "upper_le", 1.5))
    gates.append(gate("full_late_largest_B_fraction_L16_minus_L32",
                      [difference(get(16, "full", seed)["windows"]["late"]["largest_B_fraction"],
                                  get(32, "full", seed)["windows"]["late"]["largest_B_fraction"])
                       for seed in seeds], "lower_gt", 0.))
    summary = {}
    for side in (16, 32):
        for arm in protocol["arms"]:
            entries = [get(side, arm, seed) for seed in seeds]
            summary[f"L{side}-{arm}"] = {
                "late": {metric: estimate([r["windows"]["late"][metric] for r in entries])
                         for metric in (*METRICS, "large_B_mass_fraction")},
                "growth": estimate([r["growth"] for r in entries]),
                "chemical_work_per_area": estimate([r["physical_final"]["chemical_work"] / side**2 for r in entries]),
                "heat_to_bath_per_area": estimate([r["physical_final"]["heat_to_bath"] / side**2 for r in entries]),
                "turnover": [r["turnover"] for r in entries]}
    result = {"protocol_sha256": sha256((ROOT / "protocol.json").read_bytes()).hexdigest(),
              "analysis_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "unique_records": len(rows), "independent_seed_rows": len(seeds),
              "joint_claim_met": all(g["met"] for g in gates), "primary_gates": gates,
              "summary": summary, "rows": rows,
              "scope": "Finite horizon400, two finite domains, supplied chemistry/reservoirs, analysis-only tolerances. No asymptotic phase proof, copied biodiversity, self-renewed enzyme, aging, recursive individual or moral conclusion."}
    (ROOT / "analysis.json").write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"joint_claim_met": result["joint_claim_met"],
                      "gates": [{k: g.get(k) for k in ("name", "estimable", "mean", "t95", "met")}
                                for g in gates]}, indent=2), flush=True)


if __name__ == "__main__":
    main()
