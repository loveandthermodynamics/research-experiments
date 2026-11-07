"""Bounded prospective response-by-renewal experiment; no hidden rewards."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import gzip
from hashlib import sha256
import json
from pathlib import Path
import time

import numpy as np
from scipy.stats import t

from differentiated_model import ResponseParams, simulate_differentiated
from driven_model import Params, Flags, DriveParams
from lineage_readouts import measure

ROOT = Path(__file__).resolve().parent
PRIOR = ROOT.parent / "driven-spatial"
MODES = ("opposed", "aligned")
PRIMARY = ("response_effect_full_renewal", "response_by_production_interaction")
FIELDS = ("integral_R", "integral_R1", "integral_L", "both_types_present_time",
          "integral_founder_descendant_R1", "integral_focal_descendants", "integral_focal_R1",
          "focal_paid_mutant_appeared_by_deadline", "focal_established", "focal_survives_to_horizon",
          "focal_R1_uninterrupted_to_horizon", "final_focal_descendants", "final_focal_R1",
          "paid_R_births_after_founder_instance_clearance", "chemical_work", "heat_to_bath",
          "final_R", "final_R1", "final_L", "first_episode_duration")


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def write_once(path, value):
    encoded = json.dumps(value, indent=2, allow_nan=False) + "\n"
    if path.exists():
        assert path.read_text() == encoded, f"frozen file differs: {path}"
    else:
        path.write_text(encoded)


def source_hashes():
    files = [ROOT / name for name in ("differentiated_model.py", "model-spec.md",
             "run_differentiated.py", "lineage_readouts.py", "protocol-design.md")]
    files += [PRIOR / "driven_model.py", ROOT.parent / "interface-maintenance/interface_model.py"]
    return {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else "../"+str(p.relative_to(ROOT.parent)): digest(p)
            for p in files}


def base_params():
    values = load(PRIOR / "confirmation-protocol.json")["base_params"]
    return {**values, "copy_rate": .028, "initial_R1_fraction": 0.}


def freeze_pilot():
    protocol = dict(stage="pilot", source_hashes=source_hashes(), params=base_params(),
                    drive=asdict(DriveParams()), horizon=100., sample_dt=2.,
                    settings=[dict(asymmetry=epsilon, flip_rate=kappa, initial_gate=None)
                              for epsilon in (.5, .8) for kappa in (.025, .1)],
                    seeds=list(range(1050001, 1050009)),
                    confirmation_seeds=list(range(1060001, 1060129)),
                    response_modes=list(MODES), production_flags=[True, False],
                    focal_appearance_deadline=20., focal_establishment_size=4,
                    primary_outcome="integral_R", primary_contrasts=list(PRIMARY),
                    selection="Maximize minimum pilot mean/SE over both primary contrasts; tie by smaller asymmetry then flip_rate. Zero SE uses signed 1e12, or zero for zero mean.",
                    confirmation_rule="Both two-sided t(n-1) 95% lower bounds above zero in one fresh 128-seed confirmation. No pooling, early stopping, rescue scan or excluded failures.",
                    primary_bootstrap_sensitivity=dict(replicates=10000, seed=1080001,
                        rule="Fixed percentile paired bootstrap for both primary contrasts; sensitivity only, does not replace the declared t criterion."),
                    interpretation="Supplied finite reversible chemistry. Fixed copying variation, stationary exogenous gate and pure R0 founders. Appearance, ancestry expansion, functional benefit and individuality remain separate.")
    write_once(ROOT / "pilot-protocol.json", protocol)


def estimate(values):
    a = np.asarray(values, dtype=float)
    n = len(a)
    mean = float(a.mean())
    se = float(a.std(ddof=1)/np.sqrt(n)) if n > 1 else 0.
    critical = float(t.ppf(.975, n-1)) if n > 1 else 0.
    return dict(n=n, mean=mean, monte_carlo_standard_error=se,
                degrees_freedom=n-1, t_95_interval=[mean-critical*se, mean+critical*se])


def summarize(rows):
    groups, contrasts = [], []
    settings = sorted({(r["asymmetry"], r["flip_rate"]) for r in rows})
    for epsilon, kappa in settings:
        selected = [r for r in rows if (r["asymmetry"], r["flip_rate"]) == (epsilon, kappa)]
        by_seed = {}
        for r in selected:
            by_seed.setdefault(r["origin_seed"], {})[(r["response_mode"], r["produce_L"])] = r
        for mode in MODES:
            for produce in (True, False):
                arm = [r for r in selected if (r["response_mode"], r["produce_L"]) == (mode, produce)]
                groups.append(dict(asymmetry=epsilon, flip_rate=kappa, response_mode=mode,
                    produce_L=produce, statistics={field: estimate([r[field] for r in arm]) for field in FIELDS},
                    focal_appearances=sum(r["focal_paid_mutant_appeared_by_deadline"] for r in arm),
                    focal_establishments_over_both_alleles=sum(r["focal_established"] for r in arm),
                    focal_survivors=sum(r["focal_survives_to_horizon"] for r in arm),
                    final_zero_R=sum(r["final_R"] == 0 for r in arm)))
        for field in FIELDS:
            full = [r[("opposed", True)][field]-r[("aligned", True)][field] for _, r in sorted(by_seed.items())]
            interaction = [r[("opposed", True)][field]-r[("aligned", True)][field]
                           -r[("opposed", False)][field]+r[("aligned", False)][field]
                           for _, r in sorted(by_seed.items())]
            for name, values in zip(PRIMARY, (full, interaction)):
                contrasts.append(dict(asymmetry=epsilon, flip_rate=kappa, contrast=name, field=field, **estimate(values)))
    return groups, contrasts


def select_setting(pilot):
    scores = []
    primary = [r for r in pilot["paired_contrasts"] if r["field"] == "integral_R"]
    for setting in sorted({(r["asymmetry"], r["flip_rate"]) for r in primary}):
        rows = [r for r in primary if (r["asymmetry"], r["flip_rate"]) == setting]
        assert len(rows) == 2
        z = [r["mean"]/r["monte_carlo_standard_error"] if r["monte_carlo_standard_error"] > 0
             else (float(np.sign(r["mean"]))*1e12) for r in rows]
        scores.append(dict(asymmetry=setting[0], flip_rate=setting[1], score=min(z), component_scores=z))
    winner = sorted(scores, key=lambda r: (-r["score"], r["asymmetry"], r["flip_rate"]))[0]
    return {"asymmetry": winner["asymmetry"], "flip_rate": winner["flip_rate"], "initial_gate": None}, scores


def freeze_confirmation():
    pilot = load(ROOT / "pilot-results.json")
    original = load(ROOT / "pilot-protocol.json")
    assert pilot["source_hashes"] == source_hashes() == original["source_hashes"]
    selected, scores = select_setting(pilot)
    protocol = {**original, "stage": "confirmation", "settings": [selected],
                "seeds": original["confirmation_seeds"], "selection_scores": scores,
                "pilot_protocol_sha256": digest(ROOT / "pilot-protocol.json"),
                "pilot_results_sha256": digest(ROOT / "pilot-results.json")}
    write_once(ROOT / "confirmation-protocol.json", protocol)


def execute(stage):
    result_path = ROOT / f"{stage}-results.json"
    assert not result_path.exists(), "result is already frozen"
    protocol_path = ROOT / f"{stage}-protocol.json"
    protocol = load(protocol_path)
    assert protocol["source_hashes"] == source_hashes()
    if stage == "confirmation":
        assert digest(ROOT / "pilot-results.json") == protocol["pilot_results_sha256"]
        assert digest(ROOT / "pilot-protocol.json") == protocol["pilot_protocol_sha256"]
    begin = time.monotonic()
    directory = ROOT / stage
    directory.mkdir(exist_ok=True)
    rows = []
    for setting in protocol["settings"]:
        for seed in protocol["seeds"]:
            gates, initial = None, None
            for mode in protocol["response_modes"]:
                for produce in protocol["production_flags"]:
                    name = f"e{setting['asymmetry']:.2f}-k{setting['flip_rate']:.3f}-s{seed}-{mode}-p{int(produce)}.json.gz"
                    path = directory / name
                    assert not path.exists(), "physical record already exists"
                    run = simulate_differentiated(Params(**protocol["params"]), seed, protocol["horizon"],
                            flags=Flags(produce_L=produce, protect_R=True, copy_R=True),
                            drive=DriveParams(**protocol["drive"]), response=ResponseParams(**setting),
                            response_mode=mode, sample_dt=protocol["sample_dt"])
                    gate_trace = [(r["time"], r.get("gate_after", r.get("after"))) for r in run["ledger"] if r["kind"] == "gate_flip"]
                    if gates is None:
                        gates, initial = gate_trace, run["initial_state"]
                    else:
                        assert gates == gate_trace and initial == run["initial_state"]
                    with gzip.open(path, "wt", encoding="utf-8") as stream:
                        json.dump(run, stream, separators=(",", ":"), allow_nan=False)
                    observations = measure(run, protocol["focal_appearance_deadline"], protocol["focal_establishment_size"])
                    rows.append(dict(origin_seed=seed, asymmetry=setting["asymmetry"], flip_rate=setting["flip_rate"],
                                     response_mode=mode, produce_L=produce,
                                     trajectory_file=str(path.relative_to(ROOT)), trajectory_sha256=digest(path), **observations))
            print(f"{stage}: e={setting['asymmetry']} k={setting['flip_rate']} seed={seed} complete", flush=True)
    groups, contrasts = summarize(rows)
    primary = [r for r in contrasts if r["field"] == "integral_R"]
    bootstrap = []
    if stage == "confirmation":
        paired = {}
        for row in rows:
            paired.setdefault(row["origin_seed"], {})[(row["response_mode"], row["produce_L"])] = row["integral_R"]
        vectors = np.array([[p[("opposed", True)]-p[("aligned", True)],
                            p[("opposed", True)]-p[("aligned", True)]-p[("opposed", False)]+p[("aligned", False)]]
                           for _, p in sorted(paired.items())])
        spec = protocol["primary_bootstrap_sensitivity"]
        rng = np.random.default_rng(spec["seed"])
        indices = rng.integers(0, len(vectors), size=(spec["replicates"], len(vectors)))
        means = vectors[indices].mean(axis=1)
        bootstrap = [dict(contrast=name, replicates=spec["replicates"], seed=spec["seed"],
                          percentile_95_interval=np.quantile(means[:, i], [.025, .975]).tolist())
                     for i, name in enumerate(PRIMARY)]
    result = dict(stage=stage, source_hashes=source_hashes(), protocol_file=protocol_path.name,
                  protocol_sha256=digest(protocol_path), all_preparations_included=True,
                  runs=rows, groups=groups, paired_contrasts=contrasts,
                  primary_results=primary,
                  primary_bootstrap_sensitivity=bootstrap,
                  declared_joint_claim_supported=(all(r["t_95_interval"][0] > 0 for r in primary) if stage == "confirmation" else None),
                  elapsed_seconds=time.monotonic()-begin,
                  scope="Finite end-to-end response-by-production comparison under supplied reversible chemistry. Ancestry follows paid copying without immigration. No individual, rare-type Lyapunov exponent, efficiency optimum, aging or universal biodiversity law is inferred.")
    assert result["source_hashes"] == protocol["source_hashes"]
    write_once(result_path, result)
    print(json.dumps({"stage": stage, "primary_results": primary,
                      "declared_joint_claim_supported": result["declared_joint_claim_supported"]}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("freeze-pilot", "pilot", "freeze-confirmation", "confirmation"))
    action = parser.parse_args().action
    if action == "freeze-pilot":
        freeze_pilot()
    elif action == "freeze-confirmation":
        freeze_confirmation()
    else:
        execute(action)
