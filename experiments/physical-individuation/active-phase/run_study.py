"""Execute only the prospectively specified bounded active-phase study."""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, replace
import gzip
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
from active_phase_model import Params, simulate

SOURCES = ("active_phase_model.py", "model-spec.md", "morphology.py", "study-design.md", "run_study.py", "analyze_study.py")
SEEDS = list(range(720100, 720108))
ARMS = ("full", "no_drive", "no_affinity")


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def dump(path: Path, value: object) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def params_for(side: int, arm: str) -> Params:
    base = Params(side=side, material_count=77 * side * side // 256,
                  enzyme_count=8 * side * side // 256)
    if arm == "no_drive":
        return replace(base, mu=0.)
    if arm == "no_affinity":
        return replace(base, J_EB=0.)
    if arm != "full":
        raise ValueError(arm)
    return base


def freeze() -> None:
    target = ROOT / "protocol.json"
    if target.exists():
        raise FileExistsError("frozen protocol already exists")
    audit_path = ROOT / "preflight-verification.json"
    audit = json.loads(audit_path.read_text())
    if not audit.get("checks_passed"):
        raise RuntimeError("independent physical preflight must pass first")
    assert audit["model_sha256"] == digest(ROOT / "active_phase_model.py")
    assert audit["morphology_sha256"] == digest(ROOT / "morphology.py")
    assert audit["verifier_sha256"] == digest(ROOT / "verify_active_phase.py")
    dump(target, {"study": "bounded active phase, prospective operational arrest test",
                  "seeds": SEEDS, "sides": [16, 32], "arms": list(ARMS),
                  "horizon": 400., "sample_dt": 5.,
                  "early_window": [100., 200.], "late_window": [300., 400.],
                  "parameters": {f"{side}-{arm}": asdict(params_for(side, arm))
                                 for side in (16, 32) for arm in ARMS},
                  "source_sha256": {name: digest(ROOT / name) for name in SOURCES},
                  "independent_preflight_sha256": digest(audit_path),
                  "no_outcome_parameter_selection": True,
                  "no_rescue_or_optional_extension": True})
    print(json.dumps({"frozen": True, "protocol_sha256": digest(target)}), flush=True)


def one(job: tuple[int, str, int, str]) -> dict:
    side, arm, seed, protocol_sha = job
    destination = ROOT / "records" / f"L{side}-{arm}-seed{seed}.json.gz"
    start = time.monotonic()
    if destination.exists():
        with gzip.open(destination, "rt") as stream:
            result = json.load(stream)
        assert result["protocol_sha256"] == protocol_sha
        assert result["params"] == asdict(params_for(side, arm))
        assert result["seed"] == seed and result["horizon"] == 400.
        reused = True
    else:
        result = simulate(params_for(side, arm), seed, 400., sample_dt=5.)
        result["arm"] = arm
        result["protocol_sha256"] = protocol_sha
        temporary = destination.with_name(destination.name + ".tmp")
        with gzip.open(temporary, "wt", compresslevel=6) as stream:
            json.dump(result, stream, separators=(",", ":"), allow_nan=False)
        temporary.replace(destination)
        reused = False
    return {"side": side, "arm": arm, "seed": seed,
            "path": str(destination.relative_to(ROOT)), "sha256": digest(destination),
            "bytes": destination.stat().st_size, "events": len(result["ledger"]),
            "samples": len(result["samples"]),
            "proposals": result["proposal_statistics"]["proposals"],
            "final": result["final_summary"], "first_law_error": result["first_law_error"],
            "seconds": time.monotonic() - start, "reused": reused}


def run(workers: int = 2) -> None:
    protocol_path = ROOT / "protocol.json"
    protocol = json.loads(protocol_path.read_text())
    for name, expected in protocol["source_sha256"].items():
        assert digest(ROOT / name) == expected, ("changed frozen source", name)
    protocol_sha = digest(protocol_path)
    (ROOT / "records").mkdir(exist_ok=True)
    jobs = [(side, arm, seed, protocol_sha) for side in protocol["sides"]
            for arm in protocol["arms"] for seed in protocol["seeds"]]
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(one, job): job for job in jobs}
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            rows.sort(key=lambda r: (r["side"], r["arm"], r["seed"]))
            dump(ROOT / "execution-progress.json", {"expected_records": len(jobs),
                 "completed_records": len(rows), "protocol_sha256": protocol_sha, "records": rows})
            print(json.dumps({"completed": len(rows), "total": len(jobs),
                              "side": row["side"], "arm": row["arm"], "seed": row["seed"],
                              "events": row["events"], "seconds": round(row["seconds"], 3)}), flush=True)
    dump(ROOT / "execution.json", {"passed_execution": True, "protocol_sha256": protocol_sha,
         "records": rows, "unique_records": len(rows), "events": sum(r["events"] for r in rows),
         "samples": sum(r["samples"] for r in rows), "all_preparations_retained": True})


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--freeze":
        freeze()
    elif len(sys.argv) in (2, 3) and sys.argv[1] == "--run":
        run(int(sys.argv[2]) if len(sys.argv) == 3 else 2)
    else:
        raise SystemExit("use --freeze (after independent preflight) or --run [workers]")
