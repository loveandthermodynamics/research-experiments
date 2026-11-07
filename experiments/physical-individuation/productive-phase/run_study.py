"""Execute only the frozen, bounded joined productive-phase comparison."""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
import gzip
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
from productive_model import Flags, Params, simulate

SEEDS = list(range(730100, 730108))
ARMS = ("specialist", "generalist", "X_only", "Y_only")
SOURCES = ("productive_model.py", "model-spec.md", "balance-and-cost.md", "generalist-comparison.md",
           "study-design.md", "spatial_readouts.py", "dependencies/morphology.py",
           "renewal_readouts.py", "run_study.py", "analyze_study.py")


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    path = Path(path)
    temporary = path.with_name(path.name+".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+"\n")
    temporary.replace(path)


def configuration(side, arm):
    p = Params(side=side, material_count=77*side*side//256, initial_X=8*side*side//256, initial_Y=0)
    if arm == "specialist":
        return p, Flags(), "specialist"
    if arm == "generalist":
        return p, Flags(), "generalist"
    if arm == "X_only":
        return p, Flags(enable_Y=False), "specialist"
    if arm == "Y_only":
        return p, Flags(enable_X=False), "specialist"
    raise ValueError(arm)


def freeze():
    target = ROOT/"protocol.json"
    if target.exists():
        raise FileExistsError("protocol already frozen")
    audit_path = ROOT/"preflight-verification.json"
    audit = json.loads(audit_path.read_text())
    if not audit.get("checks_passed"):
        raise RuntimeError("independent preflight must pass before outcomes")
    # Report schemas may differ, but a present model binding must match exactly.
    if "model_sha256" in audit:
        assert audit["model_sha256"] == digest(ROOT/"productive_model.py")
    configs = {}
    for side in (16, 32):
        for arm in ARMS:
            p, flags, mode = configuration(side, arm)
            configs[f"{side}-{arm}"] = dict(params=asdict(p), flags=asdict(flags), activity_mode=mode)
    dump(target, dict(study="Joined productive copying and molecular phase: bounded prospective test",
         seeds=SEEDS, sides=[16, 32], arms=list(ARMS), horizon=400., sample_dt=5.,
         early_window=[100., 200.], late_window=[300., 400.], configurations=configs,
         source_sha256={n:digest(ROOT/n) for n in SOURCES},
         independent_preflight_sha256=digest(audit_path), no_outcome_parameter_selection=True,
         no_rescue_or_optional_extension=True, required_claim_A_gates=20, required_claim_B_gates=2))
    print(json.dumps(dict(frozen=True, protocol_sha256=digest(target))), flush=True)


def one(job):
    side, arm, seed, protocol_hash = job
    path = ROOT/"records"/f"L{side}-{arm}-seed{seed}.json.gz"
    p, flags, mode = configuration(side, arm)
    start = time.monotonic()
    if path.exists():
        with gzip.open(path, "rt") as stream:
            result = json.load(stream)
        assert result["protocol_sha256"] == protocol_hash and result["seed"] == seed
        assert result["params"] == asdict(p) and result["flags"] == asdict(flags)
        assert result["activity_mode"] == mode and result["horizon"] == 400.
        reused = True
    else:
        result = simulate(p, seed, 400., flags=flags, activity_mode=mode, sample_dt=5.)
        result.update(arm=arm, protocol_sha256=protocol_hash)
        temporary = path.with_name(path.name+".tmp")
        with gzip.open(temporary, "wt", compresslevel=6) as stream:
            json.dump(result, stream, separators=(",", ":"), allow_nan=False)
        temporary.replace(path)
        reused = False
    return dict(side=side, arm=arm, seed=seed, path=str(path.relative_to(ROOT)), sha256=digest(path),
                bytes=path.stat().st_size, events=len(result["ledger"]), samples=len(result["samples"]),
                proposals=result["proposal_statistics"]["proposals"], final=result["final_summary"],
                first_law_error=result["first_law_error"], seconds=time.monotonic()-start, reused=reused)


def run(workers=2):
    protocol_path = ROOT/"protocol.json"
    protocol = json.loads(protocol_path.read_text())
    for name, expected in protocol["source_sha256"].items():
        assert digest(ROOT/name) == expected, ("changed frozen source", name)
    assert digest(ROOT/"preflight-verification.json") == protocol["independent_preflight_sha256"]
    protocol_hash = digest(protocol_path)
    (ROOT/"records").mkdir(exist_ok=True)
    jobs = [(side, arm, seed, protocol_hash) for side in protocol["sides"]
            for arm in protocol["arms"] for seed in protocol["seeds"]]
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(one, job):job for job in jobs}
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            rows.sort(key=lambda r:(r["side"], r["arm"], r["seed"]))
            dump(ROOT/"execution-progress.json", dict(expected_records=len(jobs), completed_records=len(rows),
                 protocol_sha256=protocol_hash, records=rows))
            print(json.dumps(dict(completed=len(rows), total=len(jobs), side=row["side"], arm=row["arm"],
                                  seed=row["seed"], events=row["events"], seconds=round(row["seconds"], 3))), flush=True)
    dump(ROOT/"execution.json", dict(passed_execution=True, protocol_sha256=protocol_hash,
         unique_records=len(rows), events=sum(r["events"] for r in rows),
         samples=sum(r["samples"] for r in rows), all_preparations_retained=True, records=rows))


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--freeze":
        freeze()
    elif len(sys.argv) in (2, 3) and sys.argv[1] == "--run":
        run(int(sys.argv[2]) if len(sys.argv) == 3 else 2)
    else:
        raise SystemExit("use --freeze after independent preflight, or --run [workers]")
