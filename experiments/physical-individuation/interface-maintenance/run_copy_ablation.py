"""Previously specified paired copying-channel ablation of produced protection."""
from dataclasses import asdict
import gzip
from hashlib import sha256
import json
from pathlib import Path

from interface_model import Params, Flags, simulate
from run_interface import measurements

ROOT = Path(__file__).resolve().parent


def main():
    protocol = json.loads((ROOT/"ensemble-protocol.json").read_text())
    flags = Flags(copy_R=False)
    params = Params(**protocol["params"])
    plan = {"params": asdict(params), "flags": asdict(flags), "seeds": protocol["seeds"],
            "horizon": protocol["horizon"], "comparison": "Same-seed produced-protection arm of frozen baseline ensemble",
            "scope": "A paired copying-channel removal; rare passive thermal births remain."}
    (ROOT/"copy-ablation-protocol.json").write_text(json.dumps(plan, indent=2)+"\n")
    hashes = {name: sha256((ROOT/name).read_bytes()).hexdigest()
              for name in ("interface_model.py", "run_copy_ablation.py", "run_interface.py", "interface_readouts.py")}
    directory = ROOT/"copy-ablation"
    directory.mkdir(exist_ok=True)
    rows = []
    for seed in plan["seeds"]:
        run = simulate(params, seed, plan["horizon"], flags, sample_dt=2.)
        run["source_hashes"] = hashes
        path = directory/f"s{seed}.json.gz"
        with gzip.open(path, "wt") as f:
            json.dump(run, f, separators=(",", ":"))
        row = measurements(run)
        row["trajectory_file"] = str(path.relative_to(ROOT))
        rows.append(row)
    assert hashes == {name: sha256((ROOT/name).read_bytes()).hexdigest() for name in hashes}
    (ROOT/"copy-ablation-results.json").write_text(json.dumps({**plan, "source_hashes": hashes, "runs": rows}, indent=2)+"\n")
    print(json.dumps({"runs": len(rows), "events": sum(r["events"] for r in rows),
                      "mean_integrated_R": sum(r["integral_R"] for r in rows)/len(rows)}))


if __name__ == "__main__":
    main()
