"""Response variation in the SAME paid protective chemistry under thermal gates."""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import gzip
from hashlib import sha256
import json
from pathlib import Path
import time

import numpy as np

from autonomous_environment import GateParams, EnvironmentModel, simulate_environment
from interface_model import Params, Flags
from interface_readouts import estimate
from run_interface import measurements

ROOT = Path(__file__).resolve().parent


def hashes():
    return {name: sha256((ROOT/name).read_bytes()).hexdigest() for name in
            ("interface_model.py", "autonomous_environment.py", "run_responses.py", "run_interface.py", "interface_readouts.py")}


def readout(run):
    base = dict(run)
    base["initial_state"] = run["initial_state"]["base_state"]
    base["final_state"] = run["final_state"]["base_state"]
    base["ledger"] = [e for e in run["ledger"] if e["kind"] != "environment_gate"]
    row = measurements(base)
    first = next((e for e in base["ledger"] if e["after"] == 2), None)
    row["first_variant_channel"] = None if first is None else first["kind"]
    row["first_variant_time"] = None if first is None else first["time"]
    row["variant_paid_births"] = sum(e["kind"] == "copy_R" and e["after"] == 2 for e in base["ledger"])
    row["variant_thermal_births"] = sum(e["kind"] == "passive" and e["after"] == 2 for e in base["ledger"])
    row["total_gate_heat"] = run["final_state"]["cumulative_gate_heat"]
    row["joint_heat"] = run["final_summary"]["heat_to_bath"]
    model = EnvironmentModel(Params(**run["params"]), Flags(**run["flags"]),
                             GateParams(**run["gateparams"]), run["response_mode"])
    h0 = model.energy(model.state_from_dict(run["initial_state"]))
    h1 = model.energy(model.state_from_dict(run["final_state"]))
    row["joint_first_law_error"] = abs(h0-h1-row["joint_heat"])
    assert row["joint_first_law_error"] < 1e-8
    return row


def execute(stage, replicates):
    before = hashes()
    base_params = Params(**json.loads((ROOT/"ensemble-protocol.json").read_text())["params"])
    if stage == "response-pilot":
        seeds = list(range(2201, 2201+replicates))
        regimes = {"moderate": GateParams(bad_multiplier=(12., 12.)),
                   "severe": GateParams(bad_multiplier=(40., 40.)),
                   "extreme": GateParams(bad_multiplier=(100., 100.))}
        horizon = 120.
    else:
        protocol = json.loads((ROOT/"response-protocol.json").read_text())
        assert protocol["replicates_per_arm"] == replicates
        seeds = protocol["seeds"]
        regimes = {"selected": GateParams(**protocol["gateparams"])}
        horizon = protocol["horizon"]
    out = ROOT/stage
    out.mkdir(exist_ok=True)
    rows = []
    started = time.monotonic()
    for name, gp in regimes.items():
        assert gp.p_bad[0] == gp.p_bad[1] and gp.gamma[0] == gp.gamma[1]
        for seed in seeds:
            history = None
            for mode in ("common", "differentiated"):
                for varied in (False, True):
                    mu = .2 if varied else 0.
                    params = replace(base_params, mutation=((1-mu, mu), (mu, 1-mu)))
                    run = simulate_environment(params, seed, horizon, Flags(), gp, mode, sample_dt=2.)
                    run["source_hashes"] = before
                    gates = [(e["time"], e["i"], e["after"]) for e in run["ledger"] if e["kind"] == "environment_gate"]
                    if history is None:
                        history = gates
                    assert gates == history, "matched gate histories differ"
                    label = f"{name}-s{seed}-{mode}-v{int(varied)}"
                    path = out/(label+".json.gz")
                    with gzip.open(path, "wt") as f:
                        json.dump(run, f, separators=(",", ":"))
                    row = readout(run)
                    row.update(case=name, mode=mode, varied=varied, trajectory_file=str(path.relative_to(ROOT)))
                    rows.append(row)
            print(json.dumps({"stage": stage, "case": name, "seed": seed,
                              "elapsed_s": round(time.monotonic()-started, 2)}), flush=True)
    fields = ("integral_R", "integral_R1", "net_fuel_consumed", "first_episode_duration", "joint_heat")
    groups, interactions = [], []
    for name in regimes:
        for mode in ("common", "differentiated"):
            for varied in (False, True):
                selected = [r for r in rows if (r["case"], r["mode"], r["varied"]) == (name, mode, varied)]
                groups.append({"case": name, "mode": mode, "varied": varied,
                    "statistics": {field: estimate([r[field] for r in selected]) for field in fields},
                    "first_variant_paid_copy_fraction": sum(r["first_variant_channel"] == "copy_R" for r in selected)/len(selected),
                    "first_variant_thermal_fraction": sum(r["first_variant_channel"] == "passive" for r in selected)/len(selected)})
        paired = []
        for seed in seeds:
            a = {(r["mode"], r["varied"]): r["integral_R"] for r in rows if r["case"] == name and r["seed"] == seed}
            paired.append(a[("differentiated", True)]-a[("differentiated", False)]-
                          a[("common", True)]+a[("common", False)])
        interactions.append({"case": name, "paired_interaction_response_by_variation": estimate(paired)})
    assert before == hashes(), "source changed during comparison"
    summary = {"stage": stage, "source_hashes": before, "base_params": asdict(base_params),
               "gate_regimes": {k: asdict(v) for k, v in regimes.items()}, "replicates_per_arm": replicates,
               "seeds": seeds, "horizon": horizon, "groups": groups, "interactions": interactions,
               "runs": rows, "elapsed_seconds": time.monotonic()-started,
               "scope": "Two supplied response capabilities with paid variation and matched autonomous gate histories. No unlimited diversity, recursive individuality, biological aging or moral legitimacy follows."}
    (ROOT/(stage+"-results.json")).write_text(json.dumps(summary, indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("response-pilot", "response-ensemble"))
    parser.add_argument("--replicates", type=int, required=True)
    args = parser.parse_args()
    execute(args.stage, args.replicates)
