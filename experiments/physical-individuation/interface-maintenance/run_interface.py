"""Declared exploratory pilot and independently seeded factorial ensemble.

Uses the physical model only for simulation; measurements never affect rates.
All trajectories, including unfavorable ones, are retained for event replay.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
from hashlib import sha256
import gzip
import json
from pathlib import Path
import time

import numpy as np

from interface_model import Params, Flags, Model, simulate
from interface_readouts import spatial_readouts, estimate, difference_in_differences

ROOT = Path(__file__).resolve().parent
ARMS = [Flags(produce_L=p, protect_R=b) for p, b in
        ((True, True), (True, False), (False, True), (False, False))]


def source_hashes():
    return {name: sha256((ROOT/name).read_bytes()).hexdigest()
            for name in ("interface_model.py", "interface_readouts.py", "run_interface.py")}


def measurements(run):
    initial, final = run["initial_state"], run["final_state"]
    material = np.asarray(initial["material"])
    r = int(np.sum((material == 1) | (material == 2)))
    l = int(np.sum(material == 3))
    r0, r1 = int(np.sum(material == 1)), int(np.sum(material == 2))
    integrals = dict(R=0.0, L=0.0, R0=0.0, R1=0.0)
    previous = 0.0
    counts = {"paid_R_births": 0, "reverse_R_copies": 0,
              "paid_L_births": 0, "reverse_L_production": 0,
              "passive_R_births": 0, "passive_R_deaths": 0,
              "passive_L_births": 0, "passive_L_deaths": 0}
    born_paid_r, born_paid_l = set(), set()
    first_R_extinction = None
    for event in run["ledger"]:
        dt = event["time"]-previous
        assert dt >= 0
        for field, value in (("R", r), ("L", l), ("R0", r0), ("R1", r1)):
            integrals[field] += value*dt
        previous = event["time"]
        if event["kind"] in ("material_swap", "fuel_swap"):
            continue
        before, after = event["before"], event["after"]
        delta_r = int(after in (1, 2))-int(before in (1, 2))
        r += delta_r
        r0 += int(after == 1)-int(before == 1)
        r1 += int(after == 2)-int(before == 2)
        l += int(after == 3)-int(before == 3)
        if event["kind"] == "copy_R":
            key = "paid_R_births" if after in (1, 2) else "reverse_R_copies"
            if after in (1, 2):
                born_paid_r.add(event["instance_created"])
        elif event["kind"] == "produce_L":
            key = "paid_L_births" if after == 3 else "reverse_L_production"
            if after == 3:
                born_paid_l.add(event["instance_created"])
        elif before in (1, 2) or after in (1, 2):
            key = "passive_R_births" if after in (1, 2) else "passive_R_deaths"
        else:
            key = "passive_L_births" if after == 3 else "passive_L_deaths"
        counts[key] += 1
        if r == 0 and first_R_extinction is None:
            first_R_extinction = event["time"]
    for field, value in (("R", r), ("L", l), ("R0", r0), ("R1", r1)):
        integrals[field] += value*(run["horizon"]-previous)
    f = np.asarray(final["material"])
    tags = np.asarray(final["instance_tags"])
    is_r = (f == 1) | (f == 2)
    is_l = f == 3
    initial_tags = set(t for t in initial["instance_tags"] if t)
    founder_survivors = sum(int(t) in initial_tags for t in tags[is_r])
    net_fuel = sum(initial["fuel"])-sum(final["fuel"])
    paid_net = counts["paid_R_births"]+counts["paid_L_births"]-counts["reverse_R_copies"]-counts["reverse_L_production"]
    assert net_fuel == paid_net
    assert r == int(is_r.sum()) and l == int(is_l.sum())
    model = Model(Params(**run["params"]), Flags(**run["flags"]))
    h0 = model.energy(model.state_from_dict(initial))
    h1 = model.energy(model.state_from_dict(final))
    error = abs(h0-h1-final["cumulative_heat"])
    assert error < 1e-8
    return {"seed": run["seed"], "flags": run["flags"], **counts,
            **{"integral_"+key: value for key, value in integrals.items()},
            "final_R": r, "final_L": l, "founder_R_remaining": founder_survivors,
            "final_paid_R_instances": sum(int(t) in born_paid_r for t in tags[is_r]),
            "final_paid_L_instances": sum(int(t) in born_paid_l for t in tags[is_l]),
            "first_R_extinction": first_R_extinction,
            "first_episode_duration": run["horizon"] if first_R_extinction is None else first_R_extinction,
            "net_fuel_consumed": net_fuel, "fuel_energy_transferred": net_fuel*run["params"]["E_F"],
            "net_heat_to_bath": final["cumulative_heat"], "initial_H": h0, "final_H": h1,
            "first_law_error": error, "events": len(run["ledger"]),
            "fuel_energy_per_integral_R": net_fuel*run["params"]["E_F"]/integrals["R"] if integrals["R"] else None,
            "final_spatial": spatial_readouts(f, run["params"]["side"])}


def save_run(path, run):
    with gzip.open(path, "wt", encoding="utf-8") as f:
        json.dump(run, f, separators=(",", ":"))


def execute(stage, count):
    before_hashes = source_hashes()
    out = ROOT/stage
    out.mkdir(exist_ok=True)
    balanced = Params(side=8, theta=0.7, epsilon_R=8, epsilon_L=8, E_F=12,
                      J_LL=1.0, J_RL=1.0, copy_rate=0.06, production_rate=0.02,
                      passive_R=0.12, passive_L=0.04, protection_strength=1.0,
                      transport_protection_strength=0.25, initial_templates=8)
    if stage == "pilot":
        cases = {"balanced": balanced,
                 "high_production": replace(balanced, production_rate=0.08),
                 "weak_protection": replace(balanced, protection_strength=0.2,
                                             transport_protection_strength=0.05)}
        seeds = list(range(1101, 1101+count))
        horizon = 120.0
    else:
        protocol = json.loads((ROOT/"ensemble-protocol.json").read_text())
        assert protocol["stage"] == stage and protocol["replicates_per_arm"] == count
        cases = {"selected": Params(**protocol["params"])}
        seeds = protocol["seeds"]
        horizon = protocol["horizon"]
    rows = []
    start = time.monotonic()
    for name, params in cases.items():
        for seed in seeds:
            for flags in ARMS:
                label = f"{name}-s{seed}-p{int(flags.produce_L)}b{int(flags.protect_R)}"
                path = out/(label+".json.gz")
                run = simulate(params, seed, horizon, flags, sample_dt=2.0)
                run["source_hashes"] = before_hashes
                save_run(path, run)
                row = measurements(run)
                row.update(case=name, trajectory_file=str(path.relative_to(ROOT)))
                rows.append(row)
                print(json.dumps({"stage": stage, "run": label, "R_int": round(row["integral_R"], 3),
                                  "R_final": row["final_R"], "L_final": row["final_L"],
                                  "fuel": row["net_fuel_consumed"], "events": row["events"],
                                  "elapsed_s": round(time.monotonic()-start, 2)}), flush=True)
    fields = ["integral_R", "integral_L", "final_R", "final_L", "founder_R_remaining",
              "final_paid_R_instances", "final_paid_L_instances", "first_episode_duration",
              "net_fuel_consumed", "net_heat_to_bath", "fuel_energy_per_integral_R"]
    groups = []
    for name in cases:
        for flags in ARMS:
            arm = [r for r in rows if r["case"] == name and r["flags"] == asdict(flags)]
            groups.append({"case": name, "flags": asdict(flags),
                           "statistics": {field: estimate([r[field] for r in arm]) for field in fields}})
    interaction = [{"case": name, **difference_in_differences([r for r in rows if r["case"] == name], "integral_R")}
                   for name in cases]
    assert source_hashes() == before_hashes, "source changed during ensemble; results not a frozen-source comparison"
    result = {"stage": stage, "source_hashes": before_hashes, "params_by_case": {k: asdict(v) for k, v in cases.items()},
              "horizon": horizon, "replicates_per_arm": count, "seeds": seeds,
              "groups": groups, "factorial_interactions": interaction, "runs": rows,
              "elapsed_seconds": time.monotonic()-start,
              "scope": "Full physical trajectories retained. Protected contacts and founder turnover do not establish an enclosing membrane, aging, functional organizational memory or an individual. Heat is not total ensemble entropy production. Fuel/function ratios alone do not prove matched-function efficiency."}
    (ROOT/(stage+"-results.json")).write_text(json.dumps(result, indent=2)+"\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("pilot", "ensemble"))
    parser.add_argument("--replicates", type=int, required=True)
    args = parser.parse_args()
    execute(args.stage, args.replicates)
