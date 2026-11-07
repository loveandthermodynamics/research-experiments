"""Mechanistically directed fresh test after a retained inconclusive result.

Only the common physical copying coefficient changes, in every response arm.
The test is conditional evidence for a parameter regime, never a universal law.
"""
from dataclasses import asdict, replace
import gzip
from hashlib import sha256
import json
from pathlib import Path

from autonomous_environment import GateParams, simulate_environment
from interface_model import Params, Flags
from interface_readouts import estimate
from run_responses import readout

ROOT = Path(__file__).resolve().parent


def main():
    base = Params(**json.loads((ROOT/"ensemble-protocol.json").read_text())["params"])
    gp = GateParams(**json.loads((ROOT/"response-protocol.json").read_text())["gateparams"])
    base = replace(base, copy_rate=.12)
    seeds = list(range(148001, 148065))
    files = ("interface_model.py", "autonomous_environment.py", "run_rate_test.py", "run_responses.py", "run_interface.py", "interface_readouts.py")
    hashes = {name: sha256((ROOT/name).read_bytes()).hexdigest() for name in files}
    protocol = {"params_before_mutation_override": asdict(base), "gateparams": asdict(gp),
                "seeds": seeds, "horizon": 120., "copying_mutation_levels": [0., .2],
                "source_hashes": hashes,
                "prediction": "Increasing copying in every arm may allow the healthy response variant to replace itself during a hazard, overcoming its measured local subcriticality. Neighborhood and fuel opportunity remain endogenous and can defeat this prediction.",
                "primary_contrasts": ["paired response by variation interaction", "absolute variation effect under differentiated hazards"],
                "scope": "Mechanistically directed new regime after retaining original inconclusive confirmation. No parameter is changed only for the varied arm."}
    (ROOT/"rate-test-protocol.json").write_text(json.dumps(protocol, indent=2)+"\n")
    directory = ROOT/"rate-test"
    directory.mkdir(exist_ok=True)
    rows = []
    for seed in seeds:
        history = None
        for mode in ("common", "differentiated"):
            for varied in (False, True):
                mu = .2 if varied else 0.
                params = replace(base, mutation=((1-mu, mu), (mu, 1-mu)))
                run = simulate_environment(params, seed, 120., Flags(), gp, mode, sample_dt=2.)
                run["source_hashes"] = hashes
                gates = [(e["time"], e["i"], e["after"]) for e in run["ledger"] if e["kind"] == "environment_gate"]
                if history is None:
                    history = gates
                assert gates == history
                path = directory/f"s{seed}-{mode}-v{int(varied)}.json.gz"
                with gzip.open(path, "wt") as f:
                    json.dump(run, f, separators=(",", ":"))
                row = readout(run)
                row.update(mode=mode, varied=varied, trajectory_file=str(path.relative_to(ROOT)))
                rows.append(row)
        print(json.dumps({"seed": seed}), flush=True)
    assert hashes == {name: sha256((ROOT/name).read_bytes()).hexdigest() for name in files}
    groups = []
    for mode in ("common", "differentiated"):
        for varied in (False, True):
            arm = [r for r in rows if (r["mode"], r["varied"]) == (mode, varied)]
            groups.append({"mode": mode, "varied": varied,
                "statistics": {field: estimate([r[field] for r in arm]) for field in
                               ("integral_R", "integral_R1", "net_fuel_consumed", "first_episode_duration")},
                "first_variant_paid_copy_fraction": sum(r["first_variant_channel"] == "copy_R" for r in arm)/len(arm),
                "first_variant_thermal_fraction": sum(r["first_variant_channel"] == "passive" for r in arm)/len(arm)})
    interaction, absolute = [], []
    for seed in seeds:
        r = {(r["mode"], r["varied"]): r["integral_R"] for r in rows if r["seed"] == seed}
        absolute.append(r[("differentiated", True)]-r[("differentiated", False)])
        interaction.append(absolute[-1]-r[("common", True)]+r[("common", False)])
    result = {"protocol": protocol, "groups": groups, "runs": rows,
              "paired_interaction": estimate(interaction), "paired_absolute_variation_effect": estimate(absolute)}
    (ROOT/"rate-test-results.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"interaction": result["paired_interaction"], "absolute": result["paired_absolute_variation_effect"]}))


if __name__ == "__main__":
    main()
