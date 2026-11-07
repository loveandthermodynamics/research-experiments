"""Analysis of the fixed joined productive-phase comparison; no rate feedback."""
from __future__ import annotations

import gzip
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.stats import t

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
from spatial_readouts import observe
from renewal_readouts import measure

METRICS = ("B_fraction", "large_B_mass_fraction", "excess_neighbor_B_fraction",
           "spectral_wavelength_low_k", "largest_B_fraction", "mass_weighted_cluster_size",
           "template_density", "template_density_near_dense_B", "n_X", "n_Y",
           "Y_on_B_excess_fraction", "clusters_ge_9")


def mean(values):
    return None if any(x is None for x in values) else float(np.mean(values))


def ratio(a, b):
    return None if a is None or b is None or b == 0 else a/b


def difference(a, b):
    return None if a is None or b is None else a-b


def estimate(values):
    if any(x is None or not math.isfinite(x) for x in values):
        return dict(estimable=False, n=len(values), values=values,
                    reason="undefined values retained; no seed filtering")
    array = np.asarray(values, dtype=float)
    average = float(array.mean())
    sem = float(array.std(ddof=1)/np.sqrt(len(array)))
    width = float(t.ppf(.975, len(array)-1))*sem
    rng = np.random.default_rng(841701)
    bootstrap = array[rng.integers(len(array), size=(10000, len(array)))].mean(axis=1)
    return dict(estimable=True, n=len(values), values=values, mean=average, sem=sem,
                t95=[average-width, average+width],
                paired_seed_bootstrap95=list(map(float, np.quantile(bootstrap, [.025, .975]))))


def gate(name, values, relation, threshold):
    result = estimate(values)
    result.update(name=name, relation=relation, threshold=threshold)
    result["met"] = bool(result["estimable"] and
                         (result["t95"][0] > threshold if relation == "lower_gt"
                          else result["t95"][1] <= threshold))
    return result


def finite_ablation_bound(params, horizon):
    k = params["side"]**2
    c0 = math.exp(-params["E_X"]/params["theta"])+math.exp(-params["E_Y"]/params["theta"])
    cb = (math.exp(-(params["E_X"]-params["g_X"])/params["theta"])
          +math.exp(-(params["E_Y"]-params["g_Y"])/params["theta"]))
    assert params["background_AB_rate"] == 0 and params["E_X"] > max(0, params["g_X"])
    assert params["E_Y"] > max(0, params["g_Y"])
    stationary = (c0*k+max(0., cb-c0)*params["material_count"])/(1+c0)
    correction = ((params["material_count"]+params["initial_X"]+params["initial_Y"])
                  /(params["loss_rate"]*(1+c0)*horizon))
    return stationary+correction


def integrate_events(run):
    n = sum(x != 0 for x in run["initial_state"]["templates"])
    n_b = sum(x == 2 for x in run["initial_state"]["material"])
    previous, area_n, area_b = 0., 0., 0.
    currents = {1:0, 2:0}
    for event in run["ledger"]:
        now = event["time"]
        area_n += n*(now-previous)
        area_b += n_b*(now-previous)
        if event["instance_created"]:
            n += 1
        if event["instance_removed"]:
            n -= 1
        if event["kind"] in ("copy", "phase_passive"):
            n_b += int(event["after"]["material"][0] == 2)-int(event["before"]["material"][0] == 2)
        if event["kind"] == "copy" and 300. < now <= 400.:
            currents[event["activity"]] += event["direction"]
        previous = now
    area_n += n*(run["horizon"]-previous)
    area_b += n_b*(run["horizon"]-previous)
    assert n == run["final_summary"]["X"]+run["final_summary"]["Y"]
    assert n_b == run["final_summary"]["B"]
    k = run["params"]["side"]**2
    return dict(template_time=area_n, B_time=area_b, whole_mean_template_density=area_n/run["horizon"]/k,
                late_net_X_current_density=currents[1]/(100*k),
                late_net_Y_current_density=currents[2]/(100*k))


def read_one(path):
    with gzip.open(path, "rt") as stream:
        run = json.load(stream)
    side = run["params"]["side"]
    series = []
    for sample in run["samples"]:
        row = observe(np.asarray(sample["material"]).reshape(side, side),
                      np.asarray(sample["templates"]).reshape(side, side))
        row["time"] = sample["time"]
        series.append(row)
    windows = {}
    for name, bounds in (("early", (100, 200)), ("late", (300, 400))):
        selected = [s for s in series if bounds[0] <= s["time"] <= bounds[1]]
        assert len(selected) == 21
        windows[name] = {metric:mean([s[metric] for s in selected]) for metric in METRICS}
    integrals = integrate_events(run)
    renewal = measure(run, appearance_deadline=20., establishment_size=4, include_details=False)
    assert math.isclose(integrals["template_time"], renewal["integral_templates"], rel_tol=2e-11, abs_tol=1e-7)
    assert math.isclose(integrals["B_time"], renewal["integral_B"], rel_tol=2e-11, abs_tol=1e-7)
    return dict(side=side, arm=run["arm"], seed=run["seed"], path=str(path.relative_to(ROOT)),
                record_sha256=sha256(path.read_bytes()).hexdigest(), windows=windows, timeseries=series,
                growth=ratio(windows["late"]["spectral_wavelength_low_k"],
                             windows["early"]["spectral_wavelength_low_k"]), integrals=integrals,
                finite_ablation_bound_density=finite_ablation_bound(run["params"], run["horizon"])/(side**2),
                renewal=renewal, physical_final=run["final_summary"])


def main():
    protocol = json.loads((ROOT/"protocol.json").read_text())
    execution = json.loads((ROOT/"execution.json").read_text())
    rows = [read_one(ROOT/r["path"]) for r in execution["records"]]
    index = {(r["side"], r["arm"], r["seed"]):r for r in rows}
    seeds = protocol["seeds"]
    get = lambda side, arm, seed:index[(side, arm, seed)]
    gates_a, gates_b = [], []
    for side in (16, 32):
        for metric, threshold in (("B_fraction", .05), ("large_B_mass_fraction", .5),
                                  ("excess_neighbor_B_fraction", .15)):
            gates_a.append(gate(f"L{side}_late_{metric}",
                           [get(side,"specialist",s)["windows"]["late"][metric] for s in seeds],"lower_gt",threshold))
        gates_a.append(gate(f"L{side}_wavelength_growth",[get(side,"specialist",s)["growth"] for s in seeds],"upper_le",1.25))
        for activity in ("X", "Y"):
            metric = f"late_net_{activity}_current_density"
            gates_a.append(gate(f"L{side}_{metric}",
                           [get(side,"specialist",s)["integrals"][metric] for s in seeds],"lower_gt",0.))
        for control in ("X_only", "Y_only"):
            gates_a.append(gate(f"L{side}_late_template_density_minus_{control}",
                           [difference(get(side,"specialist",s)["windows"]["late"]["template_density"],
                                       get(side,control,s)["windows"]["late"]["template_density"]) for s in seeds],"lower_gt",0.))
        gates_a.append(gate(f"L{side}_whole_template_density_minus_finite_ablation_bound",
                       [get(side,"specialist",s)["integrals"]["whole_mean_template_density"]
                        -get(side,"specialist",s)["finite_ablation_bound_density"] for s in seeds],"lower_gt",0.))
        gates_b.append(gate(f"L{side}_specialist_minus_generalist_late_dense_associated_template_density",
                       [difference(get(side,"specialist",s)["windows"]["late"]["template_density_near_dense_B"],
                                   get(side,"generalist",s)["windows"]["late"]["template_density_near_dense_B"])
                        for s in seeds],"lower_gt",0.))
    gates_a.append(gate("specialist_late_wavelength_L32_over_L16",
                   [ratio(get(32,"specialist",s)["windows"]["late"]["spectral_wavelength_low_k"],
                          get(16,"specialist",s)["windows"]["late"]["spectral_wavelength_low_k"]) for s in seeds],"upper_le",1.5))
    gates_a.append(gate("specialist_late_largest_B_fraction_L16_minus_L32",
                   [difference(get(16,"specialist",s)["windows"]["late"]["largest_B_fraction"],
                               get(32,"specialist",s)["windows"]["late"]["largest_B_fraction"]) for s in seeds],"lower_gt",0.))
    assert len(gates_a) == 20 and len(gates_b) == 2
    summaries = {}
    for side in (16, 32):
        for arm in protocol["arms"]:
            selected = [get(side,arm,s) for s in seeds]
            summaries[f"L{side}-{arm}"] = dict(
                late={m:estimate([r["windows"]["late"][m] for r in selected]) for m in METRICS},
                growth=estimate([r["growth"] for r in selected]),
                whole_mean_template_density=estimate([r["integrals"]["whole_mean_template_density"] for r in selected]),
                chemical_work_per_area=estimate([r["physical_final"]["chemical_work"]/side**2 for r in selected]),
                heat_to_bath_per_area=estimate([r["physical_final"]["heat_to_bath"]/side**2 for r in selected]),
                net_X_current=estimate([r["integrals"]["late_net_X_current_density"] for r in selected]),
                net_Y_current=estimate([r["integrals"]["late_net_Y_current_density"] for r in selected]),
                finite_ablation_bound_density=selected[0]["finite_ablation_bound_density"],
                renewal=[r["renewal"] for r in selected])
    output = dict(protocol_sha256=sha256((ROOT/"protocol.json").read_bytes()).hexdigest(),
                  analysis_sha256=sha256(Path(__file__).read_bytes()).hexdigest(), unique_records=len(rows),
                  independent_seed_rows=len(seeds), claim_A_met=all(g["met"] for g in gates_a),
                  claim_B_met=all(g["met"] for g in gates_b), claim_A_gates=gates_a, claim_B_gates=gates_b,
                  summary=summaries, rows=rows,
                  scope="Finite-horizon domain-level phase/renewal and near-dense-B association only. No local component closure, asymptotic phase, species necessity, total availability, aging, exploiter resistance or moral agency proven.")
    (ROOT/"analysis.json").write_text(json.dumps(output,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print(json.dumps({k:output[k] for k in ("claim_A_met","claim_B_met","unique_records")}),flush=True)
    for g in gates_a+gates_b:
        print(json.dumps({k:g.get(k) for k in ("name","estimable","mean","t95","met")}),flush=True)


if __name__ == "__main__":
    main()
