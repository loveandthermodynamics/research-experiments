"""Secondary observations and scientific figure from retained ensemble paths."""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from interface_readouts import estimate

ROOT = Path(__file__).resolve().parent


def turnover(run):
    founders = {x for x in run["initial_state"]["instance_tags"] if x}
    founder_live = set(founders)
    origin = {x: True for x in founders}
    descendants_live = set()
    material = np.asarray(run["initial_state"]["material"]).reshape(run["params"]["side"], -1).copy()
    time = 0.0
    interval_after_founder_loss = inherited_R_time = 0.0
    first_founder_free_continuation = None
    protected_continuation = contact_integral_after_founders = 0.0
    l_births = set()
    paid_l_lost = 0
    for event in run["ledger"]:
        duration = event["time"]-time
        if not founder_live and descendants_live:
            interval_after_founder_loss += duration
            inherited_R_time += duration*len(descendants_live)
            if first_founder_free_continuation is None:
                first_founder_free_continuation = time
            rmask = (material == 1) | (material == 2)
            lmask = material == 3
            contacts = sum(int(np.sum(rmask*np.roll(lmask, -1, axis)+lmask*np.roll(rmask, -1, axis)))
                           for axis in (0, 1))
            if contacts:
                protected_continuation += duration
                contact_integral_after_founders += duration*contacts
        time = event["time"]
        flat = material.ravel()
        if event["kind"] == "material_swap":
            i, j = event["i"], event["j"]
            flat[i], flat[j] = flat[j], flat[i]
        elif event["kind"] != "fuel_swap":
            flat[event["i"]] = event["after"]
        removed, created = event["instance_removed"], event["instance_created"]
        if removed:
            founder_live.discard(removed)
            descendants_live.discard(removed)
            if removed in l_births:
                paid_l_lost += 1
        if created:
            inherited = event["kind"] == "copy_R" and origin.get(event["parent_instance"], False)
            origin[created] = inherited
            if inherited:
                descendants_live.add(created)
            if event["kind"] == "produce_L":
                l_births.add(created)
    if not founder_live and descendants_live:
        interval_after_founder_loss += run["horizon"]-time
        inherited_R_time += (run["horizon"]-time)*len(descendants_live)
        if first_founder_free_continuation is None:
            first_founder_free_continuation = time
        rmask = (material == 1) | (material == 2)
        lmask = material == 3
        contacts = sum(int(np.sum(rmask*np.roll(lmask, -1, axis)+lmask*np.roll(rmask, -1, axis)))
                       for axis in (0, 1))
        if contacts:
            protected_continuation += run["horizon"]-time
            contact_integral_after_founders += (run["horizon"]-time)*contacts
    return {"seed": run["seed"], "flags": run["flags"],
            "continuation_after_all_founder_instances_lost": interval_after_founder_loss,
            "founder_lineage_R_time_after_founder_instance_loss": inherited_R_time,
            "RL_contact_continuation_after_founders_lost": protected_continuation,
            "integral_RL_contacts_after_founder_loss": contact_integral_after_founders,
            "first_founder_free_continuation": first_founder_free_continuation,
            "paid_L_instances_lost": paid_l_lost,
            "paid_L_instances_created": len(l_births),
            "scope": "Copy ancestry and molecular-instance replacement, not collective heredity or functional regulatory information."}


def main():
    result = json.loads((ROOT/"ensemble-results.json").read_text())
    paths, rows = [], []
    grouped = {}
    for row in result["runs"]:
        with gzip.open(ROOT/row["trajectory_file"], "rt") as f:
            run = json.load(f)
        flag = (run["flags"]["produce_L"], run["flags"]["protect_R"])
        grouped.setdefault(flag, []).append(run)
        rows.append(turnover(run))
        paths.append(row)
    statistics = []
    for flag, runs in grouped.items():
        chosen = [r for r in rows if (r["flags"]["produce_L"], r["flags"]["protect_R"]) == flag]
        integrals = [r["integral_R"] for r in paths if (r["flags"]["produce_L"], r["flags"]["protect_R"]) == flag]
        fuels = [r["net_fuel_consumed"] for r in paths if (r["flags"]["produce_L"], r["flags"]["protect_R"]) == flag]
        statistics.append({"flags": {"produce_L": flag[0], "protect_R": flag[1]},
            "founder_free_continuation_frequency": sum(r["continuation_after_all_founder_instances_lost"] > 0 for r in chosen)/len(chosen),
            "continuation_time": estimate([r["continuation_after_all_founder_instances_lost"] for r in chosen]),
            "inherited_R_time": estimate([r["founder_lineage_R_time_after_founder_instance_loss"] for r in chosen]),
            "RL_continuation_time_after_founders": estimate([r["RL_contact_continuation_after_founders_lost"] for r in chosen]),
            "RL_continuation_frequency_after_founders": sum(r["RL_contact_continuation_after_founders_lost"] > 0 for r in chosen)/len(chosen),
            "paid_L_instances_lost": estimate([r["paid_L_instances_lost"] for r in chosen]),
            "aggregate_fuel_standard_energy_per_R_time": 12*sum(fuels)/sum(integrals)})
    by_seed = {(r["seed"], r["flags"]["produce_L"], r["flags"]["protect_R"]): r for r in paths}
    paired = []
    for seed in result["seeds"]:
        paired.append(by_seed[(seed, True, True)]["integral_R"]-by_seed[(seed, True, False)]["integral_R"]
                      -by_seed[(seed, False, True)]["integral_R"]+by_seed[(seed, False, False)]["integral_R"])
    values = np.array([[by_seed[(seed, True, True)]["integral_R"],
                        by_seed[(seed, False, False)]["integral_R"],
                        by_seed[(seed, True, True)]["net_fuel_consumed"],
                        by_seed[(seed, False, False)]["net_fuel_consumed"]]
                       for seed in result["seeds"]])
    rng = np.random.default_rng(335)
    boot = values[rng.integers(0, len(values), (10000, len(values)))].mean(axis=1)
    reduction = 1-(boot[:, 2]/boot[:, 0])/(boot[:, 3]/boot[:, 1])
    economy = {"measure": "ratio of ensemble total standard fuel energy to ensemble total active-template time",
               "reduction_fraction": float(1-(values[:, 2].sum()/values[:, 0].sum())/(values[:, 3].sum()/values[:, 1].sum())),
               "paired_bootstrap_95_interval": np.quantile(reduction, [.025, .975]).tolist(),
               "replicates": 10000, "bootstrap_seed": 335,
               "scope": "An activity-normalized standard-energy ratio, not matched-reliability thermodynamic efficiency."}
    summary = {"turnover_statistics": statistics, "paired_factorial_interaction_R_time": estimate(paired),
               "activity_normalized_fuel_economy": economy,
               "observations": rows,
               "scope": "Secondary readouts from the frozen 64-seed ensemble. Ratio of aggregate fuel standard energy to activity is not a matched-function total-entropy-production measure."}
    (ROOT/"turnover-analysis.json").write_text(json.dumps(summary, indent=2)+"\n")
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.1))
    colors = {(True, True): "#176b5b", (True, False): "#b6513b", (False, True): "#405e91"}
    labels = {(True, True): "Produced protection", (True, False): "Produced material, barriers off",
              (False, True): "Production off"}
    for flag, color in colors.items():
        runs = grouped[flag]
        times = np.array([s["time"] for s in runs[0]["samples"]])
        r = np.array([[s["R0"]+s["R1"] for s in run["samples"]] for run in runs])
        ll = np.array([[s["L"] for s in run["samples"]] for run in runs])
        for axis, y in ((axes[0], r), (axes[1], ll)):
            mean, se = y.mean(axis=0), y.std(axis=0, ddof=1)/np.sqrt(len(y))
            axis.plot(times, mean, color=color, label=labels[flag], linewidth=2)
            axis.fill_between(times, mean-1.96*se, mean+1.96*se, color=color, alpha=.15)
        arm = [x for x in paths if (x["flags"]["produce_L"], x["flags"]["protect_R"]) == flag]
        axes[2].scatter([x["net_fuel_consumed"] for x in arm], [x["integral_R"] for x in arm],
                        s=13, color=color, alpha=.42)
        axes[2].scatter(np.mean([x["net_fuel_consumed"] for x in arm]),
                        np.mean([x["integral_R"] for x in arm]), color=color, marker="X", s=100,
                        edgecolor="white", linewidth=.6)
    axes[0].set(title="Templates renewed through finite-fuel episodes", xlabel="Time", ylabel="Active templates")
    axes[1].set(title="Protective material produced from zero", xlabel="Time", ylabel="Protective units L")
    axes[2].set(title="Activity and fuel use are distinct outcomes", xlabel="Net fuel packets consumed", ylabel="Integrated template activity")
    axes[0].legend(fontsize=8, frameon=False)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(alpha=.15)
    fig.text(.01, .01, "64 new seeds per arm; shading: approximate 95% MC intervals. Both-off arm matches production-off here. No organism, membrane or senescence claim.", fontsize=8)
    fig.tight_layout(rect=(0, .045, 1, 1))
    fig.savefig(ROOT/"interface-maintenance.png", dpi=180)
    plt.close(fig)
    responses = json.loads((ROOT/"response-ensemble-results.json").read_text())
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4))
    for i, (mode, color) in enumerate((("common", "#b6513b"), ("differentiated", "#176b5b"))):
        selected = [g for g in responses["groups"] if g["mode"] == mode]
        for group in selected:
            x = int(group["varied"])+(i-.5)*.12
            stat = group["statistics"]["integral_R"]
            axes[0].errorbar(x, stat["mean"], yerr=1.96*stat["monte_carlo_standard_error"],
                             color=color, fmt="o", capsize=4,
                             label=mode.capitalize()+" hazards" if not group["varied"] else None)
    axes[0].set_xticks([0, 1], ["Zero copying mutation", "Copying mutation 0.2"])
    axes[0].set(title="Functional variants do not yet establish a gain", ylabel="Integrated template activity")
    axes[0].legend(frameon=False, fontsize=8)
    stat = responses["interactions"][0]["paired_interaction_response_by_variation"]
    axes[1].axhline(0, color="#666666", linewidth=1)
    axes[1].errorbar(0, stat["mean"], yerr=1.96*stat["monte_carlo_standard_error"],
                     fmt="o", color="#405e91", capsize=6)
    axes[1].set(xlim=(-.7, .7), xticks=[0], xticklabels=["Response by variation interaction"],
                ylabel="Difference in integrated activity", title="Fresh paired interval includes zero")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=.15)
    fig.text(.015, .01, "64 new paired seeds; approximate 95% Monte Carlo intervals. Thermal variant formation remains possible with zero copying mutation.", fontsize=8)
    fig.tight_layout(rect=(0, .05, 1, 1))
    fig.savefig(ROOT/"response-variation.png", dpi=180)
    plt.close(fig)
    print(json.dumps({"statistics": statistics, "paired_interaction": estimate(paired)}, indent=2))


if __name__ == "__main__":
    main()
