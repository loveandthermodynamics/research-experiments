"""Show the failed spatial assay beside the separately proved open mechanism."""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from verify_threshold import vector, limiting_material

ROOT = Path(__file__).resolve().parent
plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.facecolor": "white"})
fig, axes = plt.subplots(1, 3, figsize=(16, 4.9), gridspec_kw={"width_ratios": [1.05, 1., 1.15]})
data = json.loads((ROOT/"confirmation-results.json").read_text())
selected = [next(r for r in data["paired_contrasts"] if r["field"] == "integral_R" and r["contrast"] == name)
            for name in ("production_effect_when_protected", "production_protection_interaction")]
for y, row, color in zip([1, 0], selected, ["#726b65", "#236f7c"]):
    lo, hi = row["approximate_95_percent_normal_interval"]
    axes[0].errorbar(row["mean"], y, xerr=[[row["mean"]-lo], [hi-row["mean"]]],
                     fmt="o", color=color, capsize=5, markersize=7)
axes[0].axvline(0, color="#aaa29a", ls="--", lw=1)
axes[0].set(yticks=[1, 0], yticklabels=["Production benefit\nwith protection", "Production ×\nprotection interaction"],
            ylim=(-.65, 1.65), xlim=(-23, 65), xlabel="Change in integrated template activity",
            title="Spatial test: joint claim fails")
axes[0].text(.03, .04, "64 fresh preparations; paired 95% intervals\nBoth intervals had to lie above zero.",
             transform=axes[0].transAxes, fontsize=9, va="bottom")

check = json.loads((ROOT/"threshold-verification.json").read_text())
ell = np.linspace(0, .23, 401)
axes[1].plot(ell, [limiting_material(x) for x in ell], color="#236f7c", lw=2)
axes[1].axhline(.95, color="#b76b34", ls="--", label="Available material M = 0.95")
for point, fill in zip(check["limiting_positive_equilibria"], ["white", "#236f7c"]):
    axes[1].plot(point["L"], .95, "o", mec="#236f7c", mfc=fill, ms=8, zorder=3)
axes[1].plot(check["limiting_fold"]["L"], check["limiting_fold"]["M_min"], "x", color="#726b65", ms=8)
axes[1].set(xlabel="Protective material L", ylabel="Required material F(L)",
            ylim=(.89, 1.13), xlim=(0, .23), title="Open reactor: a persistence threshold")
axes[1].text(.05, .94, "Saddle", fontsize=9)
axes[1].text(.14, .96, "Stable high state", fontsize=9)
axes[1].legend(loc="upper left", fontsize=8, frameon=False)

t = np.linspace(0, 500, 1001)
high = [check["reversible_equilibria"][2][name] for name in ("R", "L")]
for produce, protect, label, color in [(True, True, "Production and protection", "#236f7c"),
                                      (False, True, "Production stopped", "#b76b34"),
                                      (True, False, "Protection removed", "#7766a5")]:
    sol = solve_ivp(lambda time, state: vector(state, produce, protect), (0, 500), high,
                    t_eval=t, rtol=1e-10, atol=1e-12)
    assert sol.success
    axes[2].semilogy(t, sol.y[0], color=color, label=label, lw=2)
axes[2].set(xlabel="Time after intervention", ylabel="Active template R (log scale)",
            title="Reversible open model: high state needs renewal", ylim=(1e-8, .7))
axes[2].legend(loc="upper right", bbox_to_anchor=(1, .80), fontsize=8,
               frameon=True, facecolor="white", framealpha=.95, edgecolor="white")
fig.suptitle("Finite-budget evidence and a separate maintained-flow construction", fontsize=15, y=.99)
fig.text(.5, .015, "The open reactor assumes copying, protection and counted fuel flow. Its threshold proves conditional maintenance, not spatial individuation or biodiversity.",
         ha="center", fontsize=10, color="#504a45")
fig.tight_layout(rect=(.01, .065, .99, .95), w_pad=2)
fig.savefig(ROOT/"maintenance-threshold.png", dpi=150)
plt.close(fig)
