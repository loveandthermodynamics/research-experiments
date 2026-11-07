"""Scientific figure for the two distinct conditional constructions."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


root = Path(__file__).resolve().parent
closure = json.loads((root / "agent-closure-results.json").read_text())
integrity = json.loads((root / "agent-integrity-results.json").read_text())
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.2), layout="constrained")
blue, orange, grey = "#285078", "#c26d2c", "#6a6b6d"
ax = axes[0, 0]
t = closure["cases"]["asymmetric_lesion"]["trajectory"]
ax.plot(t["time"], t["x"], color=blue, label="Type X: initially reduced by 50%")
ax.plot(t["time"], t["y"], color=orange, linestyle="--", label="Type Y")
ax.set(title="A. Reciprocal chemistry restores both types", xlabel="Time", ylabel="Concentration")
ax.set_xlim(0, 30)
ax.legend(fontsize=9, loc="lower right")
ax = axes[0, 1]
for name, label, color in (("above_symmetric_threshold", "2% above threshold", blue),
                           ("below_symmetric_threshold", "2% below threshold", orange)):
    t = closure["cases"][name]["trajectory"]
    ax.plot(t["time"], np.array(t["x"]) + np.array(t["y"]), color=color, label=label)
ax.set(title="B. A driven maintenance threshold", xlabel="Time", ylabel="Total active concentration", yscale="log")
ax.legend(fontsize=9)
ax = axes[1, 0]
c = integrity["contrasts"]["final_lesion_ROI_occupancy"]
rows = [c["sham_minus_wound_without_copying"], c["sham_minus_wound_with_copying"]]
ax.bar([0, 1], [r["mean"] for r in rows], color=[grey, blue], width=.55,
       yerr=[1.96*r["monte_carlo_standard_error"] for r in rows], capsize=5)
ax.set(title="C. Paid copying reduces the local deficit", ylabel="Final sham minus wound occupancy", xticks=[0, 1], xticklabels=["Copying removed", "Copying enabled"])
ax.axhline(0, color=grey, linewidth=.8)
ax.set_xlabel("64 paired seeds; approximate 95% Monte Carlo intervals", fontsize=8.5)
ax = axes[1, 1]
v = integrity["contrasts"]["copying_reduction_of_wound_excess_void_time"]
ax.errorbar(v["mean"], 0, xerr=1.96*v["monte_carlo_standard_error"], fmt="o", color=orange, capsize=7)
ax.axvline(0, color=grey, linestyle="--", linewidth=1)
ax.set(title="D. Geometric restoration is not established", xlabel="Copying improvement in excess void-time", yticks=[], ylim=(-1, 1))
ax.text(.5, .82, "Interval spans zero", ha="center", transform=ax.transAxes)
fig.suptitle("Reciprocal maintenance and local renewal\nConcentration model above; finite-fuel spatial assay below", fontsize=15)
fig.savefig(root / "individuation-tests.png", dpi=160)
print(root / "individuation-tests.png")
