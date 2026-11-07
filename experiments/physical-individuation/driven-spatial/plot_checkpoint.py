"""Render saved results without fitting or running new physical trajectories."""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent


def load(name):
    return json.loads((ROOT / name).read_text())


def interval(ax, mean, limits, y, color, marker="o"):
    ax.errorbar(mean, y, xerr=[[mean-limits[0]], [limits[1]-mean]],
                color=color, fmt=marker, capsize=4, markersize=6)


def main():
    main_result = load("confirmation-results.json")
    replication = load("locality-replication-results.json")
    founders = load("founder-results.json")["runs"]
    cases = load("diversity-invasion-algebra.json")["cases"]
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.5))
    a, b, c, d = axes.flat
    for y, row in enumerate(main_result["primary_confirmation_results"][:2]):
        interval(a, row["mean"], row["approximate_95_percent_normal_interval"], y, "#156b64")
    a.set(yticks=[0, 1], yticklabels=["Production when protected", "Production x protection"],
          xlabel="Difference in integrated template count", ylim=(-.6, 1.6),
          title="A. Driven lattice: renewing protection\n64 fresh preparations")
    a.axvline(0, color=".6", linewidth=1)
    previous = replication["previous_confirmation"]
    row = previous["normal_result_as_saved"]
    interval(b, row["mean"], row["approximate_95_percent_normal_interval"], 0, "#bd8b29")
    interval(b, row["mean"], previous["t_95_interval_sensitivity"], 1, "#a05b36")
    row = replication["replication_estimate"]
    interval(b, row["mean"], row["t_95_interval"], 2, "#354e85")
    b.set(yticks=[0, 1, 2], yticklabels=["Original normal interval (n=64)",
          "Original t sensitivity (n=64)", "Fresh replication, t interval (n=256)"],
          xlabel="Local minus mixed transport: template-time", ylim=(-.65, 2.65),
          title="B. Driven lattice: local transport benefit\nremains unestablished")
    b.axvline(0, color=".6", linewidth=1)
    for produce, protect, label, color in [(True, True, "Full renewal", "#156b64"),
            (False, True, "Production removed", "#bc6546"),
            (True, False, "Protection removed", "#8c79aa"),
            (False, False, "Both removed", "#777777")]:
        rows = sorted((r for r in founders if r["produce"] == produce and r["protect"] == protect),
                      key=lambda r: r["founder_R"])
        c.plot([r["founder_R"] for r in rows], [r["final_R"] for r in rows],
               marker="o", markersize=4, linewidth=1.2, label=label, color=color)
    c.set(xscale="log", yscale="log", xlabel="Supplied founder concentration, initial L = 0",
          ylabel="Template concentration at time 500",
          title="C. Separate well-mixed model\nFinite founder grid, 28 deterministic paths")
    c.legend(frameon=False, fontsize=8)
    values = [r["opposite_response_invasion_coefficient"] for r in cases]
    d.bar([0, 1], values, color=["#156b64", "#bc6546"], width=.48)
    d.axhline(0, color=".6", linewidth=1)
    d.set(xticks=[0, 1], xticklabels=["M = .95\nPositive invasion", "M = 1.5\nNegative invasion"],
          ylabel="Leading rare-type invasion coefficient K",
          title="D. Separate concentration theory\nDiversity advantage changes sign")
    for i, value in enumerate(values):
        d.text(i, value+(.012 if value > 0 else -.012), f"{value:+.3f}", ha="center",
               va="bottom" if value > 0 else "top")
    d.set_ylim(-.19, .25)
    fig.suptitle("Maintained organization and diversity: progress without a universal claim", fontsize=15, y=.985)
    fig.text(.5, .012, "A-B: paired spatial comparisons; 95% intervals. C: deterministic assay, not formation probabilities. "
             "D: asymptotic coefficient; no numerical forcing bound.\nThese are distinct models. "
             "No panel establishes a new individual, an energetic optimum, open-ended biodiversity or moral legitimacy.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .067, 1, .955), h_pad=2.6, w_pad=2.2)
    fig.savefig(ROOT / "driven-checkpoint.png", dpi=170)
    print(ROOT / "driven-checkpoint.png")


if __name__ == "__main__":
    main()
