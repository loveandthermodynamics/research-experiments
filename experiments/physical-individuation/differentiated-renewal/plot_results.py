"""Scientific figures from frozen saved studies and a labeled approximation."""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent


def load(name):
    return json.loads((ROOT/name).read_text())


def intervals(ax, rows, labels):
    for i, row in enumerate(rows):
        mean = row["mean"]
        lo, hi = row["t_95_interval"]
        ax.errorbar(mean, i, xerr=[[mean-lo], [hi-mean]], fmt="o", capsize=4,
                    color="#236b69" if lo > 0 else "#a75f44")
    ax.axvline(0, color=".6", linewidth=1)
    ax.set(yticks=range(len(rows)), yticklabels=labels, ylim=(-.6, len(rows)-.4))


def main():
    response, causal = load("confirmation-results.json"), load("causal-support-results.json")
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.3))
    a, b, c, d = axes.flat
    intervals(a, response["primary_results"], ["Response effect with renewal", "Response x production"])
    a.set(title="A. Functional response comparison\n128 fresh paired seeds", xlabel="Difference in integrated template count")
    colors = ["#236b69", "#8aa6a1", "#a75f44", "#d3a18b"]
    labels = [f"{g['response_mode'].capitalize()}\nP {'on' if g['produce_L'] else 'off'}" for g in response["groups"]]
    b.bar(range(4), [g["statistics"]["integral_R"]["mean"] for g in response["groups"]], color=colors)
    b.set(xticks=range(4), xticklabels=labels, ylabel="Mean integrated templates", title="B. All four arm means\nSame supplied drive, different realized work")
    x = np.arange(4)
    for offset, key, label, color in [(-.24, "focal_appearances", "Paid mutant by t20", "#8aa6a1"),
            (0, "focal_establishments_over_both_alleles", "Subtree reaches 4", "#236b69"),
            (.24, "focal_survivors", "Subtree survives t100", "#a75f44")]:
        c.bar(x+offset, [g[key]/128 for g in response["groups"]], width=.23, label=label, color=color)
    c.set(xticks=x, xticklabels=labels, ylabel="Fraction of all 128 seeds", ylim=(0, .6),
          title="C. Appearance and ancestry expansion\nBoth alleles count in the focal copying subtree")
    c.legend(frameon=False, fontsize=8, loc="upper right")
    intervals(d, causal["primary_results"], ["Sham minus near lesion", "Remote minus near lesion"])
    d.set(title="D. External fixed-location support probe\n64 preparations, all failures retained", xlabel="Difference in integrated local templates")
    fig.suptitle("Differentiated spatial renewal: declared comparisons and their limits", fontsize=14, y=.985)
    fig.text(.5, .012, "95% t intervals; paired seeds are the statistical units. Ancestry expansion is not a functional-benefit proof.\n"
             "The lesion compares fixed helper locations with matched initial total removed propensity; apparatus costs are outside the model.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .072, 1, .95), h_pad=2, w_pad=2)
    fig.savefig(ROOT/"differentiated-renewal-results.png", dpi=170)

    fig, ax = plt.subplots(figsize=(7.8, 4.5))
    r = np.linspace(.2, 2.4, 1000)
    for a_value, color, label in [(3., "#236b69", "A=3: two positive stationary radii"),
                                 (1.5, "#a75f44", "A=1.5: below fold threshold")]:
        ax.plot(r, a_value/r-1/r**2-r, label=label, color=color)
    roots = [2*np.cos(4*np.pi/9), 2*np.cos(2*np.pi/9)]
    ax.scatter(roots, [0, 0], c=["#a75f44", "#236b69"], zorder=3)
    ax.annotate("Unstable critical size", (roots[0], 0), xytext=(.5, -2.2),
                arrowprops={"arrowstyle": "->", "color": ".5"}, fontsize=9)
    ax.annotate("Radially stable size", (roots[1], 0), xytext=(1.4, 1.2),
                arrowprops={"arrowstyle": "->", "color": ".5"}, fontsize=9)
    ax.axhline(0, color=".6", linewidth=1)
    ax.set(xlabel="Radius R (consistent illustrative units)", ylabel="Radial growth rate",
           title="Conditional radial flux balance: Rdot = A/R - C/R² - BR\nB=C=1; stable size can coexist with material turnover",
           ylim=(-4, 3), xlim=(.2, 2.4))
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    fig.text(.5, .015, "An approximation assuming a phase interface and flux law. It does not certify shape stability, division or our lattice.",
             ha="center", fontsize=8)
    fig.tight_layout(rect=(0, .05, 1, 1))
    fig.savefig(ROOT/"radial-flux-bridge.png", dpi=170)
    print("Saved both figures from fixed evidence and the explicitly labeled radial approximation.")


if __name__ == "__main__":
    main()
