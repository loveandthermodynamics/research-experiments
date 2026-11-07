"""Standalone scientific figures from all fixed study seeds."""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from scipy.stats import t

ROOT = Path(__file__).resolve().parent
COLORS = {"full": "#1c7293", "no_drive": "#575757", "no_affinity": "#c16e22"}
LABELS = {"full": "Driven, binding enzyme", "no_drive": "No chemical drive",
          "no_affinity": "Driven, no enzyme binding"}


def main():
    data = json.loads((ROOT / "analysis.json").read_text())
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(3, 2, figsize=(11, 10), sharex=True, constrained_layout=True)
    metrics = [("large_B_mass_fraction", "B mass in clusters ≥9"),
               ("excess_neighbor_B_fraction", "B neighbor excess"),
               ("spectral_wavelength_low_k", "Low-k wavelength (lattice units)")]
    for col, side in enumerate((16, 32)):
        for arm in COLORS:
            rows = [r for r in data["rows"] if r["side"] == side and r["arm"] == arm]
            times = [r["time"] for r in rows[0]["timeseries"]]
            for line, (metric, ylabel) in enumerate(metrics):
                array = np.asarray([[np.nan if s[metric] is None else s[metric]
                                     for s in row["timeseries"]] for row in rows])
                # Undefined times stay blank; do not discard whole seeds.
                defined = np.all(np.isfinite(array), axis=0)
                mean, band = np.full(len(times), np.nan), np.full(len(times), np.nan)
                mean[defined] = array[:, defined].mean(axis=0)
                band[defined] = t.ppf(.975, 7) * array[:, defined].std(axis=0, ddof=1) / np.sqrt(8)
                axes[line, col].plot(times, mean, color=COLORS[arm], label=LABELS[arm], lw=1.7)
                axes[line, col].fill_between(times, mean-band, mean+band, color=COLORS[arm], alpha=.13)
                axes[line, col].axvspan(100, 200, color="black", alpha=.025)
                axes[line, col].axvspan(300, 400, color="black", alpha=.05)
                axes[line, col].set_ylabel(ylabel)
                axes[line, col].grid(alpha=.15)
        axes[0, col].set_title(f"Side {side}: 8 fixed seeds per arm")
        axes[-1, col].set_xlabel("Time")
        axes[0, col].set_ylim(-.07, 1.07)
        axes[1, col].set_ylim(-.1, .8)
    axes[0, 0].legend(fontsize=8, loc="upper left")
    fig.suptitle("Driven lattice organization: complete bounded comparison\n"
                 "Means and pointwise t7 intervals across seed rows; no finite-size claim assumed", fontsize=13)
    fig.savefig(ROOT / "phase-study.png", dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(10.5, 7), constrained_layout=True)
    color_map = ListedColormap(["#ffffff", "#d4d4d4", "#1c7293"])
    for row, side in enumerate((16, 32)):
        for col, arm in enumerate(COLORS):
            path = ROOT / "records" / f"L{side}-{arm}-seed720100.json.gz"
            with gzip.open(path, "rt") as stream:
                run = json.load(stream)
            material = np.asarray(run["final_state"]["material"]).reshape(side, side)
            enzyme = np.asarray(run["final_state"]["enzyme"]).reshape(side, side)
            axes[row, col].imshow(material, cmap=color_map, vmin=0, vmax=2,
                                  interpolation="nearest", origin="lower", extent=(0, side, 0, side))
            rr, cc = np.nonzero(enzyme)
            axes[row, col].scatter(cc+.5, rr+.5, s=20 if side == 16 else 8,
                                  facecolors="none", edgecolors="#da6c19", linewidths=.9)
            axes[row, col].set_title(f"L{side}: {LABELS[arm]}", fontsize=9)
            axes[row, col].set_xticks([0, side//2, side]); axes[row, col].set_yticks([0, side//2, side])
            axes[row, col].set_aspect("equal")
    fig.suptitle("Terminal microscopic states at time 400\n"
                 "First declared seed 720100 in every arm; blue B, gray A, orange rings E", fontsize=12)
    fig.savefig(ROOT / "phase-states.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
