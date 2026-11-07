"""Standalone scientific figures using every fixed seed of the joined study."""
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
COLORS = dict(specialist="#176c8e", generalist="#8b4d9d", X_only="#555555", Y_only="#c66a21")
LABELS = dict(specialist="Specialist", generalist="Full promiscuity", X_only="X activity only", Y_only="Y activity only")


def main():
    data = json.loads((ROOT/"analysis.json").read_text())
    plt.rcParams.update({"font.size":10, "axes.spines.top":False, "axes.spines.right":False})
    metrics = [("template_density_near_dense_B", "Templates near dense B / area"),
               ("template_density", "All templates / area"),
               ("large_B_mass_fraction", "B mass in clusters ≥9"),
               ("spectral_wavelength_low_k", "Low-k wavelength (lattice units)")]
    fig, axes = plt.subplots(4, 2, figsize=(11, 12), sharex=True, constrained_layout=True)
    for col, side in enumerate((16,32)):
        for arm in COLORS:
            rows = [r for r in data["rows"] if r["side"] == side and r["arm"] == arm]
            times = [s["time"] for s in rows[0]["timeseries"]]
            for line,(metric,label) in enumerate(metrics):
                values = np.asarray([[np.nan if s[metric] is None else s[metric]
                                      for s in row["timeseries"]] for row in rows])
                valid = np.all(np.isfinite(values),axis=0)
                means, bands = np.full(len(times),np.nan), np.full(len(times),np.nan)
                means[valid] = values[:,valid].mean(axis=0)
                bands[valid] = t.ppf(.975,7)*values[:,valid].std(axis=0,ddof=1)/np.sqrt(8)
                axes[line,col].plot(times,means,color=COLORS[arm],lw=1.7,label=LABELS[arm])
                axes[line,col].fill_between(times,means-bands,means+bands,color=COLORS[arm],alpha=.12)
                axes[line,col].axvspan(100,200,color="black",alpha=.025)
                axes[line,col].axvspan(300,400,color="black",alpha=.05)
                axes[line,col].set_ylabel(label)
                axes[line,col].grid(alpha=.15)
        axes[0,col].set_title(f"Side{side}: 8 coupled seed rows per arm")
        axes[2,col].set_ylim(-.07,1.07)
        axes[-1,col].set_xlabel("Time")
    axes[0,0].legend(fontsize=8,loc="upper left")
    fig.suptitle("Joined copying and molecular phase: complete bounded comparison\n"
                 "Seed means and pointwise t7 intervals; spatial association is not causal closure",fontsize=12)
    fig.savefig(ROOT/"productive-phase-study.png",dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(4,4,figsize=(12,12),constrained_layout=True)
    material_map = ListedColormap(["#ffffff","#dedede","#176c8e"])
    template_map = ListedColormap(["#ffffff","#cf7f27","#8b4d9d"])
    for size_row,side in enumerate((16,32)):
        for col,arm in enumerate(COLORS):
            with gzip.open(ROOT/"records"/f"L{side}-{arm}-seed730100.json.gz","rt") as stream:
                run = json.load(stream)
            for layer,field,cmap in ((0,"material",material_map),(1,"templates",template_map)):
                row = 2*size_row+layer
                array = np.asarray(run["final_state"][field]).reshape(side,side)
                axes[row,col].imshow(array,cmap=cmap,vmin=0,vmax=2,origin="lower",
                                     extent=(0,side,0,side),interpolation="nearest")
                axes[row,col].set_xticks([0,side//2,side]);axes[row,col].set_yticks([0,side//2,side])
                axes[row,col].set_title(f"L{side} {LABELS[arm]}\n{field}",fontsize=9)
                axes[row,col].set_aspect("equal")
    fig.suptitle("Terminal physical states at time400, first declared seed730100\n"
                 "Material: white S, gray A, blue B. Templates: white D, orange X, purple Y.",fontsize=12)
    fig.savefig(ROOT/"productive-phase-states.png",dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
