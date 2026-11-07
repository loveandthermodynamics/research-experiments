"""Plot the pilot without implying an ensemble advantage or reproduction."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np


root = Path(__file__).resolve().parent
data = json.loads((root / "agent-lattice-results.json").read_text())
fig = plt.figure(figsize=(11.5, 10.5), layout="constrained")
grid = fig.add_gridspec(3, 3, height_ratios=[1.0, 1.1, 1.1])
ax = fig.add_subplot(grid[0, :])
colors = {"compact": "#285078", "dispersed": "#a3541d"}
for layout in colors:
    for suffix, style, text in (("no_mutation", "--", "error-free copying"),
                                ("paid_copying_errors", "-", "copying errors 0.05")):
        record = data["runs"][layout + "_" + suffix][0]
        ax.step([row["time"] for row in record["trajectory"]],
                [row["active_templates"] for row in record["trajectory"]],
                where="post", color=colors[layout], linestyle=style,
                linewidth=1.7, label=layout.capitalize() + ": " + text)
ax.set(xlabel="Time in model units", ylabel="Active templates", xlim=(0, 120))
ax.set_title("One trajectory per condition: variation helps or hurts depending on the run", loc="left", fontsize=12)
ax.spines[["top", "right"]].set_visible(False)
ax.grid(alpha=.18)
ax.legend(fontsize=9, ncol=2, loc="upper left")
cmap = ListedColormap(["#f4f1eb", "#285078", "#d67b34"])
norm = BoundaryNorm([-1.5, -.5, .5, 1.5], 3)
for row, layout in enumerate(("compact", "dispersed"), start=1):
    record = data["runs"][layout + "_paid_copying_errors"][0]
    for col, snapshot in enumerate(record["snapshots"]):
        frame = fig.add_subplot(grid[row, col])
        frame.imshow(np.array(snapshot["material"]), cmap=cmap, norm=norm, interpolation="nearest")
        frame.set(xticks=[], yticks=[])
        active = int((np.array(snapshot["material"]) >= 0).sum())
        frame.set_title(f"{layout.capitalize()}, t = {snapshot['time']:g}: {active} templates", fontsize=10)
fig.suptitle("Local copying and binding with finite fuel\nPale = unassembled material; blue = ancestral response; orange = variant", fontsize=14)
fig.savefig(root / "local-lattice-trajectories.png", dpi=150)
print(root / "local-lattice-trajectories.png")
