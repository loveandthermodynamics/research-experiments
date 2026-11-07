# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "scipy", "matplotlib"]
# ///
"""Refine quadrature, retain adverse cases, and render matched physical results."""

from dataclasses import asdict, replace
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from finite_resource_groups import Parameters, run


def main():
    p = Parameters(damage_good=0.002, damage_bad=0.002,
                   correlated_shock_rate=0.15, permeability=3)
    labels = (("One shared vulnerability", (1, 1, 1)),
              ("Two vulnerabilities", (1, 1, -1)),
              ("Three vulnerabilities", (1, 0, -1)))
    rows, tracks, convergence = [], {}, []
    for label, responses in labels:
        coarse, _ = run(p, responses, samples=481)
        fine, tr = run(p, responses, samples=1921)
        diff = abs(fine["active_exergy_time_integral"]
                   - coarse["active_exergy_time_integral"])
        assert diff / fine["active_exergy_time_integral"] < 0.0001
        fine["label"] = label
        rows.append(fine); tracks[label] = tr
        convergence.append({"label": label,
                            "quad_481_vs_1921_absolute_difference": diff})
    # Coupling does not increase response-diversity value by definition: another
    # environment implemented in the same reversible model gives an adverse case.
    adverse = []
    for label, responses in labels[:2]:
        val, _ = run(Parameters(), responses, samples=481)
        val["label"] = label
        adverse.append(val)
    assert adverse[0]["active_exergy_time_integral"] > adverse[1]["active_exergy_time_integral"]
    assert rows[2]["active_exergy_time_integral"] > rows[0]["active_exergy_time_integral"]
    assert np.ptp([r["initial_total_available_free_energy"] for r in rows]) < 1e-12
    assert np.ptp([r["initial_active_available_free_energy"] for r in rows]) < 1e-12
    out = {"parameters": asdict(p), "comparison": rows,
           "quadrature_refinement": convergence,
           "active_availability_gain_three_vs_one": rows[2]["active_exergy_time_integral"]
                                                       / rows[0]["active_exergy_time_integral"] - 1,
           "active_function_gain_three_vs_one": rows[2]["active_site_time_integral"]
                                                   / rows[0]["active_site_time_integral"] - 1,
           "adverse_telegraph_gate_case": adverse,
           "checks_passed": True,
           "scope": "Actual common-reference nonequilibrium available free energy within the coarse model; not a universal maximization principle."}
    folder = Path(__file__).parent
    (folder / "physical-results-verified.json").write_text(json.dumps(out, indent=2) + "\n")
    np.savez_compressed(folder / "physical-trajectories.npz",
                        **{f"{i}_{key}": value for i, tr in enumerate(tracks.values())
                           for key, value in tr.items()})
    fig, ax = plt.subplots(1, 3, figsize=(14.7, 4.3), constrained_layout=True)
    for label, tr in tracks.items():
        ax[0].plot(tr["time"], tr["active_exergy"], label=label)
        ax[1].plot(tr["time"], tr["fuel_remaining"])
    ax[0].set(xlabel="Time in model units", ylabel="Active available free energy",
              title="Survivors support paid regeneration")
    ax[0].legend(fontsize=8)
    ax[1].set(xlabel="Time in model units", ylabel="Fuel packets remaining",
              title="Identical finite fuel at the start")
    ax[2].bar(["Shared", "Two classes", "Three classes"],
              [r["active_exergy_time_integral"] for r in rows],
              color=["#999999", "#4386b0", "#218e78"])
    ax[2].set(ylabel="Integral of active available free energy",
              xlabel="Failure dependencies among three sites",
              title="A physical benefit of different vulnerabilities")
    for i, r in enumerate(rows):
        ax[2].text(i, r["active_exergy_time_integral"] + 1.5,
                   f'{r["active_exergy_time_integral"]:.2f}', ha="center", fontsize=9)
    ax[2].set_ylim(0, max(r["active_exergy_time_integral"] for r in rows) * 1.14)
    fig.savefig(folder / "paid-regeneration-and-exergy.png", dpi=175)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
