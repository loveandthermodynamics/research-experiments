# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "scipy", "matplotlib"]
# ///
"""Expose the exchange/contagion tradeoff; retain every declared grid result."""

from dataclasses import asdict, replace
import json
from pathlib import Path

from finite_resource_groups import Parameters, run


def main():
    base = Parameters(damage_good=0.002, damage_bad=0.002,
                      correlated_shock_rate=0.15)
    grid = []
    for contagion in (0.0, 0.3, 0.6, 1.0, 3.0):
        row = []
        for permeability in (0, 0.1, 0.3, 1, 3, 10, 30, 100):
            p = replace(base, contagion_rate=contagion, permeability=permeability)
            val, _ = run(p, (1, 0, -1), horizon=30, samples=241)
            val.update(contagion=contagion, permeability=permeability)
            row.append(val); grid.append(val)
        print(json.dumps({"contagion": contagion,
                          "exergy_integrals": [r["active_exergy_time_integral"] for r in row],
                          "capacity_integrals": [r["active_site_time_integral"] for r in row]}),
              flush=True)
    Path(__file__).with_name("boundary-scan-results.json").write_text(json.dumps({
        "base_parameters": asdict(base), "grid": grid,
        "all_grid_results_retained": True,
        "note": "Permeability is a fixed kinetic parameter; its evolution is not demonstrated.",
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
