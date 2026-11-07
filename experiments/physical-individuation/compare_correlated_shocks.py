# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "scipy", "matplotlib"]
# ///
"""Matched marginal shock exposure; distinct versus shared failure classes.

This is a declared physical vulnerability comparison, not a species-count bonus.
Single-class shock rates per site are identical. The initially supplied response
templates differ only in which sites are affected by the same denaturation
channel. Rare thermal reverse transitions and actual finite-fuel repair remain.
"""

from dataclasses import asdict, replace
import json
from pathlib import Path

from finite_resource_groups import Parameters, run


def main():
    p = Parameters(damage_good=0.002, damage_bad=0.002,
                   correlated_shock_rate=0.15)
    rows = []
    for permeability in (0, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100):
        for label, response in (("monoculture", (1, 1, 1)),
                                ("two_failure_classes", (1, 1, -1)),
                                ("three_failure_classes", (1, 0, -1))):
            val, _ = run(replace(p, permeability=permeability), response,
                         horizon=30, samples=241)
            val.update(label=label, permeability=permeability)
            rows.append(val)
        print(json.dumps({"permeability": permeability,
                          "active_exergy_integrals": [v["active_exergy_time_integral"]
                                                     for v in rows[-3:]],
                          "active_site_integrals": [v["active_site_time_integral"]
                                                    for v in rows[-3:]]}), flush=True)
    Path(__file__).with_name("correlated-shock-results.json").write_text(json.dumps({
        "parameters": asdict(p), "rows": rows,
        "initial_resources_match": True,
        "per_site_marginal_shock_rate_is_identical": True,
        "no_ecological_reward": True,
        "scope": "Finite resource and thermal bath; no hereditary innovation.",
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
