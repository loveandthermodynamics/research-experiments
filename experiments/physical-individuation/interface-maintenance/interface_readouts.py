"""Observations only: these functions never enter the transition rates."""
from __future__ import annotations

import math
import numpy as np


def spatial_readouts(material, side: int) -> dict:
    """D=0, R0=1, R1=2, L=3; periodic square lattice, undirected NN edges."""
    material = np.asarray(material).reshape(side, side)
    r = (material == 1) | (material == 2)
    l = material == 3
    occupied = r | l
    rn = np.zeros_like(material, dtype=int)
    ln = np.zeros_like(material, dtype=int)
    rl = ll = rr = 0
    for axis in (0, 1):
        next_r, next_l = np.roll(r, -1, axis), np.roll(l, -1, axis)
        rl += int(np.sum(r*next_l+l*next_r))
        ll += int(np.sum(l*next_l))
        rr += int(np.sum(r*next_r))
        for shift in (-1, 1):
            rn += np.roll(r, shift, axis)
            ln += np.roll(l, shift, axis)
    side_sites = side*side
    neighbors = [[] for _ in range(side_sites)]
    for i in range(side_sites):
        a, b = divmod(i, side)
        neighbors[i] = [((a+1)%side)*side+b, ((a-1)%side)*side+b,
                        a*side+(b+1)%side, a*side+(b-1)%side]
    unseen = set(map(int, np.flatnonzero(occupied.ravel())))
    sizes, productive = [], []
    flat = material.ravel()
    while unseen:
        origin = unseen.pop()
        stack, sites = [origin], [origin]
        while stack:
            for j in neighbors[stack.pop()]:
                if j in unseen:
                    unseen.remove(j)
                    stack.append(j)
                    sites.append(j)
        sizes.append(len(sites))
        kinds = set(map(int, flat[sites]))
        productive.append(3 in kinds and bool(kinds & {1, 2}))
    expected_rl_random = 4*float(r.sum())*float(l.sum())/(side_sites-1)
    return {
        "R0": int(np.sum(material == 1)), "R1": int(np.sum(material == 2)),
        "R": int(r.sum()), "L": int(l.sum()), "D": int(np.sum(material == 0)),
        "RL_edges": rl, "LL_edges": ll, "RR_edges": rr,
        "R_with_L_neighbor": int(np.sum(r & (ln > 0))),
        "L_with_R_neighbor": int(np.sum(l & (rn > 0))),
        "RL_enrichment_over_random_placement":
            rl/expected_rl_random if expected_rl_random > 0 else None,
        "contact_components": len(sizes), "largest_contact_component": max(sizes, default=0),
        "mixed_RL_components": int(sum(productive)),
        "scope": "Spatial contacts and retention observations; not a test of an enclosing membrane, functional memory, or organismality.",
    }


def estimate(values) -> dict:
    arr = np.asarray(values, dtype=float)
    n = len(arr)
    mean = float(arr.mean())
    se = float(arr.std(ddof=1)/math.sqrt(n)) if n > 1 else None
    return {"n": n, "mean": mean, "monte_carlo_standard_error": se,
            "approximate_95_percent_normal_interval":
                [mean-1.96*se, mean+1.96*se] if se is not None else None}


def difference_in_differences(rows: list[dict], field: str) -> dict:
    """Between-arm Monte Carlo contrasts; no presumption of paired path coupling."""
    arms = {}
    for p, b in ((True, True), (True, False), (False, True), (False, False)):
        subset = [r[field] for r in rows if r["flags"]["produce_L"] == p
                  and r["flags"]["protect_R"] == b and r["flags"].get("copy_R", True)]
        arms[(p, b)] = estimate(subset)
    weights = {(True, True): 1, (True, False): -1,
               (False, True): -1, (False, False): 1}
    mean = sum(weights[k]*a["mean"] for k, a in arms.items())
    # Seeds can share an initial preparation. This explicitly conservative
    # sum-of-SE bound avoids assuming either positive or zero covariance.
    se_bound = sum(a["monte_carlo_standard_error"] or 0 for a in arms.values())
    return {"field": field, "mean": mean,
            "standard_error_upper_bound_without_covariance_assumption": se_bound,
            "approximate_95_percent_interval_using_SE_bound":
                [mean-1.96*se_bound, mean+1.96*se_bound],
            "interpretation": "Exploratory factorial interaction, not an external objective or reward."}
