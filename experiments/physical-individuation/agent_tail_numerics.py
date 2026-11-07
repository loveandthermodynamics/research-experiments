"""Reproduce agent_tail_bridge.md's finite gate, including clipped-tail bins.

The Gamma mass below the first edge is FOLDED INTO the first representative;
mass above the last edge is folded into the last representative. This differs
from dropping both tails and normalizing the retained bins. Although the lower
mass is tiny, stationary length bias divides it by a rate near 1e-9.

Run: python3 -W error agent_tail_numerics.py
Prints JSON to stdout; writes no result or existing document files.
"""

from __future__ import annotations

import json

import numpy as np
from scipy.linalg import eigh
from scipy.special import gammainc, gammaincc


def main() -> None:
    alpha = 1.5
    tau = 1.0
    entry_rate = 0.2
    activation_rate = 1.0
    theta = 0.05
    bins = 192
    edges = np.geomspace(1e-9, 100.0, bins + 1)
    rates = np.sqrt(edges[:-1] * edges[1:])

    # Use the survival-function difference in high bins, where subtracting two
    # CDF values very near one would lose positive Gamma mass to roundoff.
    scaled = tau * edges
    retained_weights = np.where(
        scaled[:-1] < alpha,
        np.diff(gammainc(alpha, scaled)),
        gammaincc(alpha, scaled[:-1]) - gammaincc(alpha, scaled[1:]),
    )
    lower_tail = float(gammainc(alpha, scaled[0]))
    upper_tail = float(gammaincc(alpha, scaled[-1]))
    weights = retained_weights.copy()
    weights[0] += lower_tail
    weights[-1] += upper_tail
    weights /= weights.sum()
    assert np.all(weights > 0)
    assert abs(float(weights.sum()) - 1) < 1e-14

    mean_duration = float(np.sum(weights / rates))
    pi_healthy = 1 / (1 + entry_rate * mean_duration)
    pi = np.r_[pi_healthy, pi_healthy * entry_rate * weights / rates]
    well_energy = theta * np.log(rates / (entry_rate * weights))
    db_error = float(np.max(np.abs(
        np.log(pi[0]) + np.log(entry_rate * weights)
        - np.log(pi[1:]) - np.log(rates)
    )))
    assert db_error < 1e-12

    times = np.array([100.0, 1000.0, 10000.0])

    def first_activation_survival(killing: float):
        # diag(sqrt(pi)) (Q - killing*1_H) diag(1/sqrt(pi))
        # is symmetric. From initial pi its exponential-mode weights are
        # squares of eigenvector projections onto sqrt(pi).
        symmetric = np.zeros((bins + 1, bins + 1))
        symmetric[0, 0] = -entry_rate - killing
        symmetric[np.arange(1, bins + 1), np.arange(1, bins + 1)] = -rates
        symmetric[0, 1:] = np.sqrt(entry_rate * weights * rates)
        symmetric[1:, 0] = symmetric[0, 1:]
        eigenvalues, eigenvectors = eigh(symmetric)
        coefficients = np.sum(eigenvectors * np.sqrt(pi)[:, None], axis=0) ** 2
        assert np.max(eigenvalues) < 0
        survival = np.array([
            np.sum(coefficients * np.exp(eigenvalues * time))
            for time in times
        ])
        return survival, float(-eigenvalues[-1])

    single, ultimate_rate = first_activation_survival(activation_rate)
    common, common_ultimate_rate = first_activation_survival(3 * activation_rate)
    independent = single ** 3
    entry_survival = np.sum(np.exp(-times[:, None] * rates) * weights, axis=1)
    residual_survival = np.sum(
        np.exp(-times[:, None] * rates) * (weights / rates), axis=1
    ) / mean_duration
    entry_ideal = (1 + times / tau) ** (-alpha)
    residual_ideal = (1 + times / tau) ** (-(alpha - 1))
    assert np.all(independent <= common)

    dropped = retained_weights / retained_weights.sum()
    output = {
        "construction": "Gamma bin masses, with both omitted tails folded into endpoint representatives",
        "parameters": {
            "alpha": alpha, "tau": tau, "entry_rate": entry_rate,
            "activation_rate": activation_rate, "theta": theta,
            "bins": bins, "edge_min": float(edges[0]), "edge_max": float(edges[-1]),
        },
        "folded_lower_tail_probability": lower_tail,
        "folded_upper_tail_probability": upper_tail,
        "lower_tail_contribution_to_mean_duration": lower_tail / float(rates[0]),
        "mean_entry_duration": mean_duration,
        "mean_entry_duration_if_tails_dropped_and_renormalized": float(np.sum(dropped / rates)),
        "stationary_healthy_probability": pi_healthy,
        "stationary_mean_gate_energy": float(np.sum(pi[1:] * well_energy)),
        "max_detailed_balance_log_residual": db_error,
        "ultimate_single_exponential_rate": ultimate_rate,
        "ultimate_common_exponential_rate": common_ultimate_rate,
        "rows": [
            {
                "time": float(time), "entry_survival": float(entry_survival[i]),
                "ideal_entry_survival": float(entry_ideal[i]),
                "residual_survival": float(residual_survival[i]),
                "ideal_residual_survival": float(residual_ideal[i]),
                "single_first_activation_survival": float(single[i]),
                "three_independent_first_activation_survival": float(independent[i]),
                "three_common_first_activation_survival": float(common[i]),
            }
            for i, time in enumerate(times)
        ],
        "apparent_exponents": {
            "single": (-np.diff(np.log(single)) / np.diff(np.log(times))).tolist(),
            "three_independent": (-np.diff(np.log(independent)) / np.diff(np.log(times))).tolist(),
            "three_common": (-np.diff(np.log(common)) / np.diff(np.log(times))).tolist(),
        },
    }
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
