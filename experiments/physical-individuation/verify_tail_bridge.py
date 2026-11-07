"""Independent finite reversible-gate check of the latency theorem."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.special import gammainc, gammaincc


def main():
    alpha, tau, entry, thermal, nu, n = 1.5, 1., .2, .05, 1., 3
    edges = np.geomspace(1e-9, 100., 193)
    rates = np.sqrt(edges[:-1] * edges[1:])
    lower = gammainc(alpha, tau * edges)
    upper = gammaincc(alpha, tau * edges)
    weights = np.where(lower[:-1] > .5, upper[:-1] - upper[1:], np.diff(lower))
    assert np.all(weights > 0)
    excluded_mass = 1. - float(weights.sum())
    weights /= weights.sum()
    mean_duration = float((weights / rates).sum())
    pi_healthy = 1. / (1. + entry * mean_duration)
    pi = np.r_[pi_healthy, pi_healthy * entry * weights / rates]
    wells = np.r_[0., thermal * np.log(rates / (entry * weights))]
    generator = np.zeros((len(pi), len(pi)))
    generator[0, 1:] = entry * weights
    generator[1:, 0] = rates
    np.fill_diagonal(generator, -generator.sum(axis=1))
    db_error = float(np.max(np.abs(np.log(entry * weights / rates)
                                  + (wells[1:] - wells[0]) / thermal)))
    stationarity_error = float(np.max(np.abs(pi @ generator)))
    times = np.array([100., 1000., 10000.])

    def killed_survival(hazard, query_times):
        killed = generator.copy()
        killed[0, 0] -= hazard
        symmetric = np.sqrt(pi[:, None]) * killed / np.sqrt(pi[None, :])
        assert np.max(np.abs(symmetric - symmetric.T)) < 1e-12
        eigenvalues, eigenvectors = eigh(symmetric)
        coefficients = (eigenvectors.T @ np.sqrt(pi)) ** 2
        assert eigenvalues[-1] < 0
        survival = np.exp(np.outer(query_times, eigenvalues)) @ coefficients
        return survival, float(-eigenvalues[-1])

    single, minimum_decay = killed_survival(nu, times)
    common, _ = killed_survival(n * nu, times)
    independent = single ** n
    query = np.geomspace(.001, 1e5, 200)
    test_single, _ = killed_survival(nu, query)
    test_common, _ = killed_survival(n * nu, query)
    assert np.all(test_common >= test_single ** n - 1e-12)
    residual_weights = weights / rates / mean_duration
    residual = np.exp(-np.outer(times, rates)) @ residual_weights
    ideal = (1 + times / tau) ** -(alpha - 1)
    exponent = lambda values: float(-np.log(values[-1] / values[-2]) / np.log(10.))
    eta = .03 * np.exp(-.8 / thermal)
    output = {
        "status": "passed",
        "tail_convention": "Omit both Gamma tails and renormalize; agent_tail_numerics.py separately checks folding into endpoint bins.",
        "gate_states": len(pi),
        "parameters": {"alpha": alpha, "tau": tau, "entry_rate": entry,
                       "thermal": thermal, "paid_hazard_per_healthy_site": nu, "sites_and_gates": n},
        "pi_healthy": pi_healthy,
        "mean_entry_duration": mean_duration,
        "expected_gate_energy_one_unit": float(pi @ wells),
        "expected_gate_energy_both_matched_comparisons": float(n * (pi @ wells)),
        "omitted_gamma_probability_before_renormalization": excluded_mass,
        "detailed_balance_log_residual": db_error,
        "stationarity_residual": stationarity_error,
        "rows": [{"time": float(t), "single_activation_survival": float(s),
                  "independent_three_activation_survival": float(i),
                  "common_three_activation_survival": float(c),
                  "finite_residual_vs_ideal_relative_error": float(r / ideal_r - 1)}
                 for t, s, i, c, r, ideal_r in zip(times, single, independent, common, residual, ideal)],
        "apparent_tail_exponents_1000_to_10000": {
            "single": exponent(single), "independent_three": exponent(independent), "common_three": exponent(common)},
        "eventual_smallest_decay_rate_single": minimum_decay,
        "jensen_order_checked_at_200_times": True,
        "passive_thermal_rate_per_site": float(eta),
        "passive_factor_at_time_10000": float(np.exp(-n * eta * times[-1])),
        "scope": "First-activation latency only; separate retained-template activation model, finite-mixture cutoff, no post-activation exergy result.",
    }
    root = Path(__file__).resolve().parent
    (root / "tail-bridge-verification.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
