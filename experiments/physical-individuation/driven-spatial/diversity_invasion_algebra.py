"""Algebraic check of a weak-response invasion coefficient; no trajectories.

Two declared inventories give positive and adverse cases of the same R,L law.
This does not certify a finite response amplitude, coexistence, or availability.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def bisect(function, left: float, right: float) -> float:
    fl = function(left)
    assert fl * function(right) <= 0
    for _ in range(140):
        midpoint = (left + right) / 2
        fm = function(midpoint)
        if fl * fm <= 0:
            right = midpoint
        else:
            left, fl = midpoint, fm
    return (left + right) / 2


def example(material: float) -> dict:
    # a=b=delta=gamma=1; alpha=4. Gate flips either way at kappa=.25.
    alpha, kappa = 4.0, 0.25
    function = lambda ell: ell + math.exp(-alpha * ell) + ell * math.exp(alpha * ell)
    derivative = lambda ell: (1 - alpha * math.exp(-alpha * ell)
                              + math.exp(alpha * ell) * (1 + alpha * ell))
    fold = bisect(derivative, 0.0, 1.0)
    ell = bisect(lambda v: function(v) - material, fold, material)
    precursor = math.exp(-alpha * ell)
    active = material - ell - precursor
    h = math.exp(-alpha * ell)
    flux = active * precursor
    j11, j12 = -active, active * (alpha * precursor - 1)
    j21, j22 = precursor - active, -active - 1
    trace = j11 + j22
    determinant = j11 * j22 - j12 * j21
    assert trace < 0 and determinant > 0
    q = 2 * kappa
    denominator = q*q - q*trace + determinant
    cov_r = -flux * (q - j22) / denominator
    cov_l = -flux * j21 / denominator
    # The covariance must solve (2*kappa I - J)c = (-J_copy, 0).
    residuals = [(q-j11)*cov_r - j12*cov_l + flux,
                 -j21*cov_r + (q-j22)*cov_l]
    assert max(map(abs, residuals)) < 1e-14
    coefficient_from_covariance = -2*alpha*h*cov_l
    coefficient_closed_form = 2*alpha*h*flux*j21 / denominator
    assert math.isclose(coefficient_from_covariance, coefficient_closed_form,
                        rel_tol=1e-14, abs_tol=1e-15)
    assert abs(active*precursor-ell) < 1e-14
    return dict(M=material, alpha=alpha, kappa=kappa,
                R=active, L=ell, D=precursor, h=h, paid_copy_flux=flux,
                jacobian=[[j11, j12], [j21, j22]], trace=trace,
                determinant=determinant, telegraph_resolvent_denominator=denominator,
                first_order_covariance_R_gate=cov_r,
                first_order_covariance_L_gate=cov_l,
                covariance_equation_residual=max(map(abs, residuals)),
                opposite_response_invasion_coefficient=coefficient_closed_form,
                general_continuous_trait_coefficient=coefficient_closed_form/2,
                sign="positive" if coefficient_closed_form > 0 else "negative")


def main() -> None:
    cases = [example(0.95), example(1.5)]
    assert [case["sign"] for case in cases] == ["positive", "negative"]
    result = dict(checks_passed=True, calculations="Deterministic scalar roots and 2x2 resolvent only",
                  cases=cases,
                  conclusion="The weak-forcing coefficient changes sign with D-R. No diversity reward enters.",
                  scope="Negligible-source concentration model near a stable renewing high state. Positive coefficient is an asymptotic invasion result conditional on resident persistence and a regular small-amplitude stationary expansion. No finite-amplitude simulation, global stochastic coexistence, finite-copy-number survival, spatial localization, or availability advantage is established by this script.")
    destination = ROOT / "diversity-invasion-algebra.json"
    assert not destination.exists(), "refusing to replace evidence"
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
