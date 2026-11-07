"""Outward-interval certificate for three finite-reverse steady states.

mpmath interval arithmetic encloses all operations at 50 decimal digits.
IVT signs certify existence; interval Jacobians certify local stability.
No claim that these are the only possible equilibria is needed here.
"""
from __future__ import annotations

import json
from pathlib import Path
from mpmath import iv

ROOT = Path(__file__).resolve().parent
iv.dps = 50
M, A, ETA = iv.mpf("0.95"), iv.mpf(4), iv.exp(-iv.mpf(20))


def material(r):
    ell = (r+ETA)*(M-r)/((1+ETA)*(1+r))
    return ell, M-r-ell


def scalar(r):
    ell, d = material(r)
    h = iv.exp(-A*ell)
    return r*d-ETA*r*r-h*r+ETA*h*d


def main():
    intervals = [
        ("3.91619e-8", "3.91621e-8", True),
        ("0.034562985656", "0.034562985659", False),
        ("0.232847617045", "0.232847617048", True),
    ]
    certificates = []
    for lo, hi, stable in intervals:
        lower, upper, r = iv.mpf(lo), iv.mpf(hi), iv.mpf([lo, hi])
        flo, fhi = scalar(lower), scalar(upper)
        assert flo*fhi < 0
        ell, d = material(r)
        assert r > 0 and ell > 0 and d > 0
        h = iv.exp(-A*ell)
        j11 = d-r-2*ETA*r-(1+ETA)*h
        j12 = -r+A*h*r-ETA*h*(1+A*d)
        j21 = d-r-ETA*ell-ETA
        j22 = -(1+ETA)*(1+r)
        trace, determinant = j11+j22, j11*j22-j12*j21
        assert trace < 0
        assert determinant > 0 if stable else determinant < 0
        # Scalar derivative along g=0 is det(J)/g_L, which has definite sign.
        assert j22 < 0
        certificates.append(dict(R_bracket=[lo, hi], L_enclosure=str(ell), D_enclosure=str(d),
                                  scalar_at_lower=str(flo), scalar_at_upper=str(fhi),
                                  jacobian_trace_enclosure=str(trace),
                                  jacobian_determinant_enclosure=str(determinant),
                                  unique_root_in_bracket=True, locally_asymptotically_stable=stable))
    hmin = iv.exp(-A*ETA*M)
    assert hmin > M
    no_protection_bound = ETA*M/(1-M)
    no_production_bound = ETA*M/(hmin-M)
    assert iv.mpf(intervals[2][0]) > no_protection_bound
    assert iv.mpf(intervals[2][0]) > no_production_bound
    result = dict(checks_passed=True, arithmetic="mpmath.iv outward intervals, 50 decimal digits",
                  method="IVT endpoint signs, derivative sign from det(J)/g_L, trace/determinant stability",
                  finite_reverse_coefficient_enclosure=str(ETA), equilibria=certificates,
                  no_protection_stationary_R_upper_bound=str(no_protection_bound),
                  no_production_stationary_R_upper_bound=str(no_production_bound),
                  scope="Three physical finite-reverse equilibria exist with certified local stability. The high equilibrium is impossible at steady state under either specified channel removal. This proves conditional reciprocal maintenance in a supplied open wellmixed chemistry, not spatial individuation, biodiversity or morality.")
    (ROOT/"threshold-interval-certificate.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
