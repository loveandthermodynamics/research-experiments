"""Outward interval check of one illustrative radial flux-balance example."""
from hashlib import sha256
import json
from pathlib import Path
from mpmath import iv

ROOT = Path(__file__).resolve().parent


def main():
    iv.dps = 50
    a, b, c = iv.mpf(3), iv.mpf(1), iv.mpf(1)
    margin = 4*a**3-27*b*c**2
    assert margin > 0
    g = lambda r: a*r-c-b*r**3
    cases = []
    for lower, upper, stable in [(".347296355333", ".347296355335", False),
                                  ("1.532088886237", "1.532088886239", True)]:
        lo, hi = iv.mpf(lower), iv.mpf(upper)
        if stable:
            assert g(lo) > 0 and g(hi) < 0 and lo > 1
        else:
            assert g(lo) < 0 and g(hi) > 0 and hi < 1
        radius = iv.mpf([lower, upper])
        stationary_derivative = a/radius**2-3*b
        assert stationary_derivative < 0 if stable else stationary_derivative > 0
        cases.append(dict(radius_bracket=[lower, upper], stable_radially=stable,
                          stationary_derivative_enclosure=str(stationary_derivative)))
    result = dict(checks_passed=True, arithmetic="50-decimal outward interval arithmetic",
                  model="Rdot=A/R-C/R^2-BR; illustrative A3,B1,C1 in consistent units",
                  fold_margin_enclosure=str(margin), radii=cases,
                  source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                  scope="Conditional radial flux approximation only. No phase-field derivation, nonspherical stability, spontaneous interface, division, lattice correspondence, biodiversity or individual is certified.")
    target = ROOT/"radius-certificate.json"
    encoded = json.dumps(result, indent=2)+"\n"
    if target.exists():
        assert target.read_text() == encoded
    else:
        target.write_text(encoded)
    print(encoded)


if __name__ == "__main__":
    main()
