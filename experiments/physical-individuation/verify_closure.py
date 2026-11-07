"""Independent ODE implementation and implicit-integrator cross-check."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


def main():
    root = Path(__file__).resolve().parent
    data = json.loads((root / "agent-closure-results.json").read_text())
    rows = {}
    for name, case in data["cases"].items():
        p = case["parameters"]
        n, eps, e, theta = (p[k] for k in ("material", "active_energy", "fuel_energy", "thermal"))
        delta, k, f, w = (p[key] for key in ("decay", "copy_coefficient", "fuel", "waste"))
        eta = delta * np.exp(-eps / theta)
        kr = k * np.exp(-(e - eps) / theta)
        gap = e + theta * np.log(f / w)

        def energy(xy):
            c = np.r_[xy, n - sum(xy)]
            return float(eps * sum(xy) + theta * np.sum(c * (np.log(c) - 1)))

        def ode(_, state):
            x, y = state[:2]
            d = n - x - y
            assert min(x, y, d) > 0
            plus = np.array([p["channel_x"] * k * f * x * y * d,
                             p["channel_y"] * k * f * x * y * d, delta * x, delta * y])
            minus = np.array([p["channel_x"] * kr * w * x * x * y,
                              p["channel_y"] * kr * w * x * y * y, eta * d, eta * d])
            net = plus - minus
            derivative = np.array([net[0] - net[2], net[1] - net[3]])
            affinities = np.array([gap - eps + theta * np.log(d / x),
                                   gap - eps + theta * np.log(d / y),
                                   eps + theta * np.log(x / d), eps + theta * np.log(y / d)])
            work = gap * sum(net[:2])
            dissipation = float(net @ affinities)
            assert dissipation >= -1e-9
            heat = e * sum(net[:2]) - eps * sum(derivative)
            return np.r_[derivative, work, dissipation, heat, sum(net[:2])]

        initial = np.r_[case["initial"], np.zeros(4)]
        sol = solve_ivp(ode, (0., case["horizon"]), initial, method="Radau",
                        rtol=1e-10, atol=1e-12, max_step=1.)
        assert sol.success
        final = sol.y[:2, -1]
        work, dissipation, heat, net_fuel = sol.y[2:, -1]
        endpoint_error = float(np.max(np.abs(final - case["final"])))
        g_error = abs(energy(final) - energy(initial[:2]) - work + dissipation)
        e_error = abs(eps * sum(final - initial[:2]) - e * net_fuel + heat)
        assert endpoint_error < 2e-7 and g_error < 2e-7 and e_error < 2e-7
        rows[name] = {"endpoint_difference_from_DOP853": endpoint_error,
                      "independent_G_work_balance_error": float(g_error),
                      "independent_energy_heat_balance_error": float(e_error),
                      "steps": len(sol.t), "final": final.tolist()}
    p = data["parameters"]
    eta = p["decay"] * np.exp(-p["active_energy"] / p["thermal"])
    r = p["copy_coefficient"] * p["waste"] * np.exp(-(p["fuel_energy"] - p["active_energy"]) / p["thermal"])
    q = p["copy_coefficient"] * p["fuel"]
    roots = sorted(z.real for z in np.roots([2*q+r, -q*p["material"], p["decay"]+2*eta, -eta*p["material"]]) if abs(z.imag) < 1e-10)
    assert np.max(np.abs(np.array(roots) - [row["u"] for row in data["symmetric_fixed_points"]])) < 1e-12
    output = {"status": "passed", "independent_integrator": "Radau versus source DOP853",
              "fixed_points_independently_recomputed": roots, "cases": rows,
              "sha256": {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                         for name in ("agent_closure.py", "agent-closure-results.json")}}
    (root / "closure-verification.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": "passed", "cases": len(rows),
                      "max_endpoint_difference": max(row["endpoint_difference_from_DOP853"] for row in rows.values()),
                      "max_G_work_error": max(row["independent_G_work_balance_error"] for row in rows.values())}, indent=2))


if __name__ == "__main__":
    main()
