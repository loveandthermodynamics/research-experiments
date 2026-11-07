"""Independent numerical audit of the proposed open R,L maintenance threshold.

Equations are reconstructed directly here; no study simulation is imported.
This verifies one explicitly supplied driven chemistry, not organismhood.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

ROOT = Path(__file__).resolve().parent
M, ALPHA, ETA, AFFINITY = 0.95, 4.0, math.exp(-20), 40.0


def limiting_material(lvalue: float) -> float:
    return lvalue + math.exp(-ALPHA * lvalue) + lvalue * math.exp(ALPHA * lvalue)


def vector(y: np.ndarray, produce: bool = True, protect: bool = True) -> np.ndarray:
    r, ell = y
    d = M - r - ell
    barrier = math.exp(-ALPHA * ell) if protect else 1.0
    return np.array([r*d - ETA*r*r - barrier*r + ETA*barrier*d,
                     (r*d - ETA*r*ell if produce else 0.) - ell + ETA*d])


def jacobian(y: np.ndarray) -> np.ndarray:
    r, ell = y
    d, h = M-r-ell, math.exp(-ALPHA*ell)
    return np.array([[d-r-2*ETA*r-(1+ETA)*h,
                      -r+ALPHA*h*r-ETA*h*(1+ALPHA*d)],
                     [d-r-ETA*ell-ETA, -(1+ETA)*(1+r)]])


def l_from_r(r: float) -> float:
    return (r+ETA)*(M-r)/((1+ETA)*(1+r))


def thermodynamics(y: np.ndarray) -> dict:
    r, ell = y
    d, h = M-r-ell, math.exp(-ALPHA*ell)
    forward = np.array([r*d, r*d, h*r, ell])
    reverse = np.array([ETA*r*r, ETA*r*ell, ETA*h*d, ETA*d])
    currents = forward-reverse
    entropy_production = float(np.sum(currents*np.log(forward/reverse)))
    power = AFFINITY*float(currents[0]+currents[1])
    chemical_potentials = np.array([20+math.log(r/d), 20+math.log(ell/d)])
    free_energy_change = float(np.sum(chemical_potentials*vector(y)))
    error = abs(free_energy_change-(power-entropy_production))
    assert entropy_production >= 0
    assert error < 1e-11
    return dict(currents=currents.tolist(), chemical_power=power,
                entropy_production_rate=entropy_production,
                free_energy_change_rate=free_energy_change, balance_error=error)


def main() -> None:
    fold = minimize_scalar(limiting_material, bounds=(0., .5), method="bounded",
                           options={"xatol": 1e-15})
    assert fold.success and fold.fun < M < 1
    lroots = [brentq(lambda ell: limiting_material(ell)-M, a, b, xtol=1e-15)
              for a, b in [(0., fold.x), (fold.x, M)]]
    limiting = [dict(R=ell*math.exp(ALPHA*ell), L=ell,
                     D=math.exp(-ALPHA*ell)) for ell in lroots]
    grid = np.unique(np.r_[np.geomspace(1e-14, 1e-3, 1000), np.linspace(1e-3, M, 10000)])
    brackets = [(a, b) for a, b in zip(grid[:-1], grid[1:])
                if vector(np.array([a, l_from_r(a)]))[0] * vector(np.array([b, l_from_r(b)]))[0] < 0]
    roots = [brentq(lambda r: vector(np.array([r, l_from_r(r)]))[0], a, b, xtol=1e-15)
             for a, b in brackets]
    assert len(roots) == 3
    equilibria = []
    for r in roots:
        y = np.array([r, l_from_r(r)])
        eig = np.linalg.eigvals(jacobian(y))
        assert np.max(np.abs(vector(y))) < 1e-12
        step = 1e-7
        numerical = np.column_stack([(vector(y+step*np.eye(2)[j])-vector(y-step*np.eye(2)[j]))/(2*step)
                                     for j in range(2)])
        assert np.max(np.abs(numerical-jacobian(y))) < 1e-7
        equilibria.append(dict(R=float(y[0]), L=float(y[1]), D=float(M-y.sum()),
                               eigenvalues=eig.tolist(), stable=bool(np.max(eig.real) < 0),
                               thermodynamics=thermodynamics(y)))
    assert [e["stable"] for e in equilibria] == [True, False, True]
    trajectories = []
    for initial, target in [([.01, .01], 0), ([.30, .15], 2)]:
        solved = solve_ivp(lambda t, y: vector(y), (0, 500), initial,
                           rtol=1e-10, atol=1e-12, dense_output=True)
        assert solved.success
        observations = solved.sol(np.linspace(0, 500, 2001))
        assert float(observations.min()) > 0 and float(observations.sum(axis=0).max()) < M
        expected = np.array([equilibria[target][k] for k in ("R", "L")])
        assert np.max(np.abs(solved.y[:, -1]-expected)) < 1e-8
        trajectories.append(dict(initial=initial, final=solved.y[:, -1].tolist(),
                                 reaches_equilibrium_index=target))
    balances = [thermodynamics(np.array(y))["balance_error"] for y in
                [(.1, .1), (.2, .3), (.01, .01), (.5, .2)]]
    interventions = []
    high = [equilibria[2][k] for k in ("R", "L")]
    for produce, protect in [(False, True), (True, False), (False, False)]:
        solved = solve_ivp(lambda t, y: vector(y, produce, protect), (0, 500), high,
                           rtol=1e-10, atol=1e-12)
        assert solved.success and solved.y[0, -1] < 1e-6
        interventions.append(dict(produce_L=produce, protect_R=protect,
                                  initial=high, final=solved.y[:, -1].tolist()))
    no_protection_bound = ETA*M/(1-M)
    no_production_hmin = math.exp(-ALPHA*ETA*M)
    assert no_production_hmin > M
    no_production_bound = ETA*M/(no_production_hmin-M)
    result = dict(checks_passed=True, parameters=dict(M=M, alpha=ALPHA, reverse_coefficient=ETA,
                  theta=1, epsilon_R=20, epsilon_L=20, fuel_affinity=AFFINITY),
                  limiting_fold=dict(L=float(fold.x), M_min=float(fold.fun)),
                  limiting_positive_equilibria=limiting, reversible_equilibria=equilibria,
                  basin_examples=trajectories, maximum_off_stationary_balance_error=max(balances),
                  interventions_from_high_state=interventions,
                  stationary_R_upper_bounds=dict(no_protection=no_protection_bound,
                                                  no_production=no_production_bound),
                  supplied_premises=["R-catalyzed copying and L production", "L-dependent symmetric passive R barrier", "material conservation", "fixed fuel/waste chemostats with counted work"],
                  scope="Conditional kinetic saddle-node and bistability with continually paid turnover. No spatial boundary, arbitrary-small-founder formation, new trait alphabet, biological aging, recursive individuality or moral legitimacy follows.")
    (ROOT/"threshold-verification.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
