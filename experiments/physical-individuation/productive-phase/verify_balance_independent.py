"""Independent exact balance and local-cycle certificate; no simulation.

This does not import an implementation of the productive spatial generator.
Rational arithmetic checks a concrete local cycle; high precision evaluates the
stated thermal-source bound. Irreducibility/positive-NESS reasoning is a
structural proof under the explicitly recorded hypotheses, not enumeration of
the exponentially large lattice state space.
"""
from __future__ import annotations

import copy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import mpmath as mp
import sympy as s

ROOT = Path(__file__).resolve().parent


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def hamiltonian(material, templates):
    side = 4
    energy = F(2)*material.count(2) + F(6)*sum(r != 0 for r in templates)
    for i, value in enumerate(material):
        row, col = divmod(i, side)
        for j in (row*side+(col+1) % side, ((row+1) % side)*side+col):
            energy -= F(12, 5)*int(value == 2 and material[j] == 2)
        if value == 2:
            energy -= {0:F(0), 1:F(-5, 2), 2:F(5, 2)}[templates[i]]
    return energy


def explicit_cycle():
    material, templates = [0]*16, [0]*16
    material[1] = material[6] = 2
    material[5] = 1
    templates[4], templates[6] = 1, 2
    initial = (material[:], templates[:])
    initial_h = hamiltonian(material, templates)
    target, x_parent, y_parent = 5, 4, 6
    assert x_parent != y_parent and target not in (x_parent, y_parent)
    mu, theta = F(12), F(1)
    steps, cx, cy, losses = [], 0, 0, 0
    h = initial_h
    specifications = (("X_copy", 1, 2, 1, mu, F(1, 10)*F(19, 20)),
                      ("X_loss", 0, 2, 0, F(0), F(1, 50)),
                      ("Y_copy", 2, 1, 2, mu, F(1, 10)*F(19, 20)),
                      ("Y_loss", 0, 1, 0, F(0), F(1, 50)))
    for name, role, new_material, new_template, work, mobility in specifications:
        before = dict(material=material[:], templates=templates[:])
        if role:
            assert templates[target] == 0
            assert templates[x_parent if role == 1 else y_parent] == role
            assert material[target] == (1 if role == 1 else 2)
            cx += int(role == 1)
            cy += int(role == 2)
        else:
            assert templates[target] in (1, 2)
            losses += 1
        material[target], templates[target] = new_material, new_template
        new_h = hamiltonian(material, templates)
        delta = new_h-h
        h = new_h
        assert material.count(2)-initial[0].count(2) == cx-cy
        assert sum(r != 0 for r in templates)-sum(r != 0 for r in initial[1]) == cx+cy-losses
        affinity = (work-delta)/theta
        steps.append(dict(name=name, before=before, after=dict(material=material[:], templates=templates[:]),
                          delta_H_exact=str(delta), chemical_work_exact=str(work),
                          heat_to_bath_exact=str(work-delta), reverse_symmetric_mobility_exact=str(mobility),
                          log_forward_reverse_ratio_exact=str(affinity)))
    assert (material, templates) == initial
    assert h == initial_h
    total_h = sum(F(step["delta_H_exact"]) for step in steps)
    total_w = sum(F(step["chemical_work_exact"]) for step in steps)
    affinity = sum(F(step["log_forward_reverse_ratio_exact"]) for step in steps)
    assert total_h == 0 and total_w == 24 and affinity == 24
    return dict(side=4, working_carriers=3, target=target, X_parent=x_parent, Y_parent=y_parent,
                exact_initial_final_H=str(initial_h), physical_arrays_return_exactly=True,
                copy_family_counts=dict(X=1, Y=1), passive_losses=2,
                total_delta_H_exact=str(total_h), total_work_exact=str(total_w),
                cycle_log_affinity_exact=str(affinity), formula="2*mu/theta",
                steps=steps, scope="Finite rational local-cycle example; no stochastic trajectory or phase-formation evidence.")


def algebra():
    jx, jy, j0, loss, birth = s.symbols("jX jY j0 loss birth", real=True)
    equations = [jx-jy+j0, jx+jy+birth-loss]
    solution = s.solve(equations, [jx, jy])
    assert s.simplify(solution[jx]-(loss-birth-j0)/2) == 0
    assert s.simplify(solution[jy]-(loss-birth+j0)/2) == 0
    K, M, delta, kappa, c0, diff, mean_n = s.symbols("K M delta kappa c0 diff mean_N", positive=True)
    limiting_equality = delta*mean_n-(kappa*M+delta*(c0*(K-mean_n)+diff*M))
    bound = s.solve(limiting_equality, mean_n)[0]
    expected = (kappa*M/delta+c0*K+diff*M)/(1+c0)
    assert s.simplify(bound-expected) == 0
    return dict(stationary_currents={"JX":str(solution[jx]), "JY":str(solution[jy])},
                tightened_one_family_expectation_bound=str(expected),
                finite_T_addition="(M+E[N_initial])/(delta*(1+c0)*T)",
                power_identity="P=mu*(delta*E[N]-E[I_intensity])",
                source_condition="E_r>=max(0,g_r) gives exact downhill loss delta; otherwise use saturated births and actual loss lower bound",
                c_difference="diff=max(cB-c0,0) for conservative arbitrary-affinity bound")


def numbers():
    mp.mp.dps = 70
    mp.iv.dps = 60
    c0 = 2*mp.exp(-6)
    cb = mp.exp(mp.mpf("-8.5"))+mp.exp(mp.mpf("-3.5"))
    delta, K, M = mp.mpf(".02"), 256, 77
    results = {}
    for label, kappa in (("no_background", mp.mpf(0)), ("weak_background", mp.mpf(".0002"))):
        bound = (kappa*M/delta+c0*K+(cb-c0)*M)/(1+c0)
        iv_c0 = 2*mp.iv.exp(-mp.iv.mpf(6))
        iv_cb = mp.iv.exp(mp.iv.mpf("-8.5"))+mp.iv.exp(mp.iv.mpf("-3.5"))
        iv_kappa = mp.iv.mpf("0" if label == "no_background" else ".0002")
        iv_bound = (iv_kappa*M/mp.iv.mpf(".02")+iv_c0*K+(iv_cb-iv_c0)*M)/(1+iv_c0)
        assert 0 < bound < K
        results[label] = dict(kappa=mp.nstr(kappa, 10), expected_population_upper_bound=mp.nstr(bound, 60),
                              outward_interval_enclosing_the_bound=str(iv_bound),
                              density_upper_bound=mp.nstr(bound/K, 60),
                              mean_source_upper_at_all_D=mp.nstr(delta*(c0*(K-M)+cb*M), 60))
    # The thermal-source inequality must hold for every admissible occupation,
    # not only the fixed point or an observed subset of states.
    cases = 0
    for n in range(257):
        for db in range(min(77, 256-n)+1):
            exact = delta*(c0*(256-n)+(cb-c0)*db)
            upper = delta*(c0*(256-n)+(cb-c0)*77)
            assert exact <= upper
            cases += 1
    return dict(c0=mp.nstr(c0, 60), cB=mp.nstr(cb, 60), cB_ge_c0=bool(cb >= c0),
                lattice_sites=K, working_carriers=M, delta=mp.nstr(delta, 10),
                cases_thermal_source_bound=cases, numerical_values=results,
                forward_X_max_cost_exact="21/2", forward_Y_max_cost_exact="68/5",
                forward_Y_acceptance_at_maximum=mp.nstr(mp.exp(mp.mpf("-1.6")), 60),
                numerical_scope="70-digit evaluation plus60-digit outward interval enclosure of each displayed bound; no spatial result.")


def certificate():
    return dict(checks_passed=True, exact_cycle=explicit_cycle(), symbolic_algebra=algebra(), thermal_and_energy_values=numbers(),
        positive_NESS_theorem=dict(
            hypotheses=["finite physical occupancy state space at fixed material inventory", "positive finite reverse rates and thermal D<->X/Y", "connected local material transport", "at least one working carrier and a target with two distinct neighboring helpers", "both copy activities enabled, theta>0 and mu>0", "kappa=0 and no other fueled bypass"],
            reasoning="Thermal template transitions and copy-plus-loss material flips connect the physical space. The four-step cycle has affinity2mu/theta. Irreducibility makes stationary pi positive, so zero entropy production would force detailed balance on all paired probability fluxes, contradicting this cycle. Since dE[B]/dt=0, JX=JY; since sigma=P/theta and P=mu(JX+JY), both are strictly positive.",
            probability_fluxes="F=pi(s)*q(s,s'); reverse=pi(s')*q_rev; sigma=sum_unordered_pairs(F-Frev)*log(F/Frev)",
            one_activity_effective_equilibrium=dict(X_only="H-mu*B", Y_only="H+mu*B"),
            not_certified_by_state_enumeration=True,
            limitations=["Currents can be arbitrarily rare", "Positive thermal sources are essential to stated irreducibility", "Journal counters/instance IDs are observational and excluded from thermodynamic state", "No high population, ancestry persistence, phase formation, size selection or taxonomic necessity follows"]),
        source_sha256={name:digest(ROOT/name) for name in ("balance-and-cost.md", "model-spec.md", "generalist-comparison.md", "verify_generalist.py")},
        verifier_sha256=digest(__file__),
        scope="Independent exact local-cycle and balance review before model implementation or ensembles. No trajectory, reward, phase selector or old evidence modification.")


if __name__ == "__main__":
    result = certificate()
    (ROOT/"balance-independent-certificate.json").write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"checks_passed":True, "cycle_affinity":result["exact_cycle"]["cycle_log_affinity_exact"],
                      "bounds":result["thermal_and_energy_values"]["numerical_values"]}, indent=2))
