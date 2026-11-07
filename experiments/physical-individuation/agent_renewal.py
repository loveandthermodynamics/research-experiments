"""Finite-fuel, paid-repair benchmark. Prints results; writes no artifacts.

Two labeled chemical sites have active/inactive states. Inactive sites retain
their type template. One surviving site can catalyze repair of the other, using
one fuel molecule. Diversity changes failure dependence, not per-site marginal
damage rates. No instantaneous resetting, diversity reward, or reproduction
operator is present.

The irreversible limit has exact renewal formulas. The finite-temperature
extension is closed to matter, coupled to a thermal bath, and detailed balanced.
First total inactivity is an episode endpoint, not irreversible extinction in
that extension. Linear state rewards below measure active energy, not exact
nonequilibrium exergy. Actual active-subsystem availability is also computed
from the evolving probability distribution over a finite declared horizon.

Run: python3 agent_renewal.py  (requires numpy and scipy)
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, replace

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.sparse.linalg import expm_multiply


@dataclass(frozen=True)
class Model:
    fuel_packets: int = 3
    active_energy: float = 1.0
    fuel_gap: float = 2.0
    damage_rate: float = 1.0
    repair_parameter: float = 1.0
    repair_law: str = "constant"
    damage_dependence: str = "independent"
    added_common_rate: float = 0.0
    thermal_energy: float | None = None

    def repair_rate(self, fuel: int) -> float:
        if fuel == 0:
            return 0.0
        if self.repair_law == "constant":
            # A saturated effective rate, paired edge-symmetrically below.
            return self.repair_parameter
        assert self.repair_law == "mass_action"
        return self.repair_parameter * fuel

    def damage_rates(self) -> tuple[float, float]:
        if self.damage_dependence == "independent":
            return self.damage_rate, 0.0
        if self.damage_dependence == "common":
            return 0.0, self.damage_rate
        assert self.damage_dependence == "mixed"
        assert self.added_common_rate >= 0
        return self.damage_rate, self.added_common_rate


def build(model: Model):
    assert model.fuel_gap >= model.active_energy > 0
    assert model.damage_rate > 0 and model.repair_parameter >= 0
    assert model.thermal_energy is None or model.thermal_energy > 0
    states = [(mask, f) for mask in range(4)
              for f in range(model.fuel_packets + 1)]
    index = {state: i for i, state in enumerate(states)}
    count = len(states)
    q = np.zeros((count, count))
    repair_plus = np.zeros(count)
    repair_minus = np.zeros(count)
    independent_rate, common_rate = model.damage_rates()

    def edge(source, target, rate, kind="passive"):
        i, j = index[source], index[target]
        assert rate >= 0 and i != j
        q[i, j] += rate
        if kind == "forward_repair":
            repair_plus[i] += rate
        elif kind == "reverse_repair":
            repair_minus[i] += rate

    for mask, fuel in states:
        if independent_rate:
            for site in range(2):
                if mask & (1 << site):
                    edge((mask, fuel), (mask ^ (1 << site), fuel),
                         independent_rate)
                elif model.thermal_energy is not None:
                    edge((mask, fuel), (mask | (1 << site), fuel),
                         independent_rate * math.exp(
                             -model.active_energy / model.thermal_energy))
        if common_rate:
            if mask:
                # One common clock disables all currently active sites.
                edge((mask, fuel), (0, fuel), common_rate)
            elif model.thermal_energy is not None:
                # Paired thermal reversals of those common-mode transitions.
                for target in (1, 2, 3):
                    edge((0, fuel), (target, fuel), common_rate * math.exp(
                        -model.active_energy * target.bit_count()
                        / model.thermal_energy))

        for site in range(2):
            helper = 1 - site
            if not (mask & (1 << helper)):
                continue
            if not (mask & (1 << site)) and fuel:
                edge((mask, fuel), (3, fuel - 1), model.repair_rate(fuel),
                     "forward_repair")
            elif (mask & (1 << site)) and fuel < model.fuel_packets \
                    and model.thermal_energy is not None:
                # Paired forward edge has fuel+1. The mobility coefficient is
                # chosen on that edge, rather than at two different endpoints.
                mobility = model.repair_rate(fuel + 1) / (fuel + 1)
                reverse = mobility * (model.fuel_packets - fuel) * math.exp(
                    -(model.fuel_gap - model.active_energy)
                    / model.thermal_energy)
                edge((mask, fuel), (mask ^ (1 << site), fuel + 1), reverse,
                     "reverse_repair")

    q[np.diag_indices(count)] = -q.sum(axis=1)
    assert np.max(np.abs(q.sum(axis=1))) < 1e-13
    energy = np.array([model.active_energy * mask.bit_count()
                       + model.fuel_gap * f for mask, f in states])
    return states, index, q, energy, repair_plus, repair_minus


def renewal_prediction(model: Model) -> tuple[float, float]:
    """Exact irreversible recurrence; allows fuel-dependent repair rates."""
    independent_rate, common_rate = model.damage_rates()
    marginal_rate = independent_rate + common_rate
    reach_one = 2 * independent_rate / (2 * independent_rate + common_rate)
    extra = 0.0
    expected_repairs = 0.0
    for fuel in range(1, model.fuel_packets + 1):
        mu = model.repair_rate(fuel)
        eta = reach_one * mu / (marginal_rate + mu)
        extra = eta * extra + model.active_energy / marginal_rate * eta
        expected_repairs = eta * (1 + expected_repairs)
    return 2 * model.active_energy / marginal_rate + extra, expected_repairs


def hitting_summary(model: Model) -> dict[str, float]:
    states, index, q, energy, plus, minus = build(model)
    live = np.array([i for i, (mask, _) in enumerate(states) if mask])
    dead = np.array([i for i, (mask, _) in enumerate(states) if not mask])
    a = -q[np.ix_(live, live)]
    start = list(live).index(index[(3, model.fuel_packets)])

    def occupation(reward):
        return float(np.linalg.solve(a, np.asarray(reward)[live])[start])

    active = np.array([model.active_energy * mask.bit_count()
                       for mask, _ in states])
    area = occupation(active)
    repairs = occupation(plus)
    reversed_repairs = occupation(minus)
    # Integrate bath heat released on transitions, including the terminal jump.
    heat_rate = -np.sum(q * energy[None, :], axis=1)
    heat = occupation(heat_rate)
    final_fuel_reward = np.sum(
        q[:, dead] * np.array([states[i][1] for i in dead])[None, :], axis=1)
    final_fuel = occupation(final_fuel_reward)
    assert abs(final_fuel - (model.fuel_packets - repairs + reversed_repairs)) < 1e-10
    ledger = 2 * model.active_energy + model.fuel_gap * (repairs - reversed_repairs)
    assert abs(heat - ledger) < 1e-10

    if model.thermal_energy is None:
        predicted, predicted_repairs = renewal_prediction(model)
        assert abs(area - predicted) < 1e-11
        assert abs(repairs - predicted_repairs) < 1e-11
        marginal_rate = sum(model.damage_rates())
        assert area <= (2 * model.active_energy
                        + model.fuel_gap * model.fuel_packets) / marginal_rate + 1e-12

    return {"first_inactivity_mean_time": occupation(np.ones(len(states))),
            "integrated_active_energy_until_first_inactivity": area,
            "expected_forward_repairs": repairs,
            "expected_reverse_repairs": reversed_repairs,
            "expected_terminal_fuel": final_fuel,
            "expected_net_heat_to_bath": heat,
            "first_law_error": abs(heat - ledger)}


def reference(model: Model, states):
    theta = model.thermal_energy
    assert theta is not None
    logz = 2 * math.log1p(math.exp(-model.active_energy / theta)) \
        + model.fuel_packets * math.log1p(math.exp(-model.fuel_gap / theta))
    return np.array([
        -(model.active_energy * mask.bit_count() + model.fuel_gap * f) / theta
        + math.lgamma(model.fuel_packets + 1) - math.lgamma(f + 1)
        - math.lgamma(model.fuel_packets - f + 1) - logz
        for mask, f in states])


def finite_temperature_summary(model: Model, horizon: float = 4.0):
    states, index, q, _, _, _ = build(model)
    theta = model.thermal_energy
    assert theta is not None
    logpi = reference(model, states)
    error = 0.0
    for i in range(len(states)):
        for j in range(i + 1, len(states)):
            if q[i, j] or q[j, i]:
                assert q[i, j] > 0 and q[j, i] > 0
                error = max(error, abs(logpi[i] + math.log(q[i, j])
                                       - logpi[j] - math.log(q[j, i])))
    assert error < 1e-11
    initial = np.zeros(len(states))
    initial[index[(3, model.fuel_packets)]] = 1
    times = np.linspace(0, horizon, 1001)
    ps = expm_multiply(q.T, initial, start=0, stop=horizon, num=len(times))
    assert ps.min() > -1e-12
    assert np.max(abs(ps.sum(axis=1) - 1)) < 1e-11
    active_ps = np.zeros((len(times), 4))
    for i, (mask, _) in enumerate(states):
        active_ps[:, mask] += ps[:, i]
    logza = 2 * math.log1p(math.exp(-model.active_energy / theta))
    logpia = np.array([-model.active_energy * mask.bit_count() / theta - logza
                       for mask in range(4)])
    safe = np.maximum(active_ps, np.finfo(float).tiny)
    ba = theta * np.sum(active_ps * (np.log(safe) - logpia), axis=1)
    assert ba.min() > -1e-12
    reward = model.active_energy * np.array([mask.bit_count() for mask in range(4)])
    def availability_at(time):
        distribution = expm_multiply(q.T * time, initial)
        marginal = np.zeros(4)
        for i, (mask, _) in enumerate(states):
            marginal[mask] += distribution[i]
        return float(theta * np.sum(marginal * (
            np.log(np.maximum(marginal, np.finfo(float).tiny)) - logpia)))

    availability_area, quadrature_error = quad(
        availability_at, 0, horizon, epsabs=1e-10, epsrel=1e-10)
    augmented = np.zeros((len(states) + 1, len(states) + 1))
    augmented[:-1, :-1] = q.T
    augmented[-1, :-1] = [reward[mask] for mask, _ in states]
    augmented_initial = np.append(initial, 0)
    active_energy_area = float(expm_multiply(
        augmented * horizon, augmented_initial)[-1])
    return {"horizon": horizon, "max_log_detailed_balance_error": error,
            "integral_actual_active_availability": availability_area,
            "availability_quadrature_error_estimate": quadrature_error,
            "integral_active_energy": active_energy_area,
            "terminal_actual_active_availability": float(ba[-1]),
            "reference": "same thermal Gibbs measure for both damage dependencies"}


def active_energy_integral(model: Model, horizon: float) -> float:
    states, index, q, _, _, _ = build(model)
    augmented = np.zeros((len(states) + 1, len(states) + 1))
    augmented[:-1, :-1] = q.T
    augmented[-1, :-1] = [model.active_energy * mask.bit_count()
                          for mask, _ in states]
    initial = np.zeros(len(states) + 1)
    initial[index[(3, model.fuel_packets)]] = 1
    return float(expm_multiply(augmented * horizon, initial)[-1])


def coupling_checks():
    """A rate-tradeoff optimum, not an autonomous evolutionary derivation.

    Independent damage lambda=1; new common hazard gamma*g with gamma=1;
    finite-time paid repair mu=alpha*g with alpha=4. One fuel packet.
    Every comparison has the same initial stored structure and fuel.
    """
    def formula(g):
        return (4 + 30*g + 10*g*g) / (2 + 13*g + 16*g*g + 5*g**3)

    def derivative_numerator(g):
        return 8 - 88*g - 410*g*g - 300*g**3 - 50*g**4

    optimum = brentq(derivative_numerator, 0, 1, xtol=1e-14)
    checks = []
    for g in list(np.linspace(0, 1, 21)) + [optimum]:
        model = Model(fuel_packets=1, repair_parameter=4*g,
                      damage_dependence="mixed", added_common_rate=g)
        result = hitting_summary(model)
        error = abs(result["integrated_active_energy_until_first_inactivity"]
                    - formula(g))
        assert error < 1e-11
        checks.append(error)
    model = Model(fuel_packets=1, repair_parameter=4*optimum,
                  damage_dependence="mixed", added_common_rate=optimum)
    result = hitting_summary(model)
    assert formula(optimum) > formula(0) > formula(1)
    return {"coupling_optimum": optimum,
            "integrated_active_energy_at_optimum": formula(optimum),
            "isolated_value": formula(0), "fully_coupled_value": formula(1),
            "initial_material": "two fixed sites and one conserved fuel/waste packet",
            "initial_stored_energy_every_case": 4.0,
            "same_fuel_packets_every_case": 1,
            "max_formula_vs_generator_error": max(checks),
            "optimum_resource_ledger": result,
            "scope": "trait-performance optimum; no heredity or discovery of optimum modeled"}


def availability_convergence_checks(horizon: float = 4.0):
    """A positive actual availability margin with an independent entropy bound.

    For two sites B_A=epsilon E[a]-theta H(p_A)+theta log Z_A.
    The entropy-difference term for two designs has magnitude <=theta log(4).
    A finite-time semigroup perturbation bounds the thermal energy correction.
    """
    independent0 = Model()
    common0 = Model(damage_dependence="common")
    zero_energy_gain = (active_energy_integral(independent0, horizon)
                        - active_energy_integral(common0, horizon))
    q0_independent = build(independent0)[2]
    q0_common = build(common0)[2]
    rows = []
    for theta in (.1, .05, .025):
        independent = replace(independent0, thermal_energy=theta)
        common = replace(common0, thermal_energy=theta)
        left = finite_temperature_summary(independent, horizon)
        right = finite_temperature_summary(common, horizon)
        eta_independent = float(np.max(np.sum(
            abs(build(independent)[2] - q0_independent), axis=1)))
        eta_common = float(np.max(np.sum(
            abs(build(common)[2] - q0_common), axis=1)))
        # Conservative bound: reward range <=2epsilon and generator row norm eta.
        thermal_correction_bound = independent.active_energy * horizon**2 \
            * (eta_independent + eta_common)
        entropy_correction_bound = theta * horizon * math.log(4)
        lower_bound = zero_energy_gain - thermal_correction_bound \
            - entropy_correction_bound
        actual_gain = left["integral_actual_active_availability"] \
            - right["integral_actual_active_availability"]
        assert lower_bound > 0
        assert actual_gain >= lower_bound - 1e-8
        rows.append({"thermal_energy": theta,
                     "actual_active_availability_gain": actual_gain,
                     "entropy_correction_bound": entropy_correction_bound,
                     "thermal_semigroup_correction_bound": thermal_correction_bound,
                     "derived_positive_gain_lower_bound": lower_bound,
                     "sum_quadrature_error_estimates": left[
                         "availability_quadrature_error_estimate"] + right[
                             "availability_quadrature_error_estimate"]})
    return {"horizon": horizon,
            "irreversible_finite_horizon_active_energy_gain": zero_energy_gain,
            "finite_temperature_cases": rows,
            "scope": "finite-horizon Gibbs availability, common reference; numerical matrix exponentials and quadrature"}


def main():
    # More fuel is never free: exact energy accounting is checked for every run.
    for law in ("constant", "mass_action"):
        for fuel in range(9):
            hitting_summary(Model(fuel_packets=fuel, repair_law=law))
            hitting_summary(Model(fuel_packets=fuel, repair_law=law,
                                  damage_dependence="common"))
    results = {}
    for dependence in ("independent", "common"):
        base = Model(damage_dependence=dependence)
        finite = Model(damage_dependence=dependence, thermal_energy=.05)
        results[dependence] = {"irreversible": hitting_summary(base),
                               "finite_temperature": hitting_summary(finite),
                               "actual_availability": finite_temperature_summary(finite)}
    assert abs(results["independent"]["irreversible"][
        "integrated_active_energy_until_first_inactivity"] - 2.875) < 1e-12
    assert results["independent"]["actual_availability"][
        "integral_actual_active_availability"] > results["common"][
            "actual_availability"]["integral_actual_active_availability"]
    print(json.dumps({"checks_passed": True, "results": results,
                      "coupling_optimum": coupling_checks(),
                      "availability_convergence": availability_convergence_checks(),
                      "scope": "finite fuel, paid repair, explicit bath; fixed templates; no new individuality or negentropy identity"}, indent=2))


if __name__ == "__main__":
    main()
