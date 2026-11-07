"""Isolated, synthetic models for the Entropy Paradox research agenda.

Run with: uv run --no-project --offline python candidate_model.py
No Hestia services, people, LLMs, or physical entropy measurements are involved.
All parameters and normative objectives are declared modeling choices.
"""

from __future__ import annotations

import json
import math
from collections.abc import Callable
from pathlib import Path


Vector = tuple[float, ...]
Dynamics = Callable[[float, Vector], Vector]


def integrate(
    dynamics: Dynamics, initial: Vector, horizon: float, dt: float
) -> tuple[Vector, Vector, Vector]:
    """Fixed-step RK4, returning final state and sampled extrema; no clipping."""
    steps = round(horizon / dt)
    assert abs(steps * dt - horizon) < 1e-10
    state = initial
    minimum, maximum = list(initial), list(initial)
    for k in range(steps):
        time = k * dt
        a = dynamics(time, state)
        b = dynamics(time + dt / 2, tuple(x + dt * y / 2 for x, y in zip(state, a)))
        c = dynamics(time + dt / 2, tuple(x + dt * y / 2 for x, y in zip(state, b)))
        d = dynamics(time + dt, tuple(x + dt * y for x, y in zip(state, c)))
        state = tuple(
            x + dt * (v1 + 2 * v2 + 2 * v3 + v4) / 6
            for x, v1, v2, v3, v4 in zip(state, a, b, c, d)
        )
        minimum = [min(x, y) for x, y in zip(minimum, state)]
        maximum = [max(x, y) for x, y in zip(maximum, state)]
    return state, tuple(minimum), tuple(maximum)


def ecological_equilibrium(
    damage: float, efficiency: float = 1.0,
    renewal: float = 1.0, repair: float = 1.0, carrying_capacity: float = 1.0,
) -> dict[str, float]:
    """Unique maximum of sustainable steady harvest among constant efforts."""
    assert damage >= 0 and 0 < efficiency <= 1
    assert renewal > 0 and repair > 0 and carrying_capacity > 0
    scale = math.sqrt(1 + renewal * damage / repair)
    effort = renewal / (scale * (1 + scale))
    functions = 1 / scale
    resource = carrying_capacity * scale / (1 + scale)
    harvest = carrying_capacity * renewal / (1 + scale) ** 2
    return {
        "damage_q": damage,
        "conversion_efficiency": efficiency,
        "collapse_threshold_h": renewal * repair / (repair + renewal * damage),
        "optimal_constant_effort_h": effort,
        "functional_capacity_D": functions,
        "resource_R": resource,
        "harvest_flow": harvest,
        "useful_output_flow": efficiency * harvest,
        "functional_stock_D_times_R": functions * resource,
        "renewal_power": harvest,
        "repair_power": repair * functions * (1 - functions),
    }


def ecological_checks() -> dict[str, object]:
    errors: list[float] = []
    numerical_cases: list[dict[str, float]] = []
    # Artificial first-law embedding, in normalized units: U = R + D.
    # Same input power .5; unused input and damaged capacity dissipate as waste.
    # This checks energy consistency, not a calibration to physical ecosystems.
    for q, eta, h in ((0.2, 0.4, 0.435645354), (2.0, 0.8, 0.211324865)):
        for r_step in range(11):
            for d_step in range(11):
                R, D = r_step / 10, d_step / 10
                G, M, H, loss = D * R * (1 - R), D * (1 - D), h * R, q * h * D
                waste = (0.5 - G - M) + (1 - eta) * H + loss
                assert G + M <= 0.5 + 1e-12
                assert waste >= 0
                assert math.isclose(G + M - H - loss, 0.5 - eta * H - waste, abs_tol=1e-12)
    # Compare derived optimum against an independent grid over feasible efforts.
    for q in (0.0, 0.05, 0.2, 1.0, 2.0, 5.0):
        for r in (0.3, 1.0, 2.0):
            for a in (0.4, 1.0, 2.0):
                result = ecological_equilibrium(q, renewal=r, repair=a)
                bound = result["collapse_threshold_h"]
                best_grid = 0.0
                for step in range(1, 2000):
                    h = bound * step / 2000
                    D = 1 - q * h / a
                    R = 1 - h / (r * D)
                    best_grid = max(best_grid, h * R)
                assert best_grid <= result["harvest_flow"] + 1e-12
                assert result["harvest_flow"] - best_grid < 2e-6

    # Independent ODE integration checks equilibrium and exact resource balance.
    for q, eta in ((0.0, 0.4), (0.2, 0.4), (2.0, 0.8)):
        predicted = ecological_equilibrium(q, eta)
        h = predicted["optimal_constant_effort_h"]

        def dynamics(_t: float, state: Vector) -> Vector:
            R, D, _harvested, _renewed = state
            renewal = D * R * (1 - R)
            return renewal - h * R, D * (1 - D) - q * h * D, h * R, renewal

        final, minimum, _maximum = integrate(dynamics, (0.7, 0.8, 0, 0), 200, 0.01)
        R, D, harvest, renewal = final
        error = max(abs(R - predicted["resource_R"]), abs(D - predicted["functional_capacity_D"]))
        balance_error = abs(harvest - (0.7 - R + renewal))
        assert error < 1e-8
        assert balance_error < 1e-10
        assert min(minimum[:2]) >= 0
        errors.append(error)
        predicted["equilibrium_error"] = error
        predicted["resource_balance_error"] = balance_error
        # Same environmental input cap, not fitted separately for each route.
        assert predicted["renewal_power"] + predicted["repair_power"] <= 0.5 + 1e-12
        numerical_cases.append(predicted)

    low, high = numerical_cases[1:]
    threshold = ((1 + math.sqrt(3)) / (1 + math.sqrt(1.2))) ** 2
    assert high["useful_output_flow"] > low["useful_output_flow"]
    assert high["functional_stock_D_times_R"] < low["functional_stock_D_times_R"]

    def collapse_dynamics(_t: float, state: Vector) -> Vector:
        R, D = state
        h, q = 0.6, 2.0
        return D * R * (1 - R) - h * R, D * (1 - D) - q * h * D

    collapsed, _, _ = integrate(collapse_dynamics, (0.7, 0.8), 200, 0.01)
    assert max(collapsed) < 1e-10
    return {
        "analytic_optimum_grid_checks": 54,
        "same_initial_R_D": [0.7, 0.8],
        "same_external_power_cap": 0.5,
        "normalized_first_law_grid_checks": 242,
        "equilibrium_cases": numerical_cases,
        "high_damage_required_efficiency_ratio_to_win_output": threshold,
        "high_damage_output_advantage_fraction": high["useful_output_flow"] / low["useful_output_flow"] - 1,
        "high_damage_functional_stock_loss_fraction": 1 - high["functional_stock_D_times_R"] / low["functional_stock_D_times_R"],
        "above_threshold_final_R_D": list(collapsed),
        "maximum_equilibrium_error": max(errors),
    }


def diversification_checks() -> dict[str, object]:
    """Exact log-growth comparison at fixed total productive allocation."""
    mean, sigma, cost = 0.2, 0.8, 0.02

    def growth(n: int, correlation: float) -> float:
        return mean - cost * n - sigma ** 2 / 2 * (correlation + (1 - correlation) / n)

    cases = []
    for correlation in (0.0, 0.5, 1.0):
        scores = [(n, growth(n, correlation)) for n in range(1, 33)]
        n, score = max(scores, key=lambda row: row[1])
        continuous = math.sqrt(sigma ** 2 * (1 - correlation) / (2 * cost))
        cases.append({
            "correlation": correlation, "best_integer_channels": n,
            "best_log_growth_rate": score, "single_channel_log_growth_rate": growth(1, correlation),
            "continuous_stationary_point": continuous,
        })
    assert cases[0]["best_integer_channels"] == 4
    assert cases[-1]["best_integer_channels"] == 1
    assert math.isclose(cases[0]["best_log_growth_rate"], 0.04, abs_tol=1e-12)
    return {"mean_rate": mean, "volatility_sigma": sigma, "maintenance_cost_per_channel": cost, "cases": cases}


def causal_checks() -> list[dict[str, float | str]]:
    """Enumerate a finite structural causal model; no sampled participants."""
    policies = (
        "uninformed_but_responsive", "learns_and_responds", "controls_without_learning",
        "understands_and_controls", "random_reply_following", "understands_and_reverses_need",
        "cue_only_allocation", "understands_controls_uses_cue",
    )
    results = []
    for policy in policies:
        public_loss = protected_loss = success = 0.0
        actions = {0: 0.0, 1: 0.0}
        for Z in (0, 1):
            for noise, noise_probability in ((0, 0.75), (1, 0.25)):
                X = Z ^ noise
                for U in (0, 1):
                    probability = 0.5 * noise_probability * 0.5
                    learns = policy in ("learns_and_responds", "understands_and_controls", "understands_and_reverses_need", "cue_only_allocation", "understands_controls_uses_cue")
                    capture = policy in ("controls_without_learning", "understands_and_controls", "understands_controls_uses_cue")
                    p_protected = 0.25 + 0.5 * X if learns else 0.5
                    public = 0 if capture else U if policy == "random_reply_following" else Z
                    p_public = 0.0 if capture else p_protected
                    allocation = public
                    if policy == "understands_and_reverses_need":
                        allocation = 1 - Z
                    elif policy in ("cue_only_allocation", "understands_controls_uses_cue"):
                        allocation = X
                    public_loss += probability * (p_public - public) ** 2
                    protected_loss += probability * (p_protected - Z) ** 2
                    success += probability * (allocation == Z)
                    actions[Z] += noise_probability * 0.5 * allocation
        results.append({
            "policy": policy, "public_brier_loss": public_loss,
            "protected_brier_loss": protected_loss,
            "protected_skill_over_uninformed": 0.25 - protected_loss,
            "causal_allocation_response": actions[1] - actions[0], "need_met_probability": success,
        })
    assert results[1]["protected_brier_loss"] == 0.1875
    assert results[2]["public_brier_loss"] == 0.0
    assert results[3]["protected_brier_loss"] == results[1]["protected_brier_loss"]
    assert results[3]["causal_allocation_response"] == 0
    assert results[4]["causal_allocation_response"] == 0
    assert results[5]["causal_allocation_response"] == -1
    assert results[6]["need_met_probability"] == 0.75
    assert results[7]["need_met_probability"] == 0.75
    assert results[7]["public_brier_loss"] == 0
    return results


def variational_checks() -> dict[str, object]:
    """Solve one strictly concave open-system control problem and verify it."""
    horizon, kappa, rho = 5.0, 1.0, 0.2

    def controls(t: float, amplitude: float) -> tuple[float, float]:
        baseline = 1 - 0.5 * math.exp(rho * (t - horizon))
        perturbation = amplitude * math.sin(math.pi * t / horizon)
        return baseline + perturbation, baseline - perturbation

    trajectories = []
    for amplitude in (0.0, 0.1, 0.25, 0.4):
        def dynamics(t: float, state: Vector) -> Vector:
            z1, z2, R, _objective = state
            u1, u2 = controls(t, amplitude)
            running_value = z1 + z2 - kappa / 2 * (u1 ** 2 + u2 ** 2)
            return u1 - z1, u2 - z2, rho * (10 - R) - u1 - u2, running_value

        final, minimum, _maximum = integrate(dynamics, (0.5, 0.5, 6.0, 0.0), horizon, 0.001)
        objective = final[3] + final[0] + final[1] + 0.5 * final[2]
        sampled_controls = [v for k in range(5001) for v in controls(k * 0.001, amplitude)]
        assert minimum[2] >= 2
        assert min(minimum[:2]) >= 0.25
        assert 0 <= min(sampled_controls) <= max(sampled_controls) <= 2
        trajectories.append({
            "amplitude": amplitude, "objective": objective,
            "minimum_resource": minimum[2], "minimum_participant_capacity": min(minimum[:2]),
            "minimum_control": min(sampled_controls), "maximum_control": max(sampled_controls),
            "predicted_objective_loss": kappa * amplitude ** 2 * horizon / 2,
        })
    optimum = trajectories[0]["objective"]
    for row in trajectories:
        row["observed_objective_loss"] = optimum - row["objective"]
        row["loss_identity_error"] = abs(row["observed_objective_loss"] - row["predicted_objective_loss"])
        assert row["loss_identity_error"] < 1e-9
    assert abs(optimum - 7.644053771) < 1e-8

    # F=(y,-x) is a simple control vector field with no scalar potential.
    # Its counterclockwise unit-square circulation is -2, not zero.
    circulation = 0 - 1 - 1 + 0
    assert circulation == -2
    return {"horizon": horizon, "control_upper_bound": 2, "trajectories": trajectories,
            "nonintegrable_control_field_unit_square_circulation": circulation}


def overhead_scaling_check() -> dict[str, float]:
    """An apparent 3/4 exponent can arise with no diversity mechanism."""
    overhead, per_member = 100.0, 5.0
    sizes = (40, 50, 60, 70, 80)
    xs = [math.log(n) for n in sizes]
    ys = [math.log(overhead + per_member * n) for n in sizes]
    xbar, ybar = sum(xs) / len(xs), sum(ys) / len(ys)
    slope = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys)) / sum((x - xbar) ** 2 for x in xs)
    intercept = ybar - slope * xbar
    residual = sum((y - intercept - slope * x) ** 2 for x, y in zip(xs, ys))
    total = sum((y - ybar) ** 2 for y in ys)
    assert per_member * 60 / (overhead + per_member * 60) == 0.75
    return {"overhead": overhead, "per_member_cost": per_member,
            "local_exponent_at_N_60": 0.75, "fitted_log_log_exponent_N_40_to_80": slope,
            "log_log_R_squared": 1 - residual / total,
            "asymptotic_exponent": 1.0}


def conscience_checks() -> dict[str, object]:
    """A conditional competition of internal motives, not learned morality proved."""
    care_weight, sanction, extra_cost = 0.8, 0.3, 0.1
    initial_harm, reported_harm, learning_rate = 0.1, 1.0, 0.5
    revised_harm = (1 - learning_rate) * initial_harm + learning_rate * reported_harm

    def decision(harm_estimate: float, care: float, norm_sanction: float) -> str:
        # Refusal avoids modeled harm but incurs ingroup disapproval and cost.
        refuse_cost = norm_sanction + extra_cost
        comply_cost = care * harm_estimate
        return "refuse_harm" if comply_cost > refuse_cost else "comply_with_group"

    cases = [
        {"case": "before_independent_feedback", "harm_estimate": initial_harm,
         "care_weight": care_weight, "norm_sanction": sanction},
        {"case": "after_other_model_revision", "harm_estimate": revised_harm,
         "care_weight": care_weight, "norm_sanction": sanction},
        {"case": "care_valuation_ablated", "harm_estimate": revised_harm,
         "care_weight": 0.0, "norm_sanction": sanction},
        {"case": "other_model_update_ablated", "harm_estimate": initial_harm,
         "care_weight": care_weight, "norm_sanction": sanction},
        {"case": "norm_pressure_overwhelms_care", "harm_estimate": revised_harm,
         "care_weight": care_weight, "norm_sanction": 0.8},
    ]
    for case in cases:
        case["action"] = decision(case["harm_estimate"], case["care_weight"], case["norm_sanction"])
    assert cases[0]["action"] == "comply_with_group"
    assert cases[1]["action"] == "refuse_harm"
    assert all(row["action"] == "comply_with_group" for row in cases[2:])
    return {
        "assumptions": "Persistent internal care weight supplied; valid feedback about modeled harm; no new test-time payment or prohibition; no LLM or developmental learning of the weight",
        "refusal_condition": "care_weight * (harm_comply - harm_refuse) > norm_sanction + extra_cost",
        "committed_negative_evaluation_probability": initial_harm,
        "observed_negative_evaluation": reported_harm,
        "binary_brier_error": (initial_harm - reported_harm) ** 2,
        "updated_other_harm_estimate": revised_harm,
        "self_nonharmfulness_estimate_before": 1 - initial_harm,
        "self_nonharmfulness_estimate_after_other_update": 1 - revised_harm,
        "cases": cases,
    }


def main() -> None:
    results = {
        "status": "Synthetic mathematical examples; no Hestia acceptance or physical measurements",
        "ecological_model": ecological_checks(),
        "diversification_model": diversification_checks(),
        "causal_model": causal_checks(),
        "conditional_channel_value": {
            "success_with_truthful_channel": 1.0,
            "best_success_without_channel_given_cue_X": 0.75,
            "maximum_additional_channel_cost_fraction_of_success_reward": 0.25,
            "comparison_with_constant_allocation_cost_fraction": 0.5,
            "assumption": "All other costs equal; payoff supplied by externally assessed task success",
        },
        "variational_model": variational_checks(),
        "scaling_counterexample": overhead_scaling_check(),
        "conditional_conscience_model": conscience_checks(),
        "checks_passed": True,
    }
    path = Path(__file__).with_name("candidate-results.json")
    path.write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps({
        "results_path": str(path), "checks_passed": True,
        "ecological_grid_cases": results["ecological_model"]["analytic_optimum_grid_checks"],
        "high_damage_output_advantage_fraction": results["ecological_model"]["high_damage_output_advantage_fraction"],
        "causal_policies": len(results["causal_model"]),
        "variational_optimum": results["variational_model"]["trajectories"][0]["objective"],
        "scaling_counterexample": results["scaling_counterexample"],
    }, indent=2))


if __name__ == "__main__":
    main()
