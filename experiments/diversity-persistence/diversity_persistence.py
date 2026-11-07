"""Exact illustrative checks for diversity and persistence, not a cosmic theorem.

No diversity reward is present. Population or maintained capacity is multiplied
by declared environmental returns. Constant composition between generations
requires regeneration of that composition; it is a model premise.

The hazard illustration separately assumes discrete cycles, with repair from a
surviving module before the next exposure, retained templates for every response
type, and sufficient paid repair time and material. Its productive power is net
of declared maintenance/repair expenditure. Each feasible module must remain
above minimum viable size. This effective model does not resolve repair kinetics
or derive biological species, individuality, or thermodynamic entropy.
"""

import itertools
import json
import math
from pathlib import Path


def binomial_log_growth(n: int, low: float, high: float) -> float:
    """Independent equal-probability channels with fixed equal regeneration."""
    return sum(
        math.comb(n, k) * math.log((k * high + (n - k) * low) / n)
        for k in range(n + 1)
    ) / (2**n)


def shock_survival(n: int, q: float, common: float) -> float:
    """One-step survival of a substitutable, repairable function."""
    return (1 - common) * (1 - q**n)


def net_work_before_ruin(n: int, common: float) -> float:
    """Unit cycles, work 100-10n each cycle, shock at end, local failure .2."""
    power = 100 - 10 * n
    collapse_probability = common + (1 - common) * 0.2**n
    return power / collapse_probability


def main() -> None:
    # Two equally likely environmental states, complementary responses.
    low, high = 0.5, 1.8
    mono_log = (math.log(low) + math.log(high)) / 2
    mix_log = math.log((low + high) / 2)
    assert mono_log < 0 < mix_log
    assert math.isclose(math.exp(mono_log), math.sqrt(0.9))
    assert math.isclose(math.exp(mix_log), 1.15)

    # Enumerated independent returns: no covariance approximation required.
    growth = [binomial_log_growth(n, low, high) for n in range(1, 21)]
    assert all(b > a for a, b in itertools.pairwise(growth))
    assert all(g < math.log((low + high) / 2) for g in growth)
    # Perfectly shared responses erase the diversification advantage.
    correlated_growth = [mono_log] * 20
    assert len(set(correlated_growth)) == 1

    # Exact enumeration checks q**n, with no simulation sampling error.
    for n in range(1, 9):
        q = 0.2
        survive = sum(
            math.prod(q if failed else 1 - q for failed in state)
            for state in itertools.product((False, True), repeat=n)
            if not all(state)
        )
        assert math.isclose(survive, shock_survival(n, q, 0), abs_tol=1e-12)

    # Fixed positive common-mode catastrophe creates a reliability ceiling.
    c = 0.01
    assert shock_survival(100, 0.2, c) <= 1 - c
    # Series dependence has the opposite sign from interchangeable backup.
    assert (1 - 0.2) ** 8 < (1 - 0.2) ** 2
    # Maintenance changes unlimited expansion into a finite optimum.
    optima = {
        str(c): max(range(1, 10), key=lambda n: net_work_before_ruin(n, c))
        for c in (0.0, 0.01)
    }
    assert optima == {"0.0": 9, "0.01": 4}

    # Without repair, more last-survivor time need not create more total work.
    # Split fixed stock B equally over n domains with independent failure clocks
    # lambda_local and common clock lambda_common. Expected active stock is
    # B exp(-(lambda_local+lambda_common)t) for EVERY n. Its time integral is
    # B/(lambda_local+lambda_common), independent of n. With no common shock,
    # expected last-survivor time is H_n/lambda_local, which DOES grow with n.
    last_survivor_times = [sum(1 / j for j in range(1, n + 1))
                           for n in range(1, 21)]
    assert all(b > a for a, b in itertools.pairwise(last_survivor_times))
    capacity, local_rate, common_rate = 100.0, 1.0, 0.1
    no_repair_integrated_capacity = capacity / (local_rate + common_rate)

    # Differential replication yields a convex finite-interval exposure.
    # N(T)/N(0)=exp(mu*T)*cosh(a*x*T); no ex-post max oracle is used.
    # g''(x)=a*a*T/cosh(a*x*T)**2 > 0 for finite x.
    # This exponential-growth approximation needs fuel and no active crowding
    # over the stated interval; it does not predict unbounded physical growth.
    a, interval, mu = 0.4, 1.0, 0.02
    convex_growth = [mu + math.log(math.cosh(a * s * interval)) / interval
                     for s in (0.0, 0.5, 1.0, 2.0)]
    assert all(b > a for a, b in itertools.pairwise(convex_growth))

    # Fat-tailed LOG losses with bounded gross gains. At each independent time,
    # Z_i ~ Pareto(z0, alpha), R_i=exp(r-cost*n-Z_i), M_n=mean_i R_i.
    # The minimum Z has Pareto exponent n*alpha. Exact pointwise bounds:
    # r-cost*n-min(Z)-log(n) <= log M_n <= r-cost*n-min(Z).
    # E[min(Z)] exists iff n*alpha>1 and is z0*n*alpha/(n*alpha-1).
    # For n*alpha<=1, time-average log growth is -infinity (no finite-copy
    # extinction barrier is included). Independent response channels can move
    # this model from catastrophic log-loss tails to finite positive growth.
    alpha, z0, r, cost = 0.5, 0.1, 2.0, 0.01
    n = 3
    min_loss_mean = z0 * n * alpha / (n * alpha - 1)
    fat_tail_lower_bound = r - cost * n - min_loss_mean - math.log(n)
    fat_tail_upper_bound = r - cost * n - min_loss_mean
    assert fat_tail_lower_bound > 0
    assert all(k * alpha <= 1 for k in (1, 2))
    for loss_vector in itertools.product((z0, 0.2, 1.0, 10.0, 1000.0), repeat=3):
        smallest = min(loss_vector)
        # Stable log-sum-exp avoids underflow for severe log losses.
        log_multiplier = r - cost * n - smallest + math.log(
            sum(math.exp(smallest - z) for z in loss_vector) / n
        )
        assert r - cost * n - smallest - math.log(n) <= log_multiplier + 1e-12
        assert log_multiplier <= r - cost * n - smallest + 1e-12

    result = {
        "status": "conditional analytical illustrations, not a biodiversity optimum",
        "no_diversity_reward": True,
        "complementary_two_state": {
            "common_arithmetic_mean": (low + high) / 2,
            "monoculture_log_growth": mono_log,
            "monoculture_geometric_multiplier": math.exp(mono_log),
            "regenerated_mix_log_growth": mix_log,
            "regenerated_mix_geometric_multiplier": math.exp(mix_log),
        },
        "independent_channels_log_growth_n_1_to_20": growth,
        "perfect_common_response_log_growth": mono_log,
        "independent_failure_formula": "Pr(complete loss)=q**n",
        "common_failure_formula": "Pr(complete loss)=c+(1-c)*q**n",
        "repair_and_cost_assumptions": __doc__,
        "work_before_ruin_optimal_n": optima,
        "no_repair_counterexample": {
            "fixed_total_stock": capacity,
            "integrated_expected_active_stock_all_n": no_repair_integrated_capacity,
            "last_survivor_times_no_common_clock_n_1_to_20": last_survivor_times,
            "meaning": "Continuation alone does not prove more cumulative organization.",
        },
        "convex_replication_growth_at_spreads_0_05_1_2": convex_growth,
        "fat_tailed_log_loss_construction": {
            "alpha": alpha, "pareto_lower_scale": z0,
            "baseline_log_return": r, "per_channel_log_cost": cost,
            "n_1_and_2_time_average_log_growth": "minus infinity",
            "n_3_expected_log_growth_lower_bound": fat_tail_lower_bound,
            "n_3_expected_log_growth_upper_bound": fat_tail_upper_bound,
            "supremum_gross_return_n_3": math.exp(r - cost * n - z0),
            "independent_time_and_channel_exposures": True,
            "minimum_log_loss_tail_exponent": n * alpha,
            "caveats": [
                "Equal allocation is restored at each time by assumed regeneration.",
                "Continuous capacity; no finite-copy extinction barrier.",
                "A shared infinite-mean log-loss shock defeats every finite n.",
                "Rates and response distributions are premises, not microscopic derivations.",
                "This is not a result about thermodynamic entropy or unlimited biodiversity.",
            ],
        },
        "limits": [
            "Species richness is not failure-response diversity.",
            "No stored exergy or physical negentropy identity is derived.",
            "Fat-tailed losses can defeat unrestricted pooling.",
            "Finite populations, repair dynamics and inheritance remain unresolved.",
        ],
        "checks_passed": True,
    }
    Path(__file__).with_name("diversity-persistence-results.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
