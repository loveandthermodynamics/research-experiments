"""Exact finite-state counterexample to entropy-only forward-path ranking.

No simulation or fitted parameter. Two driven three-edge cycles share state 0.
Uniform stationary probabilities make system-entropy boundary changes vanish.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main() -> None:
    size, horizon, jump_time = 5, 1.0, 0.5
    rates = [[0.0] * size for _ in range(size)]
    cycles = (([0, 1, 2], 2.0, 0.001), ([0, 3, 4], 1.0, 1.0))
    for vertices, entropy_per_step, mobility in cycles:
        for i, j in zip(vertices, vertices[1:] + vertices[:1]):
            rates[i][j] = mobility * math.exp(entropy_per_step / 2)
            rates[j][i] = mobility * math.exp(-entropy_per_step / 2)
    escape = [sum(row) for row in rates]
    stationarity_error = max(abs(sum(rates[j][i] for j in range(size)) - escape[i]) / size
                             for i in range(size))
    assert stationarity_error < 1e-14

    def path_density(i: int, j: int, t: float) -> float:
        return math.exp(-math.log(size) - escape[i] * t
                        + math.log(rates[i][j]) - escape[j] * (horizon - t))

    def one_jump_probability(i: int, j: int) -> float:
        # Integrate the exact density over jump time, including no second jump.
        difference = escape[i] - escape[j]
        integral = (-math.expm1(-difference * horizon) / difference
                    if difference else horizon)
        return rates[i][j] * math.exp(-escape[j] * horizon) * integral / size

    rows = []
    for target, entropy in ((1, 2.0), (3, 1.0)):
        density = path_density(0, target, jump_time)
        reversed_density = path_density(target, 0, horizon - jump_time)
        probability = one_jump_probability(0, target)
        reversed_probability = one_jump_probability(target, 0)
        assert math.isclose(math.log(density / reversed_density), entropy, abs_tol=1e-13)
        assert math.isclose(math.log(probability / reversed_probability), entropy, abs_tol=1e-13)
        rows.append(dict(target=target, total_entropy_production_in_kB=entropy,
                         forward_path_density=density, reverse_path_density=reversed_density,
                         forward_one_jump_event_probability=probability,
                         reverse_one_jump_event_probability=reversed_probability))
    high, low = rows
    assert high["forward_path_density"] < low["forward_path_density"]
    assert high["forward_one_jump_event_probability"] < low["forward_one_jump_event_probability"]
    result = dict(checks_passed=True, states=size, rates=rates, escape_rates=escape,
                  stationary_distribution=[1 / size] * size,
                  stationarity_error=stationarity_error, horizon=horizon, jump_time=jump_time,
                  paths=rows,
                  actual_high_to_low_density_ratio=high["forward_path_density"] / low["forward_path_density"],
                  actual_high_to_low_event_probability_ratio=high["forward_one_jump_event_probability"] / low["forward_one_jump_event_probability"],
                  entropy_only_density_ratio=math.exp((2.0 - 1.0) / 2),
                  scope="Driven CTMC with equal internal energies and two chemostatted reaction affinities. Reservoir entropy is counted. Uniform stationary boundary probabilities yield total path entropy equal to jump log-rate ratio. This disproves universal entropy-only ranking of forward histories, not fluctuation relations or conditional biological mechanisms.")
    (ROOT / "path-ranking-counterexample.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k not in ("rates", "escape_rates")}))


if __name__ == "__main__":
    main()
