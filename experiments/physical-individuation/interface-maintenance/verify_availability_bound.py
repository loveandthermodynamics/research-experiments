"""Independent eight-state check of the conditional availability theorem.

The positive analytical lower bound is derived in theorem-audit.md. Numerical
matrix evolution checks the finite-temperature example; it is not a replacement
for the analytical coupling/entropy proof or an interval-arithmetic certificate.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import simpson
from scipy.sparse.linalg import expm_multiply


def generator(common, theta=.05, epsilon=1., e=2., lam=1., mu=1.):
    states = [(mask, fuel) for mask in range(4) for fuel in range(2)]
    index = {s: i for i, s in enumerate(states)}
    q = np.zeros((8, 8))
    def add(src, dest, rate):
        i, j = index[src], index[dest]
        q[i, j] += rate
        q[i, i] -= rate
    for mask, fuel in states:
        if common:
            if mask:
                add((mask, fuel), (0, fuel), lam)
            else:
                for target in (1, 2, 3):
                    add((0, fuel), (target, fuel), lam*math.exp(-target.bit_count()*epsilon/theta))
        else:
            for bit in (1, 2):
                active = bool(mask & bit)
                add((mask, fuel), (mask ^ bit, fuel), lam if active else lam*math.exp(-epsilon/theta))
        if mask in (1, 2) and fuel:
            add((mask, fuel), (3, 0), mu*fuel)
        if mask == 3 and fuel == 0:
            for target in (1, 2):
                add((mask, 0), (target, 1), mu*math.exp(-(e-epsilon)/theta))
    pi = np.array([math.exp(-(epsilon*mask.bit_count()+e*fuel)/theta) for mask, fuel in states])
    pi /= pi.sum()
    assert np.max(np.abs(pi[:, None]*q-(pi[:, None]*q).T)) < 1e-14
    return states, q


def comparison(points):
    H, theta, epsilon, e = 4., .05, 1., 2.
    curves = []
    for common in (False, True):
        states, q = generator(common)
        initial = np.zeros(8)
        initial[states.index((3, 1))] = 1
        p = expm_multiply(q.T, initial, start=0, stop=H, num=points)
        assert p.min() > -1e-12 and np.max(np.abs(p.sum(axis=1)-1)) < 1e-11
        marginal = p.reshape(points, 4, 2).sum(axis=2)
        mean = np.sum(marginal*np.array([0., 1., 1., 2.]), axis=1)
        entropy = -np.sum(np.where(marginal > 0, marginal*np.log(np.maximum(marginal, 1e-300)), 0), axis=1)
        availability = epsilon*mean-theta*entropy+2*theta*math.log1p(math.exp(-epsilon/theta))
        curves.append(availability)
    return float(simpson(curves[0]-curves[1], x=np.linspace(0, H, points)))


def main():
    H, theta = 4., .05
    active_advantage = .5-2*math.exp(-H)+(1.5+H)*math.exp(-2*H)
    li = 4*math.exp(-1/theta)
    lc = 4*math.exp(-1/theta)+math.exp(-2/theta)
    correction = H*H*(li+lc)
    entropy_penalty = 2*theta*H*math.log(2)
    bound = active_advantage-correction-entropy_penalty
    coarse, fine = comparison(1001), comparison(4001)
    assert bound > 0 and fine > bound and abs(fine-coarse) < 2e-7
    result = {"model": "eight states: two binary sites, one finite repair packet",
              "zero_temperature_active_time_advantage": active_advantage,
              "reverse_event_coupling_correction": correction, "entropy_penalty": entropy_penalty,
              "analytical_lower_bound": bound, "matrix_integrated_availability_advantage": fine,
              "quadrature_refinement_change": abs(fine-coarse),
              "scope": "Conditional finite-fuel regeneration theorem. Does not prove spatial individuality, universal biodiversity maximization, or moral legitimacy."}
    Path(__file__).with_name("availability-bound-verification.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
