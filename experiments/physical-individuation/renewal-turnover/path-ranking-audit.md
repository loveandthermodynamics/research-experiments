# The variational bridge must include kinetics

For a continuous-time Markov jump history gamma with initial probability p0,

    log P[gamma] = log p0(x0) + sum_jumps log k(x,y) - integral_0^T lambda(x_t) dt.

Here lambda(x) is the total escape rate. Write reversible rates as
k(x,y) = a(x,y) exp[s(x,y)/2], where a is symmetric and s is antisymmetric.
The entropy contribution is one part of log probability. Symmetric mobilities
and waiting-time terms also enter. They cancel in an appropriate comparison
with the reversed history, but need not cancel between two forward histories.
This is consistent with [Crooks's fluctuation theorem](https://arxiv.org/abs/cond-mat/9901352)
and the dynamical-activity contribution studied by
[Maes](https://arxiv.org/abs/1705.07412).

## Exact counterexample

The executable `path_ranking_counterexample.py` defines two three-state cycles
sharing state 0: 0-1-2-0 and 0-3-4-0. All internal energies are equal. Cycle A
has forward/reverse rate ratio exp(2) on each edge and symmetric mobility
0.001. Cycle B has ratio exp(1) and mobility 1. Chemostatted chemical affinities
supply the corresponding work and reservoir entropy; this is an explicitly
driven process, not a closed equilibrium machine.

Incoming and outgoing rates balance at every state, so the stationary
distribution is uniform. Compare histories starting at 0, making one forward
jump at time 0.5, then waiting until time 1. The A-history produces 2 kB of
total entropy; the B-history produces 1 kB. Their system-entropy boundary terms
are zero under the uniform stationary ensemble. Nevertheless the higher-entropy
history has a much lower path density. The same ordering holds for the finite
events with exactly one such jump at any time in the interval. The script
integrates these event probabilities analytically, checks stationarity and
verifies both forward/reverse fluctuation ratios. There is no sampling error.

Thus greater entropy production cannot universally rank forward probabilities.
The mobility of the high-affinity route can approach zero without changing
its entropy per event. A proposed selection principle must supply and retain
these physical kinetic terms.

## A recent proposed shortcut does not close the gap

[Segal's 2026 preprint](https://arxiv.org/abs/2603.15230) explicitly introduces an
ansatz comparing distinct forward histories by exp[(sigma1-sigma2)/2]. That
assumption does not follow from a forward/reverse fluctuation relation; the
counterexample above violates it. The preprint also restricts linearly
improving replication rates to an early regime while using them for an
asymptotic dominance claim. These are reasons the paper does not furnish the
universal derivation sought here. This is an audit of the stated argument,
not a claim that adaptation or dissipation-driven organization cannot occur.

## What remains available

A full path-action or large-deviation formulation remains legitimate once the
generator, driving and admissible histories are specified. For a finite
irreducible time-homogeneous generator, the empirical occupation/flow rate
functional has the form

    I(rho,q) = sum_(x,y) [q_xy log(q_xy/(rho_x k_xy)) - q_xy + rho_x k_xy],

subject to normalized nonnegative occupation, nonnegative permitted flows
and flow conservation. At stationary occupation and typical flows q_xy =
rho_x k_xy, every summand vanishes. Each is nonnegative by the convexity of
u log u - u + 1. This conditional least-rate-functional statement includes
the entire generator. It does not select biodiversity or individuation until
the local chemistry actually generates, maintains and differentiates the
relevant organizations. A favorable model result would establish that bridge
under its premises, rather than remove the premises.
