# A finite fuel stock is a persistence assay, not a sustained-flow cosmos

The current lattice is a finite autonomous system. Its material, finite fuel
and binding energies are included in H; it exchanges heat with one bath and
receives no work after preparation. Each enabled reaction pair obeys

    log[k(x,y)/k(y,x)] = -(H(y)-H(x))/theta.

On a connected communicating class its equilibrium reference is
pi(x) proportional to exp[-H(x)/theta]. If interventions split the state graph,
use the reference restricted to each class and its conserved quantities.
Paired barrier changes alter relaxation kinetics without introducing a drive.

For the probability distribution p(t), define relative entropy

    D(p||pi) = sum_x p_x log(p_x/pi_x).

Set J_xy = p_x k_xy - p_y k_yx. The master equation and detailed balance give

    -dD/dt = (1/2) sum_(x,y) J_xy log[(p_x k_xy)/(p_y k_yx)] = sigma(t) >= 0.

The inequality follows term by term from (a-b) log(a/b) >= 0. With dimensionless
Shannon entropy S(p) = -sum p log p and F[p] = E_p[H] - theta S(p),

    F[p]-F[pi] = theta D(p||pi),
    integral_0^T sigma(t) dt = D(p0||pi)-D(pT||pi) <= D(p0||pi).

The initial preparation therefore supplies a finite availability budget.
A finite irreducible model of this type cannot sustain a stationary positive
entropy-production rate forever. Initial deterministic preparations are
permitted: since finite-state pi is positive, their D is finite; derivatives
at time zero can be interpreted by a positive-distribution limit.

## Consequences for this research

The turnover comparison can establish that paid renewal prolongs activity
over a specified finite horizon, including formation costs. It cannot derive
indefinite self-maintenance or open-ended recursive biodiversity. Faster L
turnover does not change that limitation, even if the comparison succeeds.

Integrated R occupancy alone also cannot be identified with accumulated
negentropy. Thermal equilibrium can have nonzero R occupancy and hence an
unbounded occupancy integral as the observation horizon grows, despite zero
stationary dissipation. Continuing nonequilibrium capacities need their own
physical and causal tests.

The standing-wave ambition eventually requires explicit free-energy flow,
including fuel supply and waste removal or a finite physical apparatus that
provides them. A chemostatted extension is a declared open-system idealization
whose chemical work must be counted. Neither a replenishment reward nor an
unaccounted reset is an acceptable substitute. That future extension follows
the present finite-budget test; it is not secretly present in its results.

Finally, this free-energy decay identity is not a law of biological aging.
Function can improve while available fuel is spent on repair. Deterioration,
repair and regulatory information must be represented by physical states and
interventions; neither elapsed time nor accumulated expenditure supplies them
by definition.
