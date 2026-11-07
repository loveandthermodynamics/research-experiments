# Local active-phase generator: prospective model specification

This is a new finite spatial Markov model, not a fit to the earlier lattice results or a reproduction of a published solver. Its immediate purpose is to test whether dispersed material can form and turn over in localized structures under a thermodynamically consistent driven chemistry. Finite size and arrested coarsening are separate hypotheses, not features guaranteed by a cluster rule.

The design follows the mechanism in [Kirschbaum and Zwicker (2021)](https://arxiv.org/html/2103.02921v1): a soluble/segregating interconversion with a second driven channel, and enzyme partitioning that can create different net reactive behavior in compositionally distinct regions. In their model the fuel-driven channel destroys segregating material. The present discrete transition law differs from their continuum transition-state kinetics. Agreement of mechanisms does not establish agreement of phase diagrams or selected radii.

## Microscopic state and supplied resources

A periodic square lattice of side `side>=4` has one material slot at each site: solvent S=0, soluble A=1, or sticky B=2. A separate dilute-enzyme layer has occupancy E_i=0/1. Enzyme molecules are conserved and mobile; the second-layer representation neglects their excluded volume in the material layer. That is a supplied coarse-graining assumption.

The exact material inventory N_A+N_B, solvent inventory and enzyme inventory are conserved. Initially a fixed number of A molecules and enzymes are distributed uniformly without replacement, independently of one another. Initial B is exactly zero. No droplet, interface, inside/outside field or favored site is supplied. The molecules, their affinities, reaction channels, lattice and heat/chemical reservoirs are supplied. This is unseeded phase formation from supplied chemistry, not the origin of chemistry.

Material molecule tags, if present in the saved state, identify the conserved A/B carriers. They move in swaps and persist through A/B conversion; solvent has tag0. They never affect rates and do not represent heritable function, a template lineage or an individual.

## One Hamiltonian

For b_i=1 when material_i=B and zero otherwise,

`H = epsilon_B*sum_i b_i - J_BB*sum_<ij> b_i*b_j - J_EB*sum_i E_i*b_i`.

Each undirected nearest-neighbor edge is counted once. Solvent and A have zero reference energy. Positive J_BB is a short-range attraction; positive J_EB supplies enzyme affinity for B. Both the material and enzyme motion include the resulting energy changes. No prescribed partition coefficient or droplet position enters the generator. Mixing entropy is implicit in the multiplicity of lattice microstates; it is not added a second time to H.

Temperature theta=k_B*T and epsilon_B, J_BB, J_EB, mu have units of energy. Mobility/rate prefactors have units inverse time; lattice spacing is the length unit. Physical calibration is not asserted.

## Elementary reversible channels

1. Material exchange across each undirected nearest-neighbor edge. Distinct S/A/B occupancies are swapped. Conserved material tags move with their carriers.
2. Enzyme exchange across each such edge when E_i differs from E_j. The material is unchanged and the enzyme-B binding contribution to H changes normally.
3. Passive A_i <-> B_i at every site, with no external work.
4. At sites containing an enzyme only, `B_i + F <-> A_i + W`. Positive mu=mu_F-mu_W supplies chemical work +mu for B->A and -mu for A->B. The enzyme is not consumed and is present in both endpoints. Each fuel-to-waste conversion is counted; the reverse synthesizes fuel.

F/W are explicitly declared ideal homogeneous chemostats, not finite locally stored fuel. Their concentrations and chemical potential difference are maintained externally. The model records their chemical work but does not build the reservoir apparatus or enzyme molecules. There is no unpriced periodic reset or replenishment of material.

For every elementary edge x->y, with energy change deltaH=H(y)-H(x), work w and a reverse-symmetric mobility m,

`q(x,y) = m*min(1, exp((w-deltaH)/theta))`.

The passive mobility is `passive_rate`; the catalytic mobility is `catalytic_rate*E_i`; swap mobilities are the respective diffusion prefactors. Thus

`log(q(x,y)/q(y,x)) = (w-deltaH)/theta`.

Every finite enabled channel has a finite reverse. At mu=0, all channels obey detailed balance with the same Gibbs reference exp(-H/theta), conditional on the conserved inventories. At nonzero mu, the cycle passive A->B followed by catalytic B->A has affinity mu/theta. These assertions concern elementary reaction channels; two channels connecting the same material endpoints retain separate identities and ledgers.

## Physical ledger and executable API

Each event records time, event index, kind, endpoints, pre/post occupancies, delta_H, chemical_work, heat_to_bath, and signed fuel_to_waste. The heat sign is positive toward the thermal bath:

`Q_event = W_event - delta_H`, and `H_final-H_initial = W_chemical-Q_bath`.

The cumulative counters and full initial/final microscopic states are saved for replay. Initial stored H and inventories are supplied resources; they are not preparation work or distribution-level exergy. Finite-trajectory bath heat need not be positive and is not alone total irreversible entropy production.

The prototype will expose `Params`, `State`, immutable `Event`, and `Model(Params)`, with `initial(seed)`, `rates(state)`, `apply(state,event)`, `energy(state)`, `delta_energy(state,event)`, `reverse_index`, `state_to_dict`, and `state_from_dict`. `simulate(params,seed,horizon,...)` returns full initial/final states, a complete event ledger and sampled microscopic states. Connected components, structure factors and apparent radii are readouts only and are never available to the natural rates.

## What a later bounded experiment must distinguish

Formation means statistically supported segregation/localization from the zero-B randomized state, compared with appropriate fixed-count spatial nulls. Multiple small clusters during nucleation alone do not establish selected size. Intrinsic finite scale requires size/structure statistics that stop growing over a declared late window and remain below the box scale as box size increases, at matched densities; transient slow coarsening is an alternative. Turnover requires forward and reverse channel counts and material residence/readout histories in a continuing pattern. It does not establish aging, organizational memory, reproduction or heredity.

A mu=0 control tests removal of chemical drive while preserving H. A J_EB=0 comparison removes binding-mediated partitioning but also changes H, so it is not a pure kinetic ablation. Catalytic-channel removal is paired in both directions and keeps H, while changing the accessible kinetic cycle. No control is interpreted as proving individuality from morphology alone.

No ensemble or broad parameter scan is authorized by this specification. Software checks must first verify local detailed balance, all energy deltas, conserved inventories, zero-B initialization, the mu=0 limit and complete first-law replay. Parameter choice and morphology decision rules require a separately frozen bounded protocol. Copying, mutation and productive differentiated roles are not yet in this model.
