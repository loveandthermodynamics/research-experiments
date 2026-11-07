# Differentiated response and protective renewal

This extension supplies two physically different response phenotypes within the
existing finite chemical alphabet. It imports the frozen driven spatial model
and leaves its Hamiltonian, reservoir channels, copying kernel, local protective
material chemistry, transport and neutral instance IDs unchanged. There is no
group, diversity, lineage or morphology term in a rate.

## Symmetric autonomous environment

The physical gate has states g=-1 and g=+1, both with well energy zero. It flips
in either direction at rate kappa. Its stationary distribution is one half in
each state. By default a separate seeded random stream draws that stationary
initial state independently of the material initialization. A gate is never
reset at a continuation boundary; a supplied joint initial state includes its
actual realized gate. A new exponential waiting time after such a boundary is
Markov-valid, but replacing its physical value by a preferred phase is not.

For epsilon in [0,1), R0's passive mobility multiplier is 1+epsilon*g.
In the opposed response map, R1's is 1-epsilon*g. In the aligned control, R1 has
the same multiplier as R0. The multiplier applies to BOTH directions of the
same passive R_a<->D channel. It changes neither state energies nor forward and
reverse copying rates. Each allele has marginal stationary mean multiplier1
in either response map, while the temporal correlation between vulnerabilities
changes. Actual finite-trajectory exposure need not average to1 and must be
reported rather than assumed. R0 follows exactly the same environmental map
and gate history in the paired controls.

Gate jumps have deltaH=work=heat=0. Every retained internal edge still satisfies
log(q_forward/q_reverse)=-deltaH/theta. Reservoir edges retain
log(q_forward/q_reverse)=(chemical_work-deltaH)/theta and the full balance
H_final-H_initial=cumulative_chemical_work-cumulative_heat_to_bath. In the mu=0
limit the joint generator obeys full-H detailed balance with a uniform gate
factor. At driven mu, local detailed balance remains while reservoir cycles
drive the material. The gate is a supplied coarse kinetic apparatus; its
physical realization is not derived, and zero well energies do not justify
arbitrary gate resets or ignoring material/gate correlations.

## Copying and source attribution

Every study cell retains the copying kernel [[.95,.05],[.05,.95]]. All founders
are R0 and no R1 is supplied. R1 can first appear through a fuel-paid R0-helper
copying error or rare passive thermal formation. The latter remains present
and its multiplier changes with the gate too; it cannot be silently excluded.
The first realized appearance and all later birth sources must be separated.

Copy events record helper_species as well as parent_instance and newly created
instance IDs. Passive thermal births have no copying parent. These observations
allow an externally analyzed copying subtree to include backmutation into R0;
R1 class occupancy is not the same as descent from one founding R1 instance.
Tags never enter propensities, and labels do not follow conserved atoms after
chemical degradation. The response functions and finite alphabet are supplied;
the model does not invent an unlimited grammar or derive human values.

## Optional fixed-site intervention

The natural model has production_site_mobility=1 at all sites. An optional
immutable nonnegative array multiplies both L-production directions by the
entry at event j, the helper/catalyst site. The same fixed factor appears on
the reverse edge, leaving H and local detailed balance unchanged. It does not
alter copying or passive L chemistry. The array is copied at construction,
stored read-only and serialized with the run. It cannot be changed by lineage,
group success or observations during a trajectory. Building the external
kinetic intervention is outside the model and uncosted; it is not a naturally
generated interface or a controller for genealogy. Root's later ownership
assay is separate from the natural response experiment.

## Proposed bounded design, pending physics/protocol review

The four settings are epsilon=.5/.8 crossed with kappa=.025/.1. All remaining
parameters match the selected driven study: eight dispersed R0, zero L,58F,
theta=.7,epsilon_R=epsilon_L=8,E_F=12,J_LL=J_RL=1,J_RR=0,copy_rate=.028,
production=.04,passive_R=.12,passive_L=.08,protection=1,transport_protection=.25,
local transport and reservoir pF=.9/exchange_rate1. Production-site mobility is
all ones. The gate starts stationary with separate initial/clock/material RNGs.

At each setting, eight pilot seeds run the 2x2 of opposed/aligned responses and
L production enabled/disabled from the same pure-R0 initial distribution to
time100. Protection and template copying remain enabled in every cell. The
production-off control removes both fuel-paid L directions; passive thermal L
formation and energetic L interactions remain. Every seed and failed state is
retained. These are end-to-end formation-plus-maintenance comparisons, not an
isolated intervention on an already established assembly.

The proposed primary contrasts are the opposed-minus-aligned R-time effect
with production and its response-by-production interaction. A deterministic
pilot rule selects one setting by their minimum standardized paired effect,
with explicit tie order, before one fresh confirmation. No additional scan or
rescue is part of this design. Root owns the final protocol, focal-lineage
establishment criterion, causal ownership interventions and study runner.

Full initial resources/H, reservoir work, heat, event histories, actual gate
occupation, thermal versus paid appearances, molecular turnover and spatial
readouts must remain available. Bath heat is not total entropy production.
More maintained R under the supplied drive would be a conditional functional
advantage, not by itself an individual, an enclosing membrane, protected
participant agency, aging, functional memory or moral legitimacy.

## API

DifferentiatedModel(params,flags,drive,response,response_mode,transport_mode,
production_site_mobility=None) exposes initial, copy_state, energy, rates,
material_rates, apply, reverse, reverse_index and state_to_dict/from_dict.
JointState contains base_state, gate and gate_event_count; time proxies the
base state's time. Serialized states are nested with those exact keys.
ResponseParams(asymmetry,flip_rate,initial_gate=None) uses physical -1/+1 values.
simulate_differentiated returns full states/parameters/intervention array,
samples with material/tags, and the complete combined ledger. A software smoke
is permitted before review; no pilot or confirmatory ensemble is run here.
