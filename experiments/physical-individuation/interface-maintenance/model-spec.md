# Local production and renewal of protective material

This experiment tests a supplied, finite-fuel local mechanism. Templates produce
material that can form protective contacts, and those contacts modify the rates
of template loss and transport. It does not identify an enclosing membrane by
definition, prescribe groups, or reward a desired morphology. Formation from
dispersed founders and maintenance of an already formed arrangement are separate
assays of the same generator.

## State and conserved inventories

The lattice is a periodic square of `side >= 3`. Each site has exactly one
material unit, `D`, `R0`, `R1`, or `L`, and one independent fuel slot, `F` or `W`.
`D` is precursor material, not an empty site. Consequently material units and
fuel-plus-waste units are separately conserved. Both layers move by nearest
neighbor exchanges. The initial state has dispersed template founders, no `L`,
precursor at all other sites, and a declared finite fuel stock. No feed enters.

`R0` and `R1` have equal energies and, by default, identical rates and
vulnerabilities. Their default distinction is not evidence of functional
biodiversity. A supplied copying kernel can change their labels. Optional
allele-specific passive multipliers permit later response experiments, but do
not invent an environmental response map.

Instance tags move with the corresponding material species and are replaced on
chemical creation. They identify molecular instances and copying provenance,
not conserved atoms, subject identity, functional memory or thermodynamic
microstates. They never enter a rate. Passive degradation and replacement are
turnover; this model has no accumulated damage or age-dependent deterioration
and therefore supplies no model of aging.

## One Hamiltonian and local reversible channels

Each undirected nearest neighbor edge is counted once. Let `R` mean either
template label. In energy units shared with `theta = k_B T`,

```text
H = epsilon_R N_R + epsilon_L N_L + E_F N_F
    - J_LL number(LL edges)
    - J_RL number(RL edges)
    - J_RR number(RR edges).
```

`J_RR` defaults to zero. `J_LL` and `J_RL` are supplied interactions that can
favor contact assembly. Neither their presence nor an observed cluster proves
an enclosing surface, selectively permeable compartment or higher individual.
Rates have units inverse time; concentrations, pair-count rewards and cluster
labels are not used. The full state explicitly records positions, so no
additional combinatorial entropy term is inserted into `H`. Entropy belongs to
the probability distribution over those positional states.

Every allowed edge uses an edge-symmetric mobility `a` and the Metropolis rule

```text
q(s -> s') = a(s,s') exp[-max(H(s')-H(s), 0)/theta].
```

Thus `q_forward/q_reverse = exp[-delta_H/theta]` on every retained edge. Zero
mobilities remove both directions. The supplied local channels are:

1. Passive `R_a <-> D` and `L <-> D` at one material site.
2. For each adjacent template helper `R_b`,
   `D_i + F_i + R_b,j <-> R_a,i + W_i + R_b,j`, with mobility
   `copy_rate * mutation[b,a]`.
3. For each adjacent template helper,
   `D_i + F_i + R_b,j <-> L_i + W_i + R_b,j`, with mobility
   `production_rate`.
4. Exchange any two different neighboring material species, or exchange
   neighboring `F/W` states, using the corresponding symmetric mobility.

The helper and target are different neighboring sites. A reverse copying edge
has the same product-dependent copying coefficient and the same unchanged
helper as its forward edge. Different helpers are separate parallel catalytic
channels; their event indices remain separate in the ledger. Nonequilibrium
initial fuel supplies the transient chemical driving. In a finite closed
system, these dynamics do not support perpetual stationary fuel consumption.

## Protection, retention and factorial controls

The passive `R_a <-> D` mobility is multiplied in both directions by
`exp(-protection_strength * n_L)`, where `n_L` counts neighboring `L`. It changes
a barrier, not the site energy or equilibrium preference.

For a material exchange involving at least one template, mobility is also
multiplied by `exp(-transport_protection_strength * n_L_union)`. The local union
contains both endpoints and both sets of nearest neighbors. Its number of `L`
units and the presence of a template endpoint are unchanged by reversing the
swap. This factor therefore gives symmetric kinetic retention. It does not
classify a move as entering or leaving an inside region.

`produce_L=False` removes both directions of fueled `L` production.
`copy_R=False` removes both directions of fueled template copying.
`protect_R=False` removes **both named passive and transport barrier effects**.
All these controls retain the same Hamiltonian, possible species and initial
state when initialized with the same seed and parameters. Thermal `L <-> D`
formation remains present in the production-off control. Energetic `R-L`
attraction also remains when protection is off. Thus this factorial comparison
isolates the combined additional barrier effects, not all physical effects of
`L`, and cannot attribute a result specifically to transport retention. A
strength can separately be set to zero for an explicitly identified follow-up.

For fixed environment multipliers, optional allele vulnerabilities multiply
both directions of the corresponding passive channel and retain the same
Gibbs invariant distribution. A prescribed changing environment is not a
modeled gate or apparatus; its physical realization and control costs would
require an extension. No such changing environment is needed for the first
interface experiment. Removing channels can make sectors disconnected; the
common Gibbs distribution remains invariant but need not be the unique
reachable equilibrium.

## Event accounting and physical scope

The full system is closed to matter and does no external mechanical work. It
exchanges heat with a bath at fixed `theta`. Each jump records `delta_H`, heat
to the bath `-delta_H`, external work zero, and fuel before and after. Therefore

```text
H_final - H_initial + cumulative_heat_to_bath = 0.
```

Fuel standard-energy change is `E_F * delta_N_F`. Material energy change is
`delta_H - E_F * delta_N_F`, including interactions. Their sum is the full
energy change. The standard energy per packet is not the entire finite-stock
chemical free energy: positional mixing entropy and correlations matter.
Counting `E_F` per net consumed packet as extracted work or full chemical work
would be incorrect. No nonequilibrium available-free-energy calculation is
claimed by this prototype.

Any maintenance injury or harvested initial arrangement must be prepared and
reported separately, retaining the changed species, inventory and energy. A
conditioned initial state is not an autonomous wound event and its selection
or preparation is not a costless controller. Harmless tags may quantify
turnover; organizational continuity and functional information require
additional causal assays.

## Stable executable API

```python
model = Model(Params(...), Flags(produce_L=True, protect_R=True, copy_R=True))
state = model.initial(seed)
events, propensities = model.rates(state)  # fixed catalogue; zero rates retained
state.time = event_time
record = model.apply(state, event_index)   # mutates physical state and tag metadata
reverse_event = model.reverse(event_index)
reverse_index = model.reverse_index[event_index]
delta_H = model.delta_energy(state_before_event, event_index)
H = model.energy(state)
```

`Event` contains `kind,i,j,before,after`; kinds are passive, copying, protective
production, material exchange and fuel exchange. `State` contains flat material
and fuel arrays, instance tags, the next tag, time, cumulative heat, external
work and event count. `copy_state`, `state_to_dict` and `state_from_dict` support
independent replay. `apply` can apply a catalogue edge for audits even if a flag
sets its mobility to zero; simulation selects positive-rate edges only.

The convenience `simulate(params,seed,horizon,flags,sample_dt=...)` uses direct
SSA and returns complete parameters/flags, initial and final states, observation
samples and the full event ledger. It does not save files or use observations
to change rates. Each record includes the event index/sites, energy, heat/work,
fuel counts, new/removed instance tags, a copying parent tag and a catalyst tag.
The copying parent applies only to template copying; `L` production records a
catalyst without implying inherited `L` identity. Reverse events restore the
physical species/fuel state, not arbitrary observer-assigned instance IDs.

Use the bundled Python with NumPy rather than the system Python:

```text
[local path omitted]
```

The module's command-line entry performs a short software smoke without saving
results. A successful experiment still requires independent edge-ratio tests,
full-energy/event replay, matched initial-state checks, replicate trajectories
and readout controls. Longer persistence alone does not establish lower cost at
matched continuing function; lower cost, accumulated damage, functional memory,
collective heredity and higher individuation remain distinct claims.
