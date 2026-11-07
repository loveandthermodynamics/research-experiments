# Driven local protection and its transport control

The experiment asks whether ongoing local production of protective material
supports template maintenance under an explicitly supplied chemical drive.
It preserves the frozen sibling model's species, Hamiltonian, local chemistry,
copying kernel, contacts and neutral instance tags. No cluster, inside/outside,
identity, shape or fitness measurement enters a rate.

Each of64 lattice sites holds one material unit D/R0/R1/L and an independent
F/W slot. Material inventory and fuel-plus-waste site inventory remain fixed.
The initial eight dispersed R0,58F, remaining precursor and fixed apparatus are
supplied resources. Their initial stored energy H is reported; H is not called
initial exergy or counted as reservoir work. R0/R1 have the same energies and
default responses; the finite label alphabet does not by itself demonstrate
functional biodiversity.

The full Hamiltonian is unchanged:

```text
H = epsilon_R N_R + epsilon_L N_L + E_F N_F
    - J_LL LL_edges - J_RL RL_edges - J_RR RR_edges.
```

All original internal channels retain Metropolis local detailed balance with
this H. The added coupled reservoir exchange at site i replaces W by F or vice
versa. With exchange scale rho and reservoir equilibrium fuel occupancy p,

```text
q(W_i -> F_i) = rho p
q(F_i -> W_i) = rho (1-p)
mu = E_F + theta log[p/(1-p)]
```

The reservoir supplies the chemical-potential difference mu per net fuel
packet. For an up-jump deltaH=E_F, W_chemical=mu, Q_bath=mu-E_F. All three signs
reverse for a down-jump. Thus every channel obeys

```text
log(q_forward/q_reverse) = (W_chemical - deltaH)/theta
H_final - H_initial = cumulative_W_chemical - cumulative_Q_bath.
```

For internal moves W=0 and Q=-deltaH. Every exchange is recorded, and the state
field `external_work` accumulates chemical reservoir work. With mu=0 the entire
enabled generator obeys full-H detailed balance. Positive drive breaks global
equilibrium detailed balance around fuel-consuming cycles while preserving
local detailed balance. Cumulative bath heat can have either sign and is not
the total stochastic entropy production; state entropy and reservoir terms
cannot be omitted. Uniform reservoir access is a supplied physical drive, not a
group reward, evolved resource-distribution system or derived apparatus.

## Locality control

Both modes retain the same reservoir and Hamiltonian. The local mode retains
nearest neighbor material transport. The mixed_L mode disables only base swaps
whose endpoints include L, replacing them by all unordered pair exchanges
including L with mobility4M/(N-1), where M is the original nearest-neighbor
mobility. All other original swaps remain local. This is a different supplied
passive transport graph, not instantaneous randomization or an externally
powered shuffling operation.

Each L has proposal-clock budget4M in either graph when null L-L exchanges are
included. Accepted flux and R displacement need not match: an R-L exchange also
moves R, and interaction energies and protection modify acceptance. Those are
part of the changed transport mechanism. The control does not claim to match a
literal well-mixed continuum or construction costs of two transport fixtures.

Every extra swap uses actual full-H energy differences. The double-counted
shared-edge correction is applied only when its endpoints are nearest
neighbors; a nonadjacent pair has no such correction. If protection is enabled,
R-containing swaps use the endpoint-inclusive neighbor union to count L. This
count is unchanged on reverse exchange, preserving symmetric barriers. All
mixed swaps conserve material and tags and have zero reservoir work.

## Declared design

The pilot keeps theta=.7,epsilon_R=epsilon_L=8,E_F=12,J_LL=J_RL=1,J_RR=0,
production=.04,passive_R=.12,passive_L=.08,protection=1,transport_protection=.25,
reservoir p=.9,rho=1, and the supplied copying-error kernel. Only copy rate
varies over .020,.025,.028. In every case4*copy_rate<passive_R, so unprotected
copying has a subcritical maximal per-template birth rate before thermal
immigration. Thermal reverse formation remains present and is counted; zero
templates is not absorbing.

Eight pilot seeds at each level undergo local full P1B1 preparation from zero L
for20 time units. Every resulting state, including failures, enters all four
production/protection flag arms in each of the two transport modes for80 more
units. Only bookkeeping clocks/counters are reset. No material, fuel, tag or
next-tag value is reset. The graph changes after local preparation, so this
tests postformation transport rather than formation under mixed transport.

Three paired primary R-time contrasts are continuing production with protection,
the local production-protection interaction, and local minus mixed_L full-arm
maintenance. The one copy rate maximizing the minimum of their pilot mean/SE
scores is selected, ties going to the lower rate. All pilot levels and signs
remain saved. A separate64-seed confirmation supports the declared joint claim
only if all three approximate95percent normal intervals are positive. No
additional parameter search is authorized by this protocol.

Preparation histories, initial resources/H, chemical work, heat and full event
ledgers are saved. Branch totals include preparation work and heat. Secondary
observations include L renewal, separately counted thermal births, copying
after founder-instance turnover, final contact components and sampled contact
enrichment/contact-time integrals. Contact excess subtracts the fixed-count
uniform-placement expectation4*N_R*N_L/(N-1). Lowest nonzero R Fourier power is
normalized by its fixed-count random expectation N_R*(N-N_R)/(N-1); its null
mean is1 and it is undefined at zero or full R occupancy. These observations
never alter rates, identify a person, or prove individuality.
Instance turnover is not aging or functional information. A locality advantage,
if found, would not by itself prove an enclosing membrane, independent
participants, collective heredity or a new individual.
