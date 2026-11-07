# Exact balances before spatial claims

This note concerns the proposed local productive-copying generator. It does not transfer the ideal-vessel fixed point to a lattice, predict a droplet radius, or establish biodiversity. The separate model specification must supply the actual elementary rates and all resources.

Let B be the global number of sticky material carriers, N the total number of templates and M the conserved number of A/B carriers. Let C_X and C_Y count signed forward-minus-reverse copies that respectively transform A→B and B→A. Offspring state can differ from parent state through a fixed paired copying-error kernel. Let C_0 count every other signed A→B transformation. I counts non-copy template births, L template losses. All counters refer to actual elementary events, not an observer's assignment to a group.

Pathwise,

    B_T-B_0 = C_X-C_Y+C_0
    N_T-N_0 = C_X+C_Y+I-L.

Material/template transport changes neither identity. If each net copy receives chemical work mu, then

    W_chemical = mu*(C_X+C_Y)
               = mu*(N_T-N_0-I+L).

At stationarity, using lower-case letters for mean signed currents,

    j_X = (l-i-j_0)/2,   j_Y = (l-i+j_0)/2.

Thus positive mean currents in both transformations require l-i>|j_0|. Their existence is a quantitative renewal condition, not a diversity score. Neither spatial localization nor separate taxa follows from these global equalities.

## Explicit birth and loss bounds

For a template layer with D precursor or template X/Y at each site, use standard energies E_X,E_Y and on-site binding energy -g_r*1(template=r)*1(material=B). A passive D↔r pair has mobility delta_r and rate delta_r*min(1,exp(-DeltaH/theta)). For E_r>=max(g_r,0), every passive template loss is downhill and has rate delta_r. With equal delta, total loss intensity is exactly delta*N.

For the analytical example E_X=E_Y=6, theta=1, g_X=-2.5, g_Y=2.5, define

    c0 = 2*exp(-6)
    cB = exp(-8.5)+exp(-3.5).

If there are K lattice sites and D_B precursor slots co-located with B,

    i(state) = delta*[c0*(K-N) + (cB-c0)*D_B].

This is an exact propensity sum for that passive channel law. Since D_B<=M and cB>=c0,

    E[i] <= delta*[c0*(K-E[N])+(cB-c0)*M].

If a paired passive A↔B background of mobility kappa is included, each direction's rate is at most kappa per eligible site. Therefore |j0(state)|<=kappa*M. Under either productive-family blockade, the stationary mean satisfies

    E[N] <= [kappa*M/delta+c0*K+(cB-c0)*M]/(1+c0).

For other binding values satisfying the stated downhill-loss condition, define c0 and cB from their corresponding Boltzmann birth factors and replace cB-c0 in this upper bound by its positive part. Outside that condition, use the actual loss lower bound and the saturated birth law instead. The bound is on stationary expectation, not each instantaneous count or the persistence of a particular ancestry. For K=256,M=77,delta=.02,kappa=.0002 it is about four templates. This is an analytical illustration, not a selected outcome regime or an established full-model population.

## Integrated maintenance cost

The same exact counting identity gives

    E[W_T]/T = mu*[(E[N_T]-N_0)/T
                         +delta/T*integral_0^T E[N_t]dt-E[I_T]/T]

when passive loss is delta*N. Stored initial templates and all non-copy appearances remain in the expression. At stationarity, power is mu*(delta*E[N]-E[i]). A proposed benefit that increases maintained population can therefore increase fuel use; this is not a reduction in dissipation at matched useful function. Formation, material supply, repair, transport and environmental apparatus must be included before comparing complete life-cycle costs. No damage-history or informational aging variable has yet been specified.

## Complementary-cycle theorem when kappa=0

With only X-family copying, work on every transition is mu*DeltaB; swaps and passive template births/losses have DeltaB=0. Local detailed balance then becomes equilibrium detailed balance for H_eff=H-mu*B. With only Y-family copying, H_eff=H+mu*B. Thus a one-family model can have thermal occupancy and a B-rich structure without a stationary net copying-work current. At kappa>0 the passive A/B channel is an explicitly bounded bypass; do not transfer the exact equilibrium result to that different generator.

With both productive families, a local four-step cycle is possible. Prepare a target with A/D and two distinct neighboring parents X and Y. X-assisted copying creates B plus a template at the target; passive template loss restores D there. Y-assisted copying returns the target to A plus a template; passive loss restores A/D. The parents and all other physical occupancies are unchanged, H returns to its original value, and two fuel units have been consumed. The forward/reverse cycle-rate ratio is exp(2*mu/theta), using the same paired mutation factors in both directions.

For a finite irreducible chain with positive reverse rates and stationary distribution pi, define elementary probability flux F=pi(s)*q(s,s') and its paired reverse F_rev=pi(s')*q_rev(s',s). Stationary entropy production is a sum over unordered elementary pairs of (F-F_rev)*log(F/F_rev), each nonnegative; equivalently use half the directed sum. Raw rates alone are insufficient. This sum vanishes only if each paired probability flux balances. That would force every cycle ratio to one, contradicting this cycle for mu≠0. At mu>0, stationarity and kappa=0 therefore imply j_X=j_Y>0 and power=2*mu*j_X. This proves a genuinely dissipative cycle in the supplied grammar, not a sufficiently large maintained population, a phase unit, ancestry, aging or separate taxa.

Irreducibility must be checked for the actual finite state space: connected local transport can rearrange material; positive thermal D↔X/Y channels can prepare parents and targets; productive copy followed by loss changes A/B while restoring template occupancy. Nontrivial material inventory and at least two distinct neighbors are premises. Neutral instance IDs are observer bookkeeping and are not part of this finite physical state or its thermodynamic return condition.

## Local integrity remains another test

For any fixed spatial region U, the corresponding B balance includes net material flux through its boundary, as well as the local X/Y/background reactions. A region can be supported by remote production. None of the global equations identifies the owner of a flux, proves closure or defines an individual. A moving or history-selected region needs its own transport/accounting convention. Molecular partitioning, dense formation, causal renewal dependence, continuity through turnover and adverse generalist/bypass comparisons must be measured in the same model before making the stronger claim.
