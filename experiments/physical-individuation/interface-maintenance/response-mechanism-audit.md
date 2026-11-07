# Why response diversity has not yet funded collective rebuilding

This audit reads the frozen ensemble summaries and event ledgers. It creates no new trajectories and changes no prior model or result file. `audit_response_mechanism.py` reproduces the ledger diagnostics below using only JSON/gzip reads; it prints JSON and writes nothing.

## What the results establish

The baseline 64-seed factorial supports a causal production/protection loop: integrated R is 677.77 with both mechanisms, 420.47 without L production, and 202.60 when L production remains but its protective barriers are removed. That result does not automatically establish an advantage of differentiated response types.

In the response experiment, integrated R is 387.80 with differentiated exposure and copying variation, versus 395.47 without copying variation. The response-by-variation interaction is 20.90 with Monte Carlo SE 17.13 and approximate interval [-12.67,54.47]. Preserve this inconclusive primary result. It establishes neither a positive response-diversity advantage nor a demonstrated harmful effect. First-episode lifetime rises from 52.70 to 61.32 in that comparison while integrated activity does not rise: longer thin persistence is distinct from rebuilding more active capacity.

The varied response arms use mutation probability **.2**, overriding the .05 mutation matrix in `base_params`; the actual per-trajectory matrix is [[.8,.2],[.2,.8]].

## A neutral-label limitation, not an assumed mutation tax

R0 and R1 have identical energy, binding, L production, and copying capabilities. Each can maintain R and produce L without the other. Their functional difference is which environmental gate accelerates passive turnover. There is no obligate complementary chemistry between the two types.

Under common exposure, mutation does not change aggregate forward copy propensity: row-stochastic copying coefficients sum to `copy_rate` for every eligible R–D–fuel edge. With reverse copying neglected, the aggregate R/D/L/fuel process is lumpable and independent of the mutation matrix. R0/R1 swaps are invisible internal transitions of that aggregate. Paired reverse copying is the small exception, because its coefficient depends on parent/product labels; no reverse-copy event occurred in the three ledger groups reconstructed here.

Consequently the common-varied mean 365.67 versus common-unvaried 394.23 is not evidence of a built-in energetic mutation tax. A variant copy consumes fuel, but the corresponding clonal copy consumes the same packet. The observed common-arm difference can be Monte Carlo variability. The missing mechanism must be diagnosed from differential exposure and regeneration, rather than supplied by that apparent penalty.

## The observed healthy insurer is not self-renewing fast enough

First paid R1 appears at mean time 4.86 with about 53.97 of 64 fuel packets remaining. Novelty is therefore generally early enough to matter. At the first gate-0 bad episode, however, both types remain present in only 35 of the 58 runs that still have any R; the mean surviving R1 count is 1.38. Rare founding and stochastic loss remain real bottlenecks.

For every interval with gate 0 bad and gate 1 healthy, reconstruct R1's local contexts. Let c be its time-weighted number of adjacent D sites carrying fuel, and a its time-weighted protective factor exp(-n_L). The frozen ledgers give

\[
c=.78182650,\qquad a=.39093611.
\]

The rates at the original parameters are then

\[
b_{11}=k_{copy}\mu_{11}c
=.06\times.8\times.78182650=.03752767,
\]
\[
d_1=\delta_R a=.12\times.39093611=.04691233.
\]

Thus the healthy R1 lineage's observed self-copy intensity is below its death intensity. Copying into R0 and incoming R0-to-R1 copies still occur; this calculation isolates the insurer's ability to maintain itself after the exposed population becomes unavailable. These conditional averages are **not** an exact invasion eigenvalue and are not a counterfactual forecast after parameters change.

Across all 64 differentiated-varied trajectories, during these exposure intervals R1 catalyzes only 25 R copies, including six R0 births, versus 39 L births. R0 is reintroduced from an R1 catalyst after its previous loss ten times over whole trajectories. The return path exists, but its realized throughput is small. Mean fuel remaining at first total-R extinction is 18.03 packets. This is not simply a globally exhausted reservoir: fuel-bearing precursor access and allocation fail before all stock is spent.

## Exact allocation constraint in the current grammar

Every eligible forward R or L production edge is energetically downhill at these parameters: E_F=12 exceeds onsite energies 8, and binding only lowers product energy. Therefore the two forward propensities are exactly

\[
J_R^+=k_{copy}C,\qquad J_L^+=k_L C,
\]

where C counts directed R–D–fuel contacts. Their forward birth split is .06:.08: about 57.14 percent of paid births make L. Ignoring the unobserved reverse events, expected fuel expenditure per R birth is 1+k_L/k_copy=7/3 packets, because L production competes for the same precursor/fuel opportunities.

L turnover is also real: its passive death rate is .04 per L because all possible L-removal energy changes are downhill. A repair episode must simultaneously satisfy local precursor access, template replacement, and roughly J_L^+ >= .04 L to maintain its protective material. Saving a variant without funding those return flows need not rebuild the organization.

## The next bounded test

Test **copy_rate=.12 versus .06**, retaining the same chemistry, energies, production rate, mutation matrix, finite inventories, initial monomorphic preparation, and autonomous gate regime. Use all four response-by-variation arms and the newly frozen disjoint seeds 148001–148064. This changes a paired physical kinetic coefficient; it adds no reaction, local rule, diversity reward, or free fuel.

The readout-derived self-renewal boundary under the OLD contexts is

\[
k_{copy}>\frac{\delta_R a}{\mu_{11}c}=.07500439.
\]

At .12, the old-context diagnostic would become b11-d1=.02814301 rather than -.00938466. The expected fuel burden per R birth also falls from 7/3 to 5/3 as the R:L birth split changes. These provide a mechanistic direction before seeing new results. New kinetics will change contacts, L abundance, crowding, and local fuel; all contexts and the inequality must be reconstructed anew. Merely crossing .075 in a parameter file does not verify the mechanism.

This dose does not simply make exposed types invulnerable. For a template with n_L neighbors, its total offspring rate from its own catalyst is bounded above by .12(4-n_L), while severe-exposure death is 4.8 exp(-n_L). For n_L=0,...,4 the former remains below the latter. The exposed class retains a negative own-parent growth bound, leaving a kinetic opportunity for healthy-class insurance. Common-response populations can nevertheless survive finite outages probabilistically; their failure is not assumed.

Predeclare the mechanistic checks alongside the original primary endpoint: variant establishment before exposure; healthy-survivor self-copy minus death intensity; paid template and L rebuilding during and after damage; local available fuel/D contacts; founder turnover; and continued protective organization. A positive interaction without these causal flows is insufficient. Conversely, a confirmed healthy-lineage growth margin without a response benefit would show that self-renewal alone does not supply the proposed insurance return.

## A future conditional-kinetics alternative

If the bounded copy-rate test identifies excess L allocation rather than insufficient raw copying, a minimal SAME-catalogue alternative is physical product inhibition: multiply both directions of each L-production channel by h_i=exp(-kappa n_L(i)). Neighbors of target i are unchanged by D_i <-> L_i, so this prefactor is identical on the paired edge and preserves its affinity. Production is unchanged where n_L=0, while existing nearby L slows further local L production. No group, type count, or observed success enters the rule.

That would be an added physical kinetic premise, not an already verified feedback. It must retain enough L renewal and formation, count its fuel use, and survive the same factorial ablations. It is not part of the next test and should not be introduced simultaneously with the copy-rate change.

Even a favorable next result would establish conditional response insurance in a self-produced protective field. Since either R type remains chemically self-sufficient, further work is still needed on variable collective inheritance, reproduction, repeated individuation, informational regulation and aging. Energetic persistence would not establish the moral legitimacy of coercive convergence; the intrinsic standing of participants remains independently relevant.
