# Joined productive-phase generator: design for review

This is a new proposed chemistry. It does not modify or rescue the frozen active-phase study, whose full regime failed dense formation. The purpose is to join productive copying to material supply and consumption in one spatial generator, while retaining reversible chemistry, thermal appearance and alternative explanations. No simulation, parameter ensemble or claimed spatial result accompanies this design.

## Two local layers and one Hamiltonian

Each site of a periodic square lattice has a material slot S/A/B and a template slot D/X/Y. S is solvent, A soluble working material, B sticky working material, D template precursor, and X/Y are supplied alternative catalytic capabilities. The layers are a declared coarse-graining: their molecules can share a cell, with excluded volume within each layer. Their conserved inventories are `N_A+N_B=M`, `N_S=Nsite-M`, and `N_D+N_X+N_Y=Nsite`.

Write b_i=1 for B, x_i=1 for X and y_i=1 for Y. Supply

`H = epsilon_B*sum b_i - J_BB*sum_<ij> b_i*b_j + E_X*N_X + E_Y*N_Y - g_X*sum x_i*b_i - g_Y*sum y_i*b_i`.

The bond sum counts each undirected nearest-neighbor edge once. A, S and D have zero reference energy. Positive J_BB attracts sticky material. A candidate g_Y>g_X makes Y favor B more than X; negative g_X means repulsion. These molecular preferences are supplied, not deduced from least action. No observer-detected phase, inside/outside label, component size, genealogy or group reward enters H or a rate.

Material and template molecules each exchange across nearest-neighbor edges, accounting for the complete change in H. The old conserved enzyme and its independent fueled B-destruction channel are absent. This is an explicitly different reaction alphabet, not an unnoticed deletion from old evidence.

## Copying and its reverse occur at the same target

For every oriented nearest-neighbor parent j and target i, j!=i, and offspring r in {X,Y}, include

`X_j + D_i + A_i + F <-> X_j + Rr_i + B_i + W`,

`Y_j + D_i + B_i + F <-> Y_j + Rr_i + A_i + W`.

The parent at j is unchanged. The target precursor becomes an actual offspring while its material is transformed. The fixed mutation matrix has diagonal 1-m and off-diagonal m, multiplying both directions of the corresponding local channel. Every reverse uses the actual matching parent, offspring and product-material occupancies. It may consume any chemically matching target molecule, not just a historically recorded child. Ancestry must never prohibit a physical reverse event.

Same-type parent and offspring occupy distinct cells. The oriented microscopic channels therefore identify two distinct reactants without a special count-level factor 1/2. Summing over equivalent parent/target arrangements gives the usual ordered falling-factorial count. Parallel helper channels remain separately identified in the event ledger.

Also supply passive `X_i <-> D_i` and `Y_i <-> D_i`, plus an optional weak paired `A_i <-> B_i` background with mobility kappa. The background remains in every material balance; zero kappa is a separate structural limit, not a silently imposed claim of biological isolation.

F/W are homogeneous ideal chemostats with chemical-potential gap mu. Every forward copy receives work +mu and converts one fuel to waste; every reverse receives -mu and synthesizes one fuel. All other channels have work0. Given an elementary edge with reverse-symmetric mobility a, use

`q(s,s') = a*min(1, exp[(work-(H(s')-H(s)))/theta])`.

Thus the exact local rate ratio is `exp[(work-deltaH)/theta]`. Copy mobility is k_X or k_Y times the paired mutation weight. There is no independent per-event energy exception for a mutant. Mobility0 disables both paired directions. Every enabled finite-energy pair retains a reverse; representational underflow must raise an error rather than silently remove one. At mu=0 every channel shares the conditional Gibbs reference `exp(-H/theta)`.

The event ledger must record target, parent, parent/offspring phenotype, direction, microscopic before/after state, deltaH, chemical work, heat `Q=work-deltaH`, and signed fuel conversion. Initial material/templates, preparation, lattice and reservoir apparatus are supplied resources. The work ledger counts chemical input, not apparatus construction or initial molecular synthesis. Tags can record copying-instance provenance and conserved material carriers separately, but neither changes rates or establishes functional identity. The finite thermodynamic state is the occupancy configuration. Unbounded cumulative counters and newly allocated copying-instance IDs are observational histories, not additional thermodynamic microstates whose multiplicity changes H or entropy production.

## Exact balances: turnover and material density are different

Let C_X and C_Y be signed forward-minus-reverse copy counts classified by the parent's productive transformation. C_0 counts every non-copy A->B conversion, including the reverse sign for B->A. Let I be all non-copy template births and L passive template losses. Pathwise,

`B_T-B_0 = C_X-C_Y+C_0`,

`N_T-N_0 = C_X+C_Y+I-L`,

`W_chemical = mu*(C_X+C_Y)`.

Here N is the total template count; diffusion changes neither global inventory. Mutation changes offspring identity without changing these balances. At a stationary law, with mean currents J_X,J_Y,J_0 and mean intensities Lbar,Ibar,

`J_X=(Lbar-Ibar-J_0)/2`, `J_Y=(Lbar-Ibar+J_0)/2`,

`P_chemical=mu*(Lbar-Ibar)`.

Both net productive roles are positive exactly when `Lbar-Ibar>|J_0|`. This is a statement about net currents under the declared stationary law, not every individual event or a guarantee of taxonomic diversity.

If kappa=0 and one whole copy family is disabled, the remaining chemical work is a bounded state difference: X-only gives `W=mu*(B_T-B_0)` and Y-only gives its negative. All enabled edges then satisfy detailed balance for `H_eff=H-mu*B` or `H_eff=H+mu*B`, respectively. Sustained stationary work and copy currents vanish, although thermal templates can remain. A positive passive A/B background reopens a one-family fuel cycle; it must be bounded rather than ignored.

The proposed default is explicitly kappa=0; positive kappa remains an optional, separately bounded background. With both families active, choose a target initially D/A and two distinct nearest-neighbor parents, one X and one Y. X copying at the target, passive offspring loss, Y copying at the same target, and passive offspring loss return the physical occupancy to its start while receiving work2mu. The cycle log-rate affinity is2mu/theta, independently of intermediate energies. For the nondegenerate finite occupancy model with M>0, both positive copy mobilities, positive template birth/loss and local diffusion, these states are accessible; the chain is irreducible on the stated inventories. Thus mu!=0 breaks detailed balance and yields strictly positive stationary entropy production. For mu>0 and kappa=0, the stationary net family currents are equal and strictly positive. This proves a sustained chemical cycle under the supplied grammar, not a dense spatial unit. At mu=0 the full physical generator satisfies Gibbs detailed balance.

## Explicit source and bypass bound

With a common passive template mobility delta and `E_r>max(0,g_r)`, every template-loss step is downhill and its rate is exactly delta. Its reverse at a D site is `delta*exp[-(E_r-g_r*b_i)/theta]`. Define

`c0=exp(-E_X/theta)+exp(-E_Y/theta)`,

`cB=exp[-(E_X-g_X)/theta]+exp[-(E_Y-g_Y)/theta]`.

The total thermal birth intensity is exactly

`I(s)=delta*[c0*(Nsite-N(s))+(cB-c0)*N_(D on B)(s)]`.

It is not a uniform negligible source: affinity can substantially increase thermal Y appearance in B. Paired Metropolis A/B background gives `|j0(s)|<=kappa*M`. For cB>=c0, disabling either whole copying family therefore yields the stationary expectation bound

`E[N] <= [kappa*M/delta+c0*Nsite+(cB-c0)*M]/(1+c0)`.

For arbitrary affinities, replacing `(cB-c0)` by its positive part gives a conservative valid bound. An instantaneous count can exceed an expectation bound. For a finite observation interval, add `(M+N_0)/[delta*(1+c0)*T]` to the conservative time-average version. A random initial N_0 is replaced by its expectation. Additional non-copy births or A/B channels must be added to the relevant bounds.

The illustrative values E_X=E_Y=6, theta=1, g_X=-2.5, g_Y=2.5, delta=.02, kappa=.0002, Nsite=256 and M=77 give a bound of about 3.98 templates. This is a constraint for a possible new regime, not a success criterion selected from outcomes or evidence of a high population. None of these candidate physical values has been frozen for a study.

## A physical supply/consumption criterion

Net renewal alone does not supply a dense phase. For any region Omega used only for analysis,

`d E[B_Omega]/dt = E[j_X,Omega-j_Y,Omega+j_0,Omega+J_diffusive_in,Omega]`.

Reactions are assigned by their target site; diffusive flux counts B exchanges across the region boundary. A stationary B-rich region needs its actual internal net sink to be balanced by material influx, while the exterior supplies that influx. A globally balanced, rapidly cycling dilute state also satisfies the global equations. Formation cannot be inferred from them alone.

There is an exact residence accounting test. Starting at zero B, let every material carrier's B episode end on its next B->A conversion, irrespective of how it moves. Then `integral_0^T N_B(t)dt` equals the sum of all B-episode durations truncated at T. At stationarity, if the episode intensity and mean duration exist, `E[N_B]=J_AtoB*tau_B`, by occupation/entry accounting; independence or exponential episode lengths is unnecessary. A target mean B fraction phi therefore requires enough formation throughput AND residence: `J_AtoB*tau_B>=phi*Nsite`. Large fuel use or repeated A/B conversions alone can leave this product small, as the earlier dilute candidate illustrates. A dense phase additionally requires spatial segregation, not merely this global inventory condition.

For equal copy mobilities k, a B target is consumed by a Y forward event only if its template slot is D, or an X reverse event only if its slot contains the matching offspring. Those target conditions are mutually exclusive. With four neighbors and normalized mutation weights, its total B->A conversion hazard is therefore at most `kappa+4k`, not kappa+8k. For positive kappa a conservative instantaneous lower bound is `kappa*exp[-max(0,-epsilon_B+4J_BB+max(0,g_X,g_Y))/theta]`, supplied by passive B->A alone. No positive instantaneous lower bound from that channel is asserted at kappa=0. Residence concerns a carrier's chemical B state, not replacement of the carrier or continuity of a droplet.

For interpretation only, a local factorized closure with slowly varying template fractions gives `b_dot≈Kplus*(c-b)-Kminus*b`, where c=a+b, Kplus includes X-driven production, reverse Y and passive production, and Kminus includes Y-driven destruction, reverse X and passive destruction. Its retained-material level is `b*=c*Kplus/(Kplus+Kminus)`. Reaching b_min<c requires `Kplus/Kminus>=b_min/(c-b_min)`. These are effective coefficients computed from actual local pair statistics; factorization is not an exact theorem about the lattice, and the ratio does not prove phase separation or a selected size.

The new coupling avoids one specific structural defect: B destruction is no longer performed by a permanently supplied conserved enzyme independently of template renewal. It requires a Y parent, a B substrate and a D target, creates a real offspring, and is countered by X copying that returns material to B. Template loss replenishes D. This permits a balanced supply/sink cycle; it does not guarantee one. Excessive sink current, insufficient D availability, insufficient source current, mixing or slow nucleation can still prevent dense units.

## Local energetics and partitioning must be used, not presumed

For a forward target conversion with n_B neighboring B sites, the exact copying energy changes are

`deltaH_X(r)=E_r+epsilon_B-g_r-J_BB*n_B`,

`deltaH_Y(r)=E_r-epsilon_B+J_BB*n_B`.

The Y formula has no offspring binding gain because its target product is A. These expressions include every affected term of H. Parent binding at another cell is unchanged. A sufficient condition for every forward copy to have unit Metropolis acceptance is mu at least the maximum of these expressions over r and n_B=0..4; otherwise actual acceptance must be retained.

For E=6, epsilon_B=2, J_BB=2.4 and the opposed affinities above, X's maximum cost is 10.5, whereas Y's is 13.6. At mu=12, X forward copying is saturated but the most costly Y step has acceptance `exp(-1.6)` and its reverse can be fully accepted. Replacing this chemistry with uniformly saturated Y copying would be a different generator.

Exchanging X on B with Y on A changes binding energy by g_X-g_Y. Thus the specified affinities energetically favor Y on B and X outside without a phase detector. A fast-diffusion, local-equilibrium approximation would give a relative Y/X sorting factor `exp[(g_Y-g_X)/theta]` between B and non-B environments. Reactive stationary distributions need not equal this equilibrium prediction. Measured partitioning, available D and local net currents are required to show that the intended source/sink arrangement actually occurs.

## Alternatives and next review

Complementary transformations are not equivalent to distinct taxa. A promiscuous generalist can implement both; developmental states can assign the two transformations within one genotype. Because A/B targets are mutually exclusive, two full-k pathways can still satisfy the same 4k per-parent instantaneous forward-hazard bound. That bound alone does not justify halving a generalist's rates. A shared active-site occupancy, extra catalytic-domain synthesis or switching cost requires an additional explicit physical premise and corresponding states/energy accounting. Bypasses and exploiter reproduction must likewise be modeled rather than ruled out verbally.

The ideal homogeneous coexistence theorem and the earlier three-dimensional radial reduction are consistency examples, not a solved version of this two-dimensional model. To demonstrate finite units would still require unseeded spatial formation, late-time and domain-size tests, continuing local currents, turnover, and prospective paired controls. Inheritance, regulatory memory, aging, recursive individuality and the normative participant criterion remain separate. No observer reward for diversity, phase size or a whole's propagation has been introduced.

This specification is for independent physics and theory review before code or ensembles. Parent/target reaction coordinates, reverse-edge multiplicity, thermal origins, all bypass currents and the meaning of generalist controls must be agreed before a frozen implementation and bounded prospective protocol.
