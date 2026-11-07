# Complementary chemistry does not require separate taxa

This note compares the proposed productive-role grammar, not an executed spatial study. The earlier [two-role proof](../active-phase/productive-roles-bridge.md) remains a homogeneous ideal-concentration result. Its positive reverse rates, thermal appearance and fixed copying variation are retained throughout the comparisons below. No previous source or report is modified.

## 1. A pathwise equivalence the chemistry cannot distinguish

Take the exact finite spatial count process with template states X/Y, their energies and binding interactions, source/loss channels, parent-dependent copying-error matrix and paid copy reactions. Relabel X as phenotype 0 of genotype G, and Y as phenotype 1 of the same genotype G. Keep every physical state and rate. The state-space bijection Phi gives

\[
Q_G(\Phi s,\Phi s')=Q_{XY}(s,s').
\]

Initial laws, work, heat, material counts and ancestry readouts map identically. Consequently every path event has the same probability. The phenotype inherited at a paid birth follows the same matrix; no free post-birth switch has been added. This is a compatible one-genotype interpretation of the supplied chemical states, not a derivation of a real developmental program. The model contains no genome that could distinguish these interpretations. If a later model adds interconversion between phenotypes, both directions and their actual energy/work must be supplied.

Thus separate chemical roles, even heritable phenotypes, do not establish separate species or genotypes. A new level of individuality could retain developmental differentiation within one lineage. This equivalence survives finite population size, positive immigration, finite reverse rates, spatial localization and any supplied binding energies.

## 2. A generalist control with the same alphabet, cost and mutation

An especially conservative comparison keeps BOTH X/Y physical labels, the same ideal standard energy E, template/working-material inventories, thermal source eta*d, passive loss delta, and matrix Mii=1-m, Mij=m. The labels now have identical catalytic repertoires. Every parent i can execute both transformations, each at HALF the original mobility:

\[
J^{A}_{ij}=\tfrac12M_{ij}c_i(kda-\ell c_jb),\qquad
J^{B}_{ij}=\tfrac12M_{ij}c_i(kdb-\ell c_ja).
\]

Both are signed forward-minus-reverse currents. The same factor multiplies the corresponding microscopic reverse propensity, including n_i(n_i-1) when parent and removed child have the same type. The extra direction is not a free extra full-strength catalyst. The finite chemical alphabet and equilibrium degeneracy are unchanged: replacing two labels with one undeclared chemical state would change mixing entropy and thermal source multiplicity, and is not this comparison.

Let d=M-x-y and a=S-b. Generalist dynamics are

\[
\dot c_j=\sum_i(J^A_{ij}+J^B_{ij})-\delta c_j+\eta d,
\qquad \dot b=\sum_{ij}(J^A_{ij}-J^B_{ij}).
\]

For positive M,S,k,ell,delta,eta and 0<=m<=1/2, these generalists have EXACTLY the specialist fixed point

\[
x=y=u,\ a=b=S/2,\quad
u=\frac{V+\sqrt{V^2+4U\eta M}}{2U},\quad
U=\tfrac12S(2k+\ell),\quad V=\tfrac12kMS-\delta-2\eta.
\]

The fixed point is locally asymptotically stable. Put r=1-2m, c0=ell*u*S/2, t=eta*(M-2u)/u and Q=2u[k(M-2u)+ell*u]. In coordinates x+y, x-y, b, the Jacobian at this point is diagonal, with eigenvalues

\[
-\sqrt{V^2+4U\eta M},\quad
-[c_0+(1-r)\delta+rt],\quad -Q.
\]

All are negative. The symmetric fixed-point concentration, net current of each material direction, template-renewal current, ideal free energy and chemical power equal the specialists'. Along the symmetric submanifold the concentration dynamics agree as well. This is an exact counterexample to necessity of metabolically specialized types for that maintained steady output, even under the stricter half/half allocation budget. It is NOT equality of the finite stochastic path laws: correlations, fluctuations, establishment and extinction statistics can differ. Local stability is not global attraction or spatial formation.

For theta=1, E=6, fuel gap Delta=10, k=1, delta=.1, M=S=2, m=.02, retain ell=exp(-4) and eta=.1*exp(-6). Both controls have x=y=.9413943121, total net copying .1882207548 and chemical power 1.8822075483 per normalized volume/time. Generalist eigenvalues are -1.9005574755, -.0212718663 and -.2531477350. Positive thermal appearance is not called lineage establishment.

## 3. An instantaneous hazard cap is not a synthesis tradeoff

The proposed nearest-neighbor lattice has four target neighbors per parent. At a D target with A material, only the A-to-B forward family is eligible. At a D target with B, only the B-to-A forward family is eligible. At an occupied target the matching reverse family is likewise unique. Solvent does not supply either substrate. Because mutation weights sum to one and the acceptance factor is at most one, giving a promiscuous parent BOTH pathways mobility k still gives

\[
\Lambda^{\rm forward}_{parent}\le4k.
\]

There is no simultaneous two-substrate reaction at one target. Thus the existing four-k forward-hazard cap alone DOES NOT imply kA+kB<=k. A full-strength promiscuous comparator is admissible unless an additional physical catalytic-resource constraint excludes it. Equal standard energy and local detailed balance do not by themselves supply that exclusion.

In the homogeneous example full-strength generalists replace U,V by U_h=h*S*(2k+ell), V_h=h*k*M*S-delta-2eta, with h=1 instead of h=1/2. They retain finite reverse/source/mutation rates and the same equilibrium reference. For positive fuel affinity their positive symmetric fixed point exceeds the half/half value. This is greater maintained template abundance, not automatically greater physical availability or a finite phase unit.

## 4. The physical condition a specialization advantage would need

A genuine fixed allocation budget could follow from a conserved number of constitutively specific active groups: a fraction alpha catalyzes one transformation, the rest the other, and changing specificity requires an explicitly modeled structural transition. This is a physical premise to demonstrate, not a convex cost assumed because differentiation is desired. A promiscuous shared active site, substrate-dependent occupancy or sufficiently fast conformational switching can invalidate it. Binding affinities and switching barriers must enter the same energy/rate scheme, including reverse channels.

Here is an exact, relaxed capacity comparison under that premise. In a specified physical context let a_G,b_G be positive NET productive opportunities per generalist and unit mobility, computed from forward and reverse rates; let a_X,b_Y be the corresponding specialist opportunities after actual partitioning. Hold these opportunities fixed in the allocation problem. With N templates, optimizing a fixed catalytic allocation or a freely allocatable specialist composition gives

\[
J_G^{\max}=\frac{kN}{1/a_G+1/b_G},\qquad
J_S^{\max}=\frac{kN}{1/a_X+1/b_Y}.
\]

These are upper bounds on a balanced A/B cycle current, since the cycle is limited by its slower leg. In this fixed-opportunity allocation problem the specialist optimum exceeds the generalist optimum exactly when

\[
\boxed{1/a_X+1/b_Y<1/a_G+1/b_G.}
\]

The opportunities include finite reverse flux and the supplied variation kernel. They are context dependent; these formulas do not solve the self-consistent phase dynamics. Changing the population composition can change the reverse opportunities, so a global bound needs opportunity envelopes valid throughout the admissible physical states. Copying variation and loss rates can also prevent the optimizing specialist composition. To make a necessity claim, one must bound every admissible generalist below a required maintenance current and construct a specialist regime that actually meets it. For a background-free stationary cycle, that requirement is Jreq=(template losses-spontaneous births)/2. Background A/B conversion or bypass copying modifies this balance and cannot be omitted.

For illustration only, suppose two preexisting spatial contexts offer rates a0 outside and b0 inside. If dilute templates partition faster than renewal, an inner binding advantage g gives p_inside=Vi*exp(g/theta)/(Vo+Vi*exp(g/theta)). A fixed-allocation generalist has a_G=a0*(1-p), b_G=b0*p. Even optimizing p gives

\[
J_G^{\max}=\frac{kN a_0b_0}{(\sqrt{a_0}+\sqrt{b_0})^2},
\]

whereas the perfectly localized specialist bound is kN*a0*b0/(a0+b0), attained in the relaxed problem if trait numbers can realize the optimizing ratio. Symmetric mutation and loss may prevent that ratio; equal opportunities avoid this particular obstruction. This ideal separated-opportunity calculation identifies a physical localization/allocation penalty, at most a factor two in this example. It is not a theorem about our driven lattice. Finite leakage and thermal appearance must be included in actual opportunities and maintenance requirements. Binding energy is not free; creation, transport and heat remain counted.

Crucially, a generalist that reallocates its active capacity to the locally available substrate instead has balanced capacity kN*min[a0*(1-p),b0*p]. Optimizing p reaches the same kN*a0*b0/(a0+b0) as the specialists. Developmental or spatial specialization of one genotype can do the same. The proposed advantage therefore rests on a real inability to retarget capacity, not on a species-count reward.

## 5. Exploiters are not ruled out by omitting their chemistry

In the exact two-family grammar every paid copy necessarily performs one material transformation. Disabling both activities also disables paid reproduction; persistence through thermal or mutant arrival is not evidence that an exploiter can self-renew. This grammar has not demonstrated resistance to a copier that bypasses A/B cycling.

Such a copier needs an explicit additional paired channel, for example Z+D+F reversible to 2Z+W, with its energy, copy cost, finite reverse and positive source retained. In a frozen resident context its Z-count equation has the local form

\[
\dot z=I_Z(z)+[ (1-m_Z)\beta_Z-\delta_Z-r_Z]z+O(z^2).
\]

Here beta_Z is its actual paid copying opportunity, r_Z includes linear reverse-consumption hazards from resident channels, and I_Z(z) retains thermal births and resident-to-Z copying errors, including any source dependence on available precursor. A linearization of TOTAL occupancy must include I_Z'(0); it is not ancestral amplification. The inequality (1-m_Z)*beta_Z>delta_Z+r_Z identifies positive local self-type amplification after separating immigration. A full copying subtree needs the corresponding multitype birth-minus-loss operator; offspring leaving Z are still descendants. The frozen-background diagnostic is not a finite-population survival theorem. In the finite process actual provenance and birth/loss compensators must be measured. If a cheaper bypass exists, it can defeat the productive organization; neither thermodynamics nor the present grammar forbids it.

## Consequence for the next physical comparison

Before an ensemble, include the half/half same-alphabet generalist and, absent an independently modeled specificity budget, the full-strength promiscuous generalist. Preserve node energies, mutation, sources, reverse edges, resources and initial laws. Interpret one-genotype phenotypic relabeling as exact observational equivalence. A positive spatial result would still need to establish self-consistent localization, paid turnover and functional contribution; it would not require separate taxa, imply unlimited biodiversity or establish participant agency.

Run `python3 -W error verify_generalist.py` for the symbolic certificate and numerical example. It runs no trajectories and writes only its new certificate here. The chemistry, kinetic specificity and capacity constraints remain supplied premises, while the displayed equivalence and inequalities are mathematical deductions from them.
