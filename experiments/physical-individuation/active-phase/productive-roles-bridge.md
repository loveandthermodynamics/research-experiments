# Productive roles, paid variation and a possible phase-unit bridge

This is a proposed chemical extension, not a result of the current conserved-enzyme spatial generator. It supplies an exact homogeneous coexistence example and a separate finite-population balance-law obstruction. Neither proves spatial formation, lineage-level biodiversity or an ethical individual.

## A fixed reaction grammar with two productive transformations

Let A and B be soluble and phase-separating forms of one working material. D supplies template material; X and Y are templates with different catalytic activities. For each offspring j in {X,Y}, supply the reversible channels

\[
X+D+A+F\rightleftharpoons X+R_j+B+W,\qquad
Y+D+B+F\rightleftharpoons Y+R_j+A+W.
\]

Also retain X reversible to D and Y reversible to D. A fixed mutation matrix has diagonal entries 1-m and off-diagonal entries m, with 0<m<=1/2. Multiply BOTH directions of each corresponding channel by that entry. A wrong-type offspring is an actual synthesis consuming D and fuel; it is not a replacement supplied by an observer. X copies consume A and supply B; Y copies consume B and return A. The catalytic functions and finite two-type alphabet are premises. This is metabolic cross-feeding, unlike the previous two types that differed only in vulnerability.

For a same-type channel the reverse propensity contains n_X(n_X-1), or n_Y(n_Y-1), because two distinct template instances are reactants. For a cross-type channel it contains n_X*n_Y. Use the matching parent, product, precursor and working-material factors in both directions; mutation does not create an additional unpaired transition. The indistinguishable same-type template instances require a declared observer convention for assigning ancestry, but do not change the physical count process.

There are three conserved material inventories in a closed vessel: D+X+Y, A+B, and F+W. Every copy adds one template, consumes one D, interconverts one A/B and one F/W. Every template loss returns its unit to D. An open fuel bath keeps F/W activities fixed and must carry an explicit work ledger.

## One thermodynamics, including partitioning

For an ideal dilute example choose temperature energy theta, template standard energy E, zero D/A/B standard energies, and fuel chemical-potential gap Delta. With equal standard X/Y energies and equal fuel/waste activities, let

\[
\eta/\delta=e^{-E/\theta},\qquad
\ell/k=e^{-(\Delta-E)/\theta}.
\]

Here delta is passive template loss, eta its thermal reverse; k and ell are forward and reverse copying coefficients after fixing bath activities. Both are positive. Standard-energy heat per forward copy is Delta-E, and passive template loss releases E. Mixing contributions enter chemical potentials and free energy; they must not be omitted from a coarse thermodynamic balance.

For a phase-separating extension, use one nonideal free energy G for A/B interactions, interfaces and template binding, then derive every chemical potential from G. An elementary reaction can use paired transition-state fluxes proportional to exp(sum reactant chemical potentials/theta) and exp(sum product chemical potentials/theta), with the same catalytic prefactor. Alternatively, on the microscopic spatial lattice use a symmetric edge mobility and a rate ratio exp[(work-deltaH)/theta]. Template binding and the phase-material rearrangement then contribute to copying deltaH too. Ideal mass-action rates cannot simply be transplanted into a nonideal phase.

Diffusive fluxes proportional to minus concentration times chemical-potential gradient, together with these paired reactions, give

\[
\dot G=P_{\rm chemical}-\theta\sum_r\int(J_r^+-J_r^-)
 \log(J_r^+/J_r^-)\,dV
 -\sum_i\int M_i c_i|\nabla\mu_i|^2\,dV.
\]

The last two terms are nonnegative dissipation. No diversity or group reward appears. Finite fuel eventually runs down; a maintained stationary state below assumes the explicitly supplied fuel bath.

[Kirschbaum and Zwicker (2021), equations 6, 8, 10-11 and the segregated-enzyme analysis](https://arxiv.org/html/2103.02921v1), provide the primary thermodynamic precedent: fuel-driven conversion and molecular partitioning can sustain different reactions inside and outside a condensate. Their model does not contain the copying grammar proposed here.

## Exact finite-reverse coexistence in an ideal reduction

This subsection assumes a homogeneous ideal vessel, not a droplet. Write x,y,a,b,d for concentrations, with d=M-x-y and a=S-b. The four net channel currents are

\[
\begin{aligned}
J_{XX}&=(1-m)x(kda-\ell xb),&J_{XY}&=mx(kda-\ell yb),\\
J_{YX}&=my(kdb-\ell xa),&J_{YY}&=(1-m)y(kdb-\ell ya).
\end{aligned}
\]

Thus

\[
\dot x=J_{XX}+J_{YX}-\delta x+\eta d,\quad
\dot y=J_{XY}+J_{YY}-\delta y+\eta d,\quad
\dot b=J_{XX}+J_{XY}-J_{YX}-J_{YY}.
\]

For any positive M,S,k,ell,delta,eta and 0<=m<=1/2 there is a symmetric physical fixed point x=y=u, a=b=S/2, d=M-2u>0. Define

\[
U=\frac{S(2k+\ell)}2,\quad
V=\frac{kMS}2-\delta-2\eta,\quad
u=\frac{V+\sqrt{V^2+4U\eta M}}{2U}.
\]

**This fixed point is locally asymptotically stable.** The symmetric concentration eigenvalue is -sqrt(V^2+4U eta M). For the antisymmetric modes x-y and b-S/2, put r=1-2m, c0=ell*u*S/2, t=eta*d/u and Q=2u(kd+ell*u). Their Jacobian is

\[
J_-=
\begin{pmatrix}
-[c_0+(1-r)\delta+rt]&-rQ\\
\delta-rc_0-t&-Q
\end{pmatrix}.
\]

Its trace is strictly negative and

\[
\det J_-=Q[\delta+4m(1-m)c_0]>0.
\]

This proves local stability with a positive reverse copying rate, positive thermal appearance and fixed copying error. It does not claim uniqueness or global attraction. In the negligible-source limit eta=0 a positive symmetric state requires kMS/2>delta; positive eta removes an exact extinction threshold. It is therefore essential to compare a driven high state to the thermal baseline, rather than counting mere presence as establishment.

One reproducible dimensionless example is theta=1, E=6, Delta=10, k=1, delta=.1, M=S=2 and m=.02. Then ell=exp(-4), eta=.1*exp(-6), u=.9413943121 and d=.1172113758. Thermal equilibrium without the driven copy channels has each template concentration .0049330487. Setting Delta=0 while retaining the reversible copy channels gives that same thermal fixed point and zero net copying. Total net copying in the driven example is .1882207548 per normalized volume/time, and chemical power is 1.8822075483 per normalized volume/time. The eigenvalues are -1.9005574755 and -.1372098006 +/- .0826465187 i. These are algebraic fixed-point calculations, not a spatial simulation or a finite-population persistence proof.

Starting with X but no Y, an actual forward X-to-Y copy has propensity proportional to m*k*x*d*a, whereas thermal Y appearance has propensity eta*d. Their instantaneous forward-source ratio is m*k*x*a/eta. Both sources remain present. A lineage experiment must record which channel actually created Y and its paid descendants; a positive stationary Y count does not identify its origin.

## A finite-state role-necessity bound that permits thermal sources

This separate result applies to any finite continuous-time count process with the stated stoichiometry, whether spatial or well mixed. It needs no ideal rates or zero immigration. Let C_X,C_Y be signed forward-minus-reverse copy counts in the two reaction families; C_0 is every OTHER signed A-to-B conversion, including the current spatial model's passive A/B channel and its conserved-E driven channel if retained. Neither may be omitted from the background bound. Let N be total templates, I spontaneous/external template births and L template losses. Pathwise,

\[
B_T-B_0=C_X-C_Y+C_0,\qquad
N_T-N_0=C_X+C_Y+I-L.
\]

Diffusion and partitioning do not affect these global balances. Mutation changes offspring identity but not either identity. Suppose the instantaneous background net current obeys |j0(state)|<=Lmax, source intensity i(state)<=Imax, and passive template-loss intensity l(state)>=delta_min*N, with delta_min>0. All rates and bounds refer to the actual finite model, including its reverse chemistry and binding energies.

Block BOTH directions of either productive reaction family. Stationarity gives E[jX-jY+j0]=0, so the mean net total copying cannot exceed Lmax. The template balance then proves

\[
\boxed{\quad E[N]\le\frac{L_{\max}+I_{\max}}{\delta_{\min}}.\quad}
\]

Consequently P(N>=n*) is at most (Lmax+Imax)/(delta_min*n*). Large source/background terms make this bound uninformative, which is a real limitation, not a reason to omit them. For a finite observation time T, A+B=S_count and initial template count N0 give the conservative bound

\[
\frac1T\int_0^T E[N_t]dt\le
\frac{L_{\max}+I_{\max}}{\delta_{\min}}
+\frac{S_{\rm count}+N_0}{\delta_{\min}T}.
\]

This follows from the same path identities and compensators; it accounts for stored A/B material and initial founders. No absorbing state, infinite lifetime or source-free invasion approximation is used. With no background A/B conversion, neither one-role ablation can sustain a positive net copying current, although thermal templates remain. In a closed large finite population eventual rare losses still occur, and positive immigration can restore occupancy; neither fact invalidates the balance law or proves ancestral continuity.

## What would actually connect this to finite phase units

Molecular interactions must physically enrich the B-destroying role inside the B-rich phase and generate B supply outside it. In the prior thin-interface radial reduction, use ACTUAL net currents to obtain

\[
\dot R=A_R/R-C_R/R^2-B_RR,\qquad
B_R=s_{\rm net\ destruction}^{\rm in}/(3\Delta c).
\]

The productive role concentrations affect A_R and B_R through local rates and molecular partitioning, not an observer's inside label. A positive radially stable branch requires

\[
0<B_R<4A_R^3/(27C_R^2),\qquad A_R,C_R>0.
\]

Zero turnover has no stable finite radius in that reduction, while excessive turnover destroys the finite branch. Neither the sink activity nor diversity is maximized. This condition presupposes phases, scale separation and the flux law; it does not prove nucleation, shape stability or multiple-droplet stability. The coexistence calculation and radius condition have NOT been joined into a solved spatial model. In particular, spatially averaged coexistence could be supported by a globally shared bath rather than by each renewing unit. Actual coupled formation, local turnover and paired role ablations remain necessary evidence.

The strongest conclusion here is necessity of complementary productive TRANSFORMATIONS within this supplied grammar. A generalist performing both transformations, or one genotype with developmental specialization, can satisfy it without separate taxa. An exploiter with a cheaper bypass copying pathway can violate the assumptions and invade. These are required comparison models, not exceptions ruled out by thermodynamics. Even a verified physical phase unit can coercively absorb other agents. The Unitology critique still needs the independent participant criterion: the whole must preserve consequential capacities of its constituents, while permitting voluntary fusion or sacrifice. Metabolic efficiency alone does not supply that premise.

## Reproduction and scope

Run `python3 verify_productive_roles.py` in this directory. The certificate reconstructs the four-channel ODE and Jacobian symbolically, verifies material conservation and the positive determinant identity, and recalculates the numeric example. It does not run trajectories or modify previous evidence. SymPy is required for the symbolic checks. Mathematical results above are conditional derivations; primary literature supports the thermodynamic framework, not this new role model or a universal biodiversity law.
