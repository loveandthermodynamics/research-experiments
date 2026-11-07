# A conditional route from reciprocal renewal to functional diversity

The existing R,L threshold supplies a renewing high state, not biodiversity. Here is a specific additional mechanism: a resident's vulnerability changes its own protective production, creating an environmental covariance that can let a differently responding lineage grow when rare. The sign, size and failure cases follow from reaction rates. This note supplies an invasion calculation and a conditional coexistence theorem; it reports no new spatial simulation or availability advantage.

## Physical extension and paid appearance

Retain fixed material \(D+R_++R_-+L=M\), continuing counted fuel/waste supply, paid template copying and paid L production. Both template classes have the same standard energy and productive coefficients a,b. Replace their passive loss mobilities by

\[
d_\pm(L,\xi)=\delta e^{-\alpha L}(1\pm\varepsilon\xi),
\qquad 0<\varepsilon<1.
\]

Here \(\xi\in\{-1,+1\}\) is one physical gate, switching in either direction at rate \(\kappa>0\). Its wells can be isoenergetic and its flip rates symmetric. Every comparison retains this same gate and initial distribution. Each template class has the same marginal exposure, but opposite responses. The multiplier applies to both passive reaction directions. It is a paired barrier change, not a changed well energy or an unpriced destructive intervention. A full joint thermodynamic account retains the gate and its correlations; a symmetric gate's zero heat per flip does not imply zero joint entropy production. Fuel driving remains the source of maintained chemical work.

Actual paid copying with error can produce a previously absent class through

\[
R_i+D+F\rightleftharpoons R_i+R_j+W,
\qquad k_{ij}^+=k_c\mu_{ij},\quad
k_{ij}^-=k_c\mu_{ij}e^{-(E_F-\epsilon_R)/\theta}.
\]

The row-stochastic kernel \(\sum_j\mu_{ij}=1\) preserves total forward copying mobility for each parent. Each realized offspring uses precursor and fuel; no replacement seed or external diversity reward is supplied. The alternative response chemistry and finite mutation alphabet remain supplied possibilities. Actual emergence of a realized variant is different from deriving its response capability from unspecified chemistry.

With rare mutation, a mutant density z initially obeys an equation of the form

\[
\dot z=g_v(t)z+s_\mu(t),\qquad
s_\mu(t)=\mu_{uv}k_c fD_u(t)R_u(t)>0.
\]

The source itself can maintain a label even when its immigration-free growth rate is negative. Appearance, invasion, coexistence and improved renewal are therefore separate tests. Invasion below means growth after removing that source term in the rare-lineage linearization; it is a diagnostic of the local generator, not a controller that restores desirable types. Outbound copying error reduces heritable birth rate and must remain in g_v when it is appreciable. At finite temperature spontaneous births and reverse productive consumption must also be retained. The clean analytical calculation below uses the negligible-source/reverse concentration limit, with a consistent reversible completion available from the prior threshold model.

## Exact invasion identity

Let a resident type u maintain an ergodic positive regime with concentrations \(R_u,L_u,D_u\). Assume its logarithmic abundance has zero long-run growth and the required time averages exist. In the source-free approximation, a rare mutant with productive coefficient \(a_v\) and passive loss \(d_v\) has exponent

\[
\lambda(v\mid u)=a_v\langle D_u\rangle-\langle d_v(L_u,\xi)\rangle
=\frac{a_v}{a_u}\langle d_u(L_u,\xi)\rangle
 -\langle d_v(L_u,\xi)\rangle.                 \tag{1}
\]

The second equality uses the resident's own zero-growth identity. It exposes actual copying or response costs rather than assigning fitness independently of the chemistry.

For opposite responses with equal a,

\[
\boxed{\lambda(-\mid+)
       =2\varepsilon\delta\langle e^{-\alpha L_+}\xi\rangle.}       \tag{2}
\]

The same environmental marginal, \(\langle\xi\rangle=0\), does not force this expression to vanish: L is produced by the resident, so protection and environment can correlate. Conversely, if L is fixed or the protective return is absent, this covariance vanishes and the opposite response alone has zero invasion advantage in this equal-rate model. Equal responses are neutral labels, not functional differentiation.

## A solvable weak-fluctuation criterion

Take a stable positive fixed point \((R_0,L_0,D_0)\) of the unforced renewing model

\[
\dot R=R[aD-\delta e^{-\alpha L}],\qquad
\dot L=bRD-\gamma L,\qquad D=M-R-L.
\]

Write \(h_0=e^{-\alpha L_0}\), \(J_c=aR_0D_0=\delta h_0R_0\), and

\[
J=\begin{pmatrix}
-aR_0&aR_0(\alpha D_0-1)\\
b(D_0-R_0)&-bR_0-\gamma
\end{pmatrix},\qquad \operatorname{tr}J<0,\quad\det J>0.
\]

Assume the resident remains in a persistent high-state neighborhood and its stationary response admits a regular small-\(\varepsilon\) expansion. These assumptions exclude collapse into the low basin. The first-order fluctuation \(y=(r_1,\ell_1)\) satisfies

\[
\dot y=Jy+(-J_c,0)^T\xi_t.
\]

Since a symmetric telegraph gate has correlation \(e^{-2\kappa t}\),

\[
c\equiv\langle y\xi\rangle=(2\kappa I-J)^{-1}(-J_c,0)^T,
\]
\[
c_L=-\frac{J_c\,b(D_0-R_0)}
 {(2\kappa)^2-2\kappa\operatorname{tr}J+\det J}.
\]

Expanding (2) gives the conditional invasion proposition

\[
\boxed{\lambda(-\mid+)=K\varepsilon^2+o(\varepsilon^2),\quad
K=\frac{2\alpha\delta h_0J_c\,b(D_0-R_0)}
 {(2\kappa)^2-2\kappa\operatorname{tr}J+\det J}.}                 \tag{3}
\]

Gate symmetry gives the same coefficient for the opposite invasion. The denominator is positive. Consequently opposite types mutually invade for sufficiently weak nonzero fluctuation if \(D_0>R_0\); the sign reverses if \(D_0<R_0\). At equality the second-order test is inconclusive.

The physical distinction is precise. When \(D_0>R_0\), loss of resident R reduces L production because \(\partial_R(bRD)=b(D_0-R_0)>0\). Protection is weakest when the resident is vulnerable, giving its differently responding counterpart an invasion opportunity. When \(D_0<R_0\), losing R frees enough precursor to increase L production; that compensatory feedback can favor the resident instead. Renewal does not universally fund diversity.

For the existing benchmark \(a=b=\delta=\gamma=1,\alpha=4,M=.95\), the limiting high state has \(R_0=.2328476111,D_0=.5817040018\). At \(\kappa=.25\), equation (3) gives \(K=.18921510545>0\). The same reaction family at M=1.5 has \(R_0=.8953560623,D_0=.3190133768\) and \(K=-.12055996515\). These are algebraic coefficients, not measured finite-amplitude growth rates. [diversity_invasion_algebra.py](diversity_invasion_algebra.py) checks the scalar roots, Jacobian, covariance resolvent and both coefficients; [diversity-invasion-algebra.json](diversity-invasion-algebra.json) retains the results. No trajectory or parameter scan was run.

Finite upkeep is already present through \(\gamma L\), actual productive currents, material competition and counted supply. Additional trait costs or copying leakage can overwhelm (3). For example, a small leakage probability m from the mutant's own copies reduces its exponent by approximately \(m a\langle D_u\rangle\), before other rate changes. A positive covariance term must exceed such physical losses; more mutation is not automatically better.

## The precise gate from mutual invasion to coexistence

**Conditional coexistence proposition.** For the source-free concentration process with the two classes and the gate, suppose a compact set \(\mathcal K\) is forward invariant under both gate vector fields and satisfies \(R_++R_-\geq s_{\min}>0\). Assume smooth bounded per-capita rates, an irreducible finite gate chain, and positive missing-type invasion rate for **every** ergodic stationary law on either pure-type boundary in \(\mathcal K\). Then both initially positive types are stochastically persistent in the occupation-measure sense: time averages place arbitrarily little mass sufficiently near either type's extinction boundary. This does not assert a pointwise abundance floor or unique mixing law. The inference uses [Benaïm, Stochastic Persistence, §3.3 and Theorems 6.1 and 4.4](https://arxiv.org/pdf/1806.08450), which permits auxiliary feedback variables such as L.

For this model the average-Lyapunov function can be taken as

\[
V=-\log(R_+/M)-\log(R_-/M).
\]

Its generator is \(\mathcal LV=-(g_++g_-)\), continuously extending to the boundary; gate jumps do not change V. On a pure-type stationary boundary law the resident's mean g is zero, while the missing type's mean g is positive. The high-state lower bound excludes the all-extinct boundary, where the current subcritical chemistry would instead fail the condition. Positive invasion must be checked on all relevant boundary laws, not one favorably selected resident trajectory. The broader weighted-invasion criterion also appears in [Schreiber, Benaïm and Atchadé, Journal of Mathematical Biology 62, 655–683 (2011)](https://link.springer.com/article/10.1007/s00285-010-0349-5).

A route to certifying the missing high region is available. With \(S=R_++R_-\), the aggregate equations have the unforced R,L vector field plus bounded forcing \(-\varepsilon\delta hS m\xi\), where \(m=(R_+-R_-)/S\in[-1,1]\). Around a stable high equilibrium solve \(J^TP+PJ=-I\). A sufficiently small ellipsoid of \((S-R_0,L-L_0)^TP(S-R_0,L-L_0)\) has inward derivative for all \(|m\xi|\leq1\) if the forcing bound is smaller than its restoring margin. This is a local proof strategy with computable constants; no numerical admissible amplitude or spatial invariant-region certificate is supplied here. Likewise, a contraction certificate or direct analysis must establish the required resident boundary laws and the validity of the response expansion.

This theorem is usable once those conditions are proved; mutual invasion need not remain a permanently inconclusive observation. It nevertheless concerns large-copy concentrations with environmental randomness. A finite population can lose either rare type through demographic events, even when its invasion exponent is positive. A weak selection coefficient may coexist with very small founder numbers or long excursions near extinction. Thermal/paid immigration then changes the question again: recurrent labels are not proof of self-sustaining lineages.

## Transfer to actual local rates

For a rare spatial lineage, use the resident's physical state \(X_t\), including local D, fuel, L and gates. Its linearization has the form

\[
\dot z=A_v(X_t)z+s_\mu(X_t).
\]

Off-diagonal entries contain actual local copying into adjacent precursor/fuel targets and transport. Diagonal entries contain passive loss, reverse productive loss and transport outflow. Source-free typical growth is the top Lyapunov exponent of \(\dot z=A_v(X_t)z\). An eigenvalue of the averaged matrix or an annealed mean-lineage calculation need not equal this typical exponent. If every column sum of A is at most \(-\eta<0\), then \(d\sum_i z_i/dt\leq-\eta\sum_i z_i\): local kinetics rule out invasion despite continuing mutation supply.

This gives a concrete next diagnostic for the same chemistry: reconstruct local birth opportunities and paired loss hazards conditioned on gate state, and test whether their renewal-induced covariance gives positive rare-lineage growth in both resident backgrounds. Deleting the response difference should remove that advantage. Deleting the L return or the variant's own paid L production tests whether the result is coupled regeneration or a lineage living from everybody else's common protective pool. A global L subsidy can support a variant without its improving the organization, and it does not identify a localized individual. Those alternatives must be retained as causal outcomes.

Positive invasion and conditional coexistence establish functional retention only within the supplied response grammar. To claim improved maintained organization as well, the same local model must show that differentiated survivors actually replenish shared/neighboring renewing capacities, and that this benefit survives matched work, material and gate accounting. The existing exact finite-fuel availability theorem supplies a benchmark for that stronger return; it is not automatically inherited by equation (3).

## Continuous traits and the obstruction to unlimited diversity

If response traits fill an interval, let \(d_s=\delta h(1+\varepsilon s\xi)\). With equal productive coefficients, the same expansion gives

\[
\lambda(v\mid u)=C\varepsilon^2u(u-v)+o(\varepsilon^2),\qquad C=K/2.
\]

The response-insensitive generalist u=0 maintains the unforced high state, and \(\lambda(v\mid0)=0\) exactly in this equal-rate source-free model. For C>0 the generalist can invade a responding resident to second order. The two-opposite-type result therefore does **not** derive continuous-trait evolutionary branching; the finite grammar omitted that intermediate capability.

A genuine catalytic tradeoff would have to be measured or derived. To expose the required sign rather than hide it in a diversity reward, suppose actual paired copying mobilities were \(a(s)=a_0[1+\chi\varepsilon^2s^2]\). Then the local expansion becomes

\[
\lambda(v\mid u)=\varepsilon^2\{A(v^2-u^2)+Cu(u-v)\}+o(\varepsilon^2),
\quad A=\chi\delta h_0.
\]

The usual local conditions of attraction toward u=0 but invasibility there would require \(0<A<C/2\). This is a precise possible branching test in kinetic coefficients, not an established biochemical tradeoff or proof of eventual polymorphism. If no such supplied/measured tradeoff exists, the current derivation must stop at its finite-response conditional result.

No unbounded biodiversity optimum follows. Redundant response maps add labels without new insurance; finite material limits viable copy numbers and creates founder/extinction costs. A model requiring at least n copies of each retained type in a population of N can sustain at most \(\lfloor N/n\rfloor\) such types. Additional distinct maintenance functions can also consume finite supply. These are conditional resource bounds, not grounds for assigning moral worth by productivity.

## What the readout would mean

Robustness means retaining a stated function under a specified disturbance range. Response insurance means different vulnerabilities leave some regenerative capacity when another class fails. Antifragility requires a specified outcome to improve as a stated measure of disturbance dispersion increases; that stronger claim is not supplied by survival or coexistence alone. This distinction follows the dispersion-sensitivity definition in [Taleb and Douady, Mathematical Definition, Mapping, and Detection of (Anti)Fragility](https://arxiv.org/abs/1208.1189). Equation (3) can make local invadability increase with weak fluctuation amplitude; it does not establish that total organization performance or available energy does so.

Neither R occupancy nor type count is negentropy by definition. Actual available free energy uses a declared common equilibrium reference and the full state distribution: \(B=E[H]-\theta\mathsf H(p)-F_{\rm eq}\). Increased maintained energetic capacity must exceed any entropy difference, with fuel, protective material, gate correlations and interactions included. A finite-stock system has a finite availability budget; an open renewing state must count supply work. No dissipation-only or Hamiltonian least-action principle replaces these conditions.

Finally, even a diverse, renewing and efficient whole can be coercive. The Unitology defense still requires independently consequential participants whose interests and development matter beyond total energetic performance. Demonstrated physical function is not evidence of subjecthood; efficiency cannot authorize incorporation. Informed voluntary fusion and sacrifice remain possible, and dependent participants and outsiders retain standing without having to win an invasion or productivity contest. This ethical premise complements the kinetic theorem rather than being smuggled into its rates.
