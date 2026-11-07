# From produced protection to a renewing organization

This is a proposed minimal mathematical bridge, not a new spatial experiment or a universal law of individuation. The present checkpoint establishes produced protection and copying benefits, while its post-formation comparison does not yet establish that ongoing protective production is necessary. The following model identifies a sufficient mechanism and an exact threshold for that next claim. It also identifies what remains between reciprocal maintenance, individuality, biodiversity and the Unitology defense.

## Supplied physical premises

Consider a homogeneous reactor of fixed volume with precursor D, active template R and protective material L. The finite material inventory is

\[
D+R+L=M.
\]

The productive reaction grammar is supplied:

\[
R+D+F\rightleftharpoons 2R+W,
\qquad R+D+F\rightleftharpoons R+L+W.
\]

R supplies the catalytic/template function; L supplies a barrier effect. Productive channels use mass action. Passive R-to-D conversion has the paired mobility \(\delta e^{-\alpha L}\), with \(\alpha\geq0\); L-to-D conversion has mobility \(\gamma>0\). The L-dependent barrier law is a supplied constitutive assumption, not a derived membrane or an independently resolved encounter mechanism. It multiplies both directions and changes kinetics without changing well energies. There is no diversity, group, care or longevity term in a rate or objective.

All concentrations share a reference concentration \(c_*\). Energies below are per particle with \(\theta=k_BT\), or consistently per mole with \(\theta=RT\). Thus \(\alpha\) has inverse-concentration units; \(\delta,\gamma\) have inverse-time units. Fixed reactor walls, mixing, reaction possibilities and rate constants remain premises. A threshold inside this reactor does not derive its walls.

Choose standard chemical potentials \(\epsilon_D=\epsilon_W=0\), \(\epsilon_R,\epsilon_L>0\), and \(\epsilon_F=E_F>\max(\epsilon_R,\epsilon_L)\). Define

\[
q_R=e^{-\epsilon_R/\theta},\quad q_L=e^{-\epsilon_L/\theta},\quad
k_c^-=k_c e^{-(E_F-\epsilon_R)/\theta},\quad
k_L^-=k_L e^{-(E_F-\epsilon_L)/\theta}.
\]

The net productive currents are

\[
J_c=k_cRDF-k_c^-R^2W,\qquad
J_L=k_LRDF-k_L^-RLW.
\]

The passive currents, oriented toward D, are

\[
J_R=\delta e^{-\alpha L}(R-q_RD),\qquad
J_d=\gamma(L-q_LD).
\]

Consequently \(\dot R=J_c-J_R\) and \(\dot L=J_L-J_d\). These reverse coefficients share one thermodynamics; they cannot be chosen independently. For example,

\[
\log(J_c^+/J_c^-)
=\frac{\mu_D+\mu_F-\mu_R-\mu_W}{\theta},
\quad \mu_i=\epsilon_i+\theta\log(c_i/c_*).
\]

The analogous identities hold for every channel. The barrier factor cancels from the R passive forward/reverse ratio. This uses the standard chemical-work and local-detailed-balance framework of [Rao and Esposito, Physical Review X 6, 041064 (2016), equations 50 and 78](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.6.041064); the particular R,L construction here is a proposal.

## Exact threshold in the negligible-reverse limit

First supply fuel and remove waste through fixed, explicitly counted chemostats \(F=f,W=w\). Let \(a=k_cf\) and \(b=k_Lf\). In the approximation where the thermally uphill and reverse productive currents are negligible,

\[
\dot R=R[a(M-R-L)-\delta e^{-\alpha L}],
\qquad
\dot L=bR(M-R-L)-\gamma L.                 \tag{1}
\]

The approximation is not exact finite-temperature irreversibility. A consistent finite-temperature completion is given below. The triangle \(R\geq0,L\geq0,R+L\leq M\) is invariant.

Set \(c=\delta/a\), a concentration, and \(\beta=b/\gamma\), an inverse concentration. At any positive stationary state,

\[
D=ce^{-\alpha L},\qquad
R=\frac{Le^{\alpha L}}{\beta c}.
\]

Thus every positive stationary state, and no others, corresponds to a positive root of

\[
\boxed{M=F(L)\equiv L+ce^{-\alpha L}
             +\frac{Le^{\alpha L}}{\beta c}.}       \tag{2}
\]

The decisive property is strict convexity:

\[
F'(L)=1-\alpha ce^{-\alpha L}
       +\frac{e^{\alpha L}(1+\alpha L)}{\beta c},
\]
\[
F''(L)=\alpha^2ce^{-\alpha L}
       +\frac{\alpha e^{\alpha L}(2+\alpha L)}{\beta c}>0
\quad(\alpha>0).                                  \tag{3}
\]

The empty stationary state is linearly asymptotically stable exactly when \(aM<\delta\), or \(M<c\). Equality is nonhyperbolic and is not classified by this linear criterion. This is the regime in which unprotected templates cannot grow even at maximal precursor availability.

**Threshold proposition.** For fixed parameters with \(M<c\), equation (1) has a stable positive stationary state, separated from extinction by a saddle, **if and only if**

\[
\boxed{\min_{L\geq0}F(L)<M<c.}                     \tag{4}
\]

A nonempty interval of material inventories satisfying (4) exists **if and only if**

\[
\boxed{\alpha c>1+\frac1{\beta c}
                 =1+\frac{\gamma}{bc}.}            \tag{5}
\]

Proof: \(F(0)=c\), \(F\to\infty\), and \(F'\) is strictly increasing. Its minimum is below c exactly when \(F'(0)=1-\alpha c+1/(\beta c)<0\). A horizontal line strictly between that minimum and c cuts the curve twice. At a positive stationary state the Jacobian is

\[
J=\begin{pmatrix}
-aR&aR(\alpha D-1)\\
b(D-R)&-bR-\gamma
\end{pmatrix},
\quad \mathrm{tr}J=-(a+b)R-\gamma<0,
\quad \det J=abRD F'(L).
\]

The lower-L root is therefore a saddle and the higher-L root is stable. At the minimum they meet in a saddle-node. The Dulac multiplier \(1/R\) gives divergence \(-a-b-\gamma/R<0\) in the positive interior, excluding periodic orbits there. The stable manifold of the saddle supplies a basin threshold, not a universal scalar minimum founder count. Both R and inherited L affect which basin is reached.

Equation (5) makes the maintenance tradeoff explicit: larger protective efficacy and sufficient paid L production can compensate for finite L turnover. At fixed M, the exact criterion is (4); (5) alone only guarantees that some suitable M exists. Production that is too slow, turnover that is too fast, or protection that is too weak removes this bistability. Taking \(\gamma=0\) would instead permit a permanent inherited shield and would defeat the intended ongoing-maintenance test.

This is an internally generated persistence threshold under a supplied reaction grammar and supply regime. It is not obtained by assigning a group reward or selecting a preferred morphology. It also does not show that an arbitrarily small founder can form the high state: a sufficiently small founder lies in the low basin. Formation from zero L under an allowed founder budget requires a separate basin/reachability test.

## A fully reversible instance

Use concentration and time units in which

\[
a=b=\delta=\gamma=1,\quad \alpha=4,\quad M=.95.
\]

Here \(c=\beta=1\). The limiting minimum is at \(L_*=.08355007308\), with \(F(L_*)=.91616510401<.95<1\). The limiting positive states are:

| State | R | L | D | Stability |
| --- | ---: | ---: | ---: | --- |
| Lower branch | .03456303266 | .03058322867 | .88485373868 | Saddle |
| Higher branch | .23284761106 | .13544838716 | .58170400178 | Stable |

Now choose \(\theta=1\), \(\epsilon_R=\epsilon_L=20\), \(E_F=40\), \(f=w=1\), and \(k_c=k_L=1\). All omitted coefficients are consistently \(\eta=e^{-20}>0\). The exact open equations are

\[
\dot R=RD-\eta R^2-e^{-4L}R+\eta e^{-4L}D,
\quad
\dot L=RD-\eta RL-L+\eta D.                         \tag{6}
\]

At stationarity the second equation gives the scalar reduction

\[
L(R)=\frac{(R+\eta)(M-R)}{(1+\eta)(1+R)}.
\]

Substitution into the first equation gives the following three certified roots. The outward-interval calculation in [interval_threshold_certificate.py](interval_threshold_certificate.py), with results in [threshold-interval-certificate.json](threshold-interval-certificate.json), certifies existence, local uniqueness and stability within three root brackets; it does not certify global exhaustiveness.

| State | R | L | Stability |
| --- | ---: | ---: | --- |
| Thermal/low baseline | 3.91619708e-8 | 3.91619650e-8 | Stable |
| Threshold | .03456298566 | .03058319180 | Saddle |
| Renewing high state | .23284761705 | .13544838977 | Stable |

The high-state eigenvalues are approximately -.134702356 and -1.330992885. The low absorbing state of the approximation becomes a small source-maintained baseline. Hyperbolic limiting branches also persist for all sufficiently small consistent reverse coefficients by continuity; the explicit numbers are not a claim about arbitrary temperatures. These are concentration fixed points, not simulation evidence for a finite-noise individual. A finite irreducible molecular-count completion has one stationary distribution and may exhibit long-lived metastable regimes rather than permanently separate deterministic attractors.

## Maintenance and finite resources are actually paid

In the limiting high state,

\[
J_c=aRD=\delta e^{-\alpha L}R>0,\qquad
J_L=bRD=\gamma L>0.
\]

R copying replaces lost active templates; L production replaces lost protection. Neither inventory is merely preserved forever. Every net productive birth consumes one fuel unit. In particular \(J_L/J_c=b/a\): the protective burden is explicit rather than a free parameter in an external fitness score.

For chemostats the fuel chemical-potential difference is

\[
\Delta\mu_F=E_F+\theta\log(f/w).
\]

Let \(G_{\rm mat}\) be the ideal-mixture free-energy density of D,R,L. With entropy-production rate density in Boltzmann units \(\sigma=\sum_\rho(J_\rho^+-J_\rho^-)\log(J_\rho^+/J_\rho^-)\),

\[
\dot G_{\rm mat}=\Delta\mu_F(J_c+J_L)-\theta\sigma.
\]

The finite-temperature high state above has \(J_c=.135448388527\), \(J_L=.135448388573\), and chemical power \(40(J_c+J_L)=10.8358710840\), equal to \(\theta\sigma\) at stationarity. These are normalized per-volume quantities. There is no inference that this state minimizes dissipation or is morally better because it dissipates.

For a closed finite batch instead set \(F+W=K\) and \(\dot F=-(J_c+J_L)\), keeping all the reverse channels. The complete ideal-mixture free-energy density

\[
G=\epsilon_RR+\epsilon_LL+E_FF
 +\theta\sum_{i\in\{D,R,L,F,W\}}c_i[\log(c_i/c_*)-1]
\]

obeys \(\dot G=-\theta\sigma\leq0\). Its equilibrium has \(R=q_RD\), \(L=q_LD\) and \(F/W=e^{-E_F/\theta}\), with the inventories fixing the concentrations. It cannot sustain an indefinitely driven high plateau. A continuous plateau consuming net fuel at least \(j_{\min}>0\) lasts at most its starting available fuel stock divided by \(j_{\min}\). Fuel standard energy alone is not the full free-energy account: mixing and correlations also matter. A maintained stationary state requires continuing supplied work, while the existing finite lattice tests transient persistence from a finite preparation budget.

## The interventions that establish the loop

Four comparisons have distinct meanings. They should use the same well energies, possible species and matched initial physical inventories, removing both directions of any ablated productive channel.

1. **Remove the barrier return, \(\alpha=0\).** For (1) with \(M<c\), \(\dot R\leq(aM-\delta)R<0\); protection-required renewal is impossible. Producing L without its protective return can still spend fuel.
2. **Stop paid L production, \(b=0\), after the same preparation.** With \(\gamma>0\), inherited L decays and the unprotected subcritical R population eventually declines. A horizon shorter than L's effective lifetime cannot establish ongoing production's necessity. At finite temperature a small thermal L baseline remains and must be retained in the comparison.
3. **Stop paid R copying, \(a=0\).** In the limiting model R loss remains at least \(\delta e^{-\alpha M}\) per template, so inherited protection can delay loss but cannot maintain R indefinitely. This separates protection from regeneration.
4. **Perturb and recover.** Convert lost R or L to D, record preparation energy/work, and preserve total material; then test recovery without resets or replacement seeds. Loss below the basin threshold and recovery within the high basin are both predicted. Repeat long enough for founder chemical instances and old L to turn over. Conserved atoms can be recycled; replacing chemical instances is not replacing every atom or demonstrating a persistent subject.

There are exact finite-temperature steady-state bounds for the example (6), so necessity need not rest on an absorbing-zero approximation. With protection removed, any stationary state satisfies \(R\leq\eta M/(1-M)\). With both directions of paid L production removed, stationarity gives \(L=\eta D\leq\eta M\). Hence \(h=e^{-\alpha L}\geq h_{\min}=e^{-\alpha\eta M}>M\), and \(R\leq\eta M/(h_{\min}-M)\). Removing both directions of R copying instead gives \(R=\eta D\leq\eta M\). These bounds are of order \(10^{-8}\), rather than the intact high state's R approximately .233. They establish reciprocal necessity for the high stationary state in this new open model, not an already demonstrated result in the finite-fuel spatial checkpoint. Transition times and legacy protection still require explicit intervention trajectories.

At the high fixed point the passive founder-loss hazards are \(\delta e^{-\alpha L}>0\) for R and \(\gamma>0\) for L, while the population profile is renewed. These are the total loss rates in the negligible-reverse model; in the finite-temperature completion reverse productive reactions add instance loss. This proves continuity of a chemical organization through role turnover. It does not prove heritable organizational memory, adaptation or biological aging.

The crucial existing-checkpoint gap is therefore specific: show a regime in which continued paid replacement of finite-lived L makes a reproducible post-formation difference beyond inherited protection, then show that R copying and the L barrier return are causally necessary for that same restoring regime. Increasing turnover can make this test discriminating; increasing it too far destroys the organization, exactly as the threshold predicts.

## What still separates this from individuation and biodiversity

A restoring R,L state does not require a crisp self-made enclosure. It can occur in a shared well-mixed bath. Consequently it establishes a reciprocal kinetic organization, not a localized individual. To attribute maintenance to a particular assembly, distinguish locally retained protection and descendant renewal from benefits supplied by everybody's common L pool. Equal-inventory isolated-assembly versus shared-bath comparisons, local production/barrier lesions, transport perturbations and founder-turnover assays must identify the causal scale. Localization, a selective interface, enclosure and reproduction are separate results.

The model contains only one active template function and one preset protective product. It has no theorem of biodiversity. A response-diversity advantage additionally needs physically differentiated capacities that can preserve or rebuild a shared renewing capacity under matched disturbances, with that return exceeding fuel, material, founder-dilution, transport and common-failure costs. The existing finite-fuel availability theorem supplies one conditional example; the spatial response comparisons remain inconclusive. Neither more labels nor the high branch here supplies that missing result.

For recursive individuation, the unresolved criterion is stronger still: a maintained organization must generate variable descendant organizations through the same local physical grammar, with a demonstrable causal inheritance of maintenance-relevant variation; a second-level association must renew its lower-level units and pay its own interfaces and regulation. A common pool, persistent cluster or favorable whole-system energy balance cannot substitute for these interventions. A full path-action formulation may describe the resulting generator, but cannot infer its kinetics from entropy production alone; see the separate path-ranking audit.

## A logically valid Unitology defense

An efficient coercive whole can satisfy every physical condition above. Physical individuality, regenerative success, archived information and lower cumulative cost therefore cannot by themselves rule it out.

A defensible conditional argument makes its ethical premise explicit: ongoing participants' interests and capacities to learn, revise commitments, affect collective purposes and pursue development matter independently of their contribution to total energetic performance. If an imposed convergence destroys or appropriates those capacities without valid authorization or adequate justification, thermodynamic success does not cancel that wrong. This premise is not derived from the R,L equations, and R or L molecules are not thereby moral participants.

The premise permits development, voluntary sacrifice and informed voluntary fusion; it does not freeze existing selves or require every constituent to remain separately embodied forever. Where participants continue, their representation and consequential influence must be tested rather than inferred from an archive or a controller's claim to speak for them. Dependent participants and outsiders retain standing even when their protection does not maximize output. The physical inquiry can establish which integrations sustain differentiated capacities; the ethical commitment explains why their fate matters. That combination answers the strongest Unitology-style challenge conditionally, without pretending that a dissipation scalar proves agency or legitimacy.
