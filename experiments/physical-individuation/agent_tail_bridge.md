# Reversible gates, heavy recovery tails, and paid first repair

This is a conditional bridge between the earlier best-survivor tail calculation and finite-resource chemistry. It proves a latency benefit. It does not establish an increase in integrated chemical availability or a biodiversity optimum.

## A finite gate with an explicit equilibrium

A gate has a healthy state H and bad states B_j, j=1,...,m. Choose positive weights w_j summing to one, recovery rates r_j>0, an entry rate lambda>0, and thermal energy theta>0:

\[
H\mathrel{\mathop{\rightleftharpoons}^{\lambda w_j}_{r_j}}B_j.
\]

Set E_H=0 and

\[
E_j=\theta\log\frac{r_j}{\lambda w_j},\qquad
\pi_H=\frac1{1+\lambda\sum_jw_j/r_j},\qquad
\pi_j=\pi_H\lambda w_j/r_j.
\]

Then pi_H lambda w_j = pi_j r_j: every gate edge satisfies detailed balance. These rates can also be represented by ordinary activation barriers. With an attempt rate a at least max_j(r_j,lambda w_j), choose barrier V_j=theta log[a/(lambda w_j)]. Both barrier heights above their adjacent wells are nonnegative and reproduce the two rates.

On entry to a bad episode, the recovery survival is

\[
S_{\rm entry}(t)=\sum_jw_je^{-r_jt}.
\]

The ideal Lomax survival (1+t/tau)^(-alpha) is exactly a mixture of exponential recoveries with rate R distributed Gamma(shape alpha, rate tau). Positive finite quadrature weights approximate that mixture. For a rate grid covering [r_min,r_max], tail probability omitted or clipped eta, and within-bin rate error at most Delta, the absolute survival error on [0,T] is at most eta+T Delta. A finer quadrature gives a sharper numerical bound. Relative tail claims require this error to be small compared with the tail probability being measured.

Every finite mixture ultimately has an exponential tail. More generally, a finite-state continuous-time chain has phase-type first-passage times; an almost-surely reached target has an eventual exponential times polynomial tail. Reversibility makes the killed generator similar to a symmetric matrix, so its survival is a finite sum of exponential modes. A finite construction therefore supports a declared approximation window, never a literal asymptotic Pareto law.

## Stationary preparation changes the exponent

An infinite gate of this form has a normalizable stationary distribution only if E[T_entry]=sum_j w_j/r_j is finite. For the ideal Lomax mixture this requires alpha>1, with E[T_entry]=tau/(alpha-1). Consequently an infinite-mean entry duration with positive healthy-to-bad flux cannot simultaneously be the positive-recurrent equilibrium of this gate.

At a stationary observation time, conditional on being bad, the sampled weights are proportional to w_j/r_j. Long episodes are sampled more often. The ideal residual recovery survival is therefore

\[
S_{\rm residual}(t)
=\frac{E[R^{-1}e^{-Rt}]}{E[R^{-1}]}
=(1+t/\tau)^{-(\alpha-1)}.
\]

Thus alpha in (1,2] permits a normalizable equilibrium gate but an infinite mean residual waiting time. An entry-conditioned preparation has exponent alpha; a stationary bad-state preparation has exponent alpha-1. These preparations must not be interchanged.

## Matched finite fuel and a paid-first-activation theorem

Use the same n physical gates, their same joint initial law, n labeled material sites, K finite fuel molecules, well energies, and thermal bath in both comparisons. In the common-response case every site reads gate 1; in the independent-response case site i reads gate i. All unused gates remain in the common-response inventory. Response changes only reaction barriers, not well energy or initial energy.

This is a separate activation model. It adds a physically retained inactive template or basal catalyst that can initiate its own activation; the helper-catalyzed repair and copying channels in the current prototypes cannot initiate from an all-inactive, catalyst-free state. It therefore does not establish rescue from irreversible extinction.

Each site has that retained template and inactive/active state a_i=0/1. All sites start inactive, with F=K and W=0. Healthy-gated activation is the reaction D_i+F+T_i <-> A_i+W+T_i, where the retained catalyst T_i is unchanged. Its forward rate is rho f when a_i=0 and its reverse rate is

\[
\rho(K-f)\exp[-(e-\epsilon)/\theta]
\]

when a_i=1. The reaction consumes one F and makes one W, raising site energy by epsilon; e>epsilon is the fuel energy gap. The gate-H indicator multiplies both directions. In a bad gate both directions of this channel are absent; a small paired leakage mobility would impose another cutoff.

Include passive damage a_i=1 -> 0 at rate delta and its paired thermal activation at eta=delta exp(-epsilon/theta). The complete finite model is reversible with the same reference in both assignments:

\[
\pi(a,f,g)\propto
\frac{\exp\{-[\epsilon\sum_i a_i+ef+\sum_\ell E_{g_\ell}]/\theta\}}
{f!(K-f)!}.
\]

The activation edge ratio includes f/(W+1), exactly matching the factorial mixing term. Gate transitions and passive transitions also obey detailed balance. No externally scheduled gate or free fuel restoration is needed.

Before the first activation, f=K and all a_i=0. Hence the paid hazard per healthy site is the constant nu=rho K. Let A_t be healthy occupation time of one gate and Z_t=exp(-nu A_t). For independent initial gates with the same distribution,

\[
P(U_{\rm independent}>t)
=e^{-n\eta t}(E[Z_t])^n,
\]
\[
P(U_{\rm common}>t)
=e^{-n\eta t}E[Z_t^n]
\ \ge\ e^{-n\eta t}(E[Z_t])^n.
\]

The inequality is Jensen's inequality. Independent responses therefore shorten first-activation latency in stochastic order, with the same finite shared fuel and material. A paid first activation consumes exactly one fuel molecule. The thermal factor explicitly distinguishes rare bath activation from that paid channel; the probability of any thermal activation before time T is at most n eta T. When K=0 only the identical thermal term remains.

For the untruncated ideal stationary gate with eta=0, let S_nu(t)=E[Z_t]. Then beta=alpha-1 and

\[
S_\nu(t)\sim\pi_B(\tau/t)^\beta,
\quad
P(U_{\rm independent}>t)\sim\pi_B^n(\tau/t)^{n\beta},
\quad
P(U_{\rm common}>t)\sim\pi_B(\tau/t)^\beta.
\]

Here pi_B=1-pi_H and nu>0. To see the first equivalence, a stationary bad gate must first finish its residual episode. All subsequent healthy waits and repeated entry episodes have a moment of some order p strictly between beta and alpha; their tail is negligible beside the residual t^(-beta) tail. This also supplies an elementary proof without assuming instantaneous repair. An entry-conditioned preparation instead has the earlier exponent alpha, and independent responses multiply that exponent by n.

In the ideal stationary model with n independent response gates, the p-th first-activation moment is finite exactly when p<n(alpha-1); equality is divergent. The common-response case remains finite only for p<alpha-1, despite its n sites. For alpha=1.5, three independent responses make mean first activation finite while one response, or two at the boundary, does not. This is a statement about latency moments. Finite gates, thermal leakage, and finite observation windows remove the literal divergence.

## Numerical finite-mixture check

Computed a 193-state single-gate killed generator: 192 logarithmic recovery bins over [10^-9,100], Gamma shape alpha=1.5 and rate tau=1, exact Gamma bin masses, geometric bin representatives, lambda=.2 and nu=1. Gamma probability below 10^-9 is folded into the first representative; probability above 100 is folded into the last. The lower mass is only 2.37883e-14, but stationary length bias divides it by a rate near 10^-9: its contribution to E[T_entry] is 2.22699e-5. Dropping and renormalizing the tails instead gives 2.00137906, a distinct valid finite approximation. [agent_tail_numerics.py](agent_tail_numerics.py) reproduces the folded-tail convention used here.

The initial gate distribution was its stationary distribution in every comparison. The gate has pi_H=.71414275, E[T_entry]=2.00140133, and mean gate energy .07139768 at theta=.05. Every gate transition's detailed-balance log residual was below 7.2e-15. Both full comparisons contain three such gate units and thus the same expected gate energy .21419303.

| t | Single first-activation survival S_nu | Three independent responses S_nu^3 | Three sites, common gate S_3nu |
|---:|---:|---:|---:|
| 100 | .02886682 | .0000240545 | .02859472 |
| 1,000 | .00904913 | .000000741005 | .00904067 |
| 10,000 | .00285446 | .0000000232581 | .00285420 |

The apparent exponents from 1,000 to 10,000 are .5011 for one response, 1.5032 for three independent responses, and .5007 for three sites sharing a gate. The finite residual recovery mixture differs from the ideal (1+t)^(-.5) by less than .153 percent at these three times. These are exact finite-generator matrix-exponential calculations, not Monte Carlo and not a simulation of post-activation regeneration. Its eventual smallest decay rate is approximately 1.06818e-9, confirming the eventual exponential cutoff.

For example K=12, rho=1/12, epsilon=.8, e=1.2, theta=.05, and delta=.03 realize these paid hazards. Then eta=3.3761e-9; for three sites the exact passive-activation factor at t=10,000 is .99989872. It has negligible effect in the displayed window but forbids extrapolating the power law indefinitely.

## Where the inference stops

- After the first activation, shared fuel f is random and depleted. Repair hazards no longer factor into independent gate occupation times. Survivor-catalyzed repair, contacts, competition, recolonization, or synchronized repair add further dependencies. The theorem concerns the first useful repairer, not completion or repeated regeneration.
- A retained template was supplied. If the only repairer or template decays during a long outage, the recovery-tail calculation does not imply its survival. In particular, a repair requiring an already active helper cannot restart from total inactivity without a separate physical activation pathway. Adding an independent loss hazard cuts off useful survival even when gate recovery remains slow.
- A common gate, contagion, or shared adverse exposure can restore the single-response exponent. If independent waiting X has survival proportional to t^(-beta), an additional common waiting C with exponent gamma gives U=max(C,min_i X_i) exponent min(gamma,n beta). Independence of response labels alone is insufficient.
- Finite material bounds occupied response classes and the number of retained templates. Splitting stock can create founder extinction, insufficient catalytic concentration, and repair costs. The best-survivor exponent omits these costs; increasing n is not automatically feasible or beneficial.
- First-activation latency is not integrated active availability. First activation can restore only a small fraction of the stock; repairing the rest costs further fuel and time. No exergy advantage follows from the latency inequality alone. In the complete finite model, active-subsystem availability is bounded by n theta log(1+exp(epsilon/theta)); its integral on [0,T] is bounded by T times this quantity. Closed-system total availability cannot be manufactured by the gates. With positive passive reversals the finite chain ultimately equilibrates, and its availability decays exponentially.

The supported bridge is: differentiated exposure can reduce long waits for a fuel-funded restart under matched resources, and stationary long-outage sampling can make this benefit especially large. Demonstrating that the earlier contact, mutation, and paid-regeneration models retain this benefit under these gates requires their joint generator and an actual availability calculation. Recursive individuation and a universal diversity optimum are not consequences of this theorem.
