# Small forcing gives an actual local persistence theorem

The missing high-neighborhood and stationary-response assumptions in the weak-forcing calculation can be proved for sufficiently small forcing. The result below is an existence theorem for the supplied open, large-copy, source-free concentration model. It gives no numerical admissible forcing amplitude and makes no claim about the separate spatial experiment, finite molecular counts or a finite-temperature immigration process.

## Model and statement

Let a,b,delta,gamma,alpha,M and kappa be positive. In equations write \(\delta,\gamma,\alpha,\kappa\) for the corresponding parameters. Set

\[
S=R_++R_-,\quad D=M-S-L,\quad h(L)=e^{-\alpha L}.
\]

Consider the piecewise deterministic Markov process

\[
\dot R_\pm=R_\pm[aD-\delta h(L)(1\pm\varepsilon\xi_t)],
\qquad
\dot L=bSD-\gamma L,                              \tag{1}
\]

where the autonomous gate \(\xi_t\in\{-1,+1\}\) switches in either direction at rate \(\kappa>0\). Productive copying and L production are paid chemical channels; fuel/waste concentrations are supplied chemostats with their work counted. The equations are the negligible-reverse/source approximation specified in the preceding theory notes. There is no diversity objective, assigned group label or population controller.

Assume the unforced aggregate system

\[
\dot x=f(x),\qquad x=(S,L),\qquad
f(S,L)=(S[a(M-S-L)-\delta h(L)],\ bS(M-S-L)-\gamma L)
\]

has a strictly physical fixed point \(x_0=(S_0,L_0)\), with \(S_0,L_0,D_0=M-S_0-L_0>0\), whose Jacobian J is Hurwitz. Assume additionally \(D_0>S_0\).

**Theorem.** There is \(\varepsilon_*>0\) and a fixed compact aggregate neighborhood E of \(x_0\), strictly inside \(S,L,D>0\), such that for every \(0<\varepsilon<\varepsilon_*\):

1. E is forward invariant for every bounded gate history and every changing type composition in (1).
2. Each pure-type boundary has a unique stationary law in E. Its stationary response has a first-order expansion with a uniform quadratic remainder.
3. The two missing-type invasion exponents are equal and strictly positive.
4. Starting with both types positive and \((S,L)\in E\), the two types are stochastically persistent in the occupation-measure sense. In particular an interior invariant probability exists.

No uniqueness or mixing claim for the full two-type invariant law is required. E is a neighborhood in concentration state space, not a self-produced physical enclosure.

## 1. One Lyapunov matrix gives a common invariant high region

Write

\[
q(S,L)=(-\delta h(L)S,0)^T.
\]

For \(S>0\), the aggregate dynamics is exactly

\[
\dot x=f(x)+u_tq(x),\qquad
u_t=\varepsilon m_t\xi_t,\quad
m_t=(R_+-R_-)/S\in[-1,1].                         \tag{2}
\]

This representation requires no slow-frequency averaging. It remains valid even when the composition fluctuates rapidly.

Since J is Hurwitz, the matrix

\[
P=\int_0^\infty e^{J^Tt}e^{Jt}\,dt
\]

is positive definite and satisfies \(J^TP+PJ=-I\). Let \(p_{\min},p_{\max}>0\) be its smallest and largest eigenvalues and \(p=\|P\|=p_{\max}\). Use the chosen concentration and time units throughout.

Choose \(\rho>0\) so that the closed Euclidean ball of radius rho around \(x_0\) is strictly physical and

\[
\sup_{|x-x_0|\leq\rho}\|Df(x)-J\|\leq\frac1{8p}.
\]

Continuity supplies such a positive radius. On this ball let

\[
Q_0=\sup|q(x)|<\infty,\qquad Q_1=\sup\|Dq(x)\|<\infty.
\]

Choose \(r>0\) with \(r/\sqrt{p_{\min}}<\rho\), and define the convex ellipsoid

\[
E=\{x:(x-x_0)^TP(x-x_0)\leq r^2\}.
\]

Take an initial forcing bound \(\varepsilon_0>0\) sufficiently small that

\[
\varepsilon_0<1,\qquad
\varepsilon_0Q_1\leq\frac1{8p},\qquad
\varepsilon_0\leq\frac{r}{8pQ_0\sqrt{p_{\max}}}.      \tag{3}
\]

These are existential but computable inequalities, not an unprovided assumption about the gate process.

For \(z=x-x_0\) and \(V=z^TPz\), the mean Jacobian between x and \(x_0\) gives

\[
\dot V\leq-\tfrac12|z|^2+2p\varepsilon_0Q_0|z|.
\]

On \(\partial E\), \(|z|\geq r/\sqrt{p_{\max}}\), so (3) implies

\[
\dot V\leq-\tfrac14|z|^2<0.
\]

E is therefore forward invariant for every measurable input \(|u_t|\leq\varepsilon_0\), including every actual composition/gate input in (2). At gate jumps concentrations do not jump. The aggregate high-state region is thus protected against all bounded gate histories of sufficiently small amplitude, not merely against typical or short histories. Its compactness yields positive lower bounds \(S\geq s_{\min},L\geq\ell_{\min},D\geq d_{\min}>0\).

## 2. Uniform contraction gives unique resident stationary laws

On a pure boundary, \(m\equiv\sigma\in\{-1,+1\}\), so the aggregate resident field is

\[
F_{\sigma,\varepsilon}(x,\xi)=f(x)+\sigma\varepsilon\xi q(x).
\]

The first two bounds in (3) and the bound on Df give, throughout E,

\[
[D_xF_{\sigma,\varepsilon}]^TP
 +P[D_xF_{\sigma,\varepsilon}]\preceq-\tfrac12 I.   \tag{4}
\]

For two resident trajectories driven by the same gate history, convexity of E permits integration of these Jacobians along their joining segment. Consequently

\[
\frac d{dt}|x_t-y_t|_P^2
\leq-\frac1{2p_{\max}}|x_t-y_t|_P^2,
\quad
|x_t-y_t|_P\leq e^{-t/(4p_{\max})}|x_0-y_0|_P.      \tag{5}
\]

Construct the gate on the whole real line in its stationary, ergodic law. Starting a resident anywhere in E at time -T, (5) shows that its state at time t converges as \(T\to\infty\), independently of the starting state. The limit \(x^\varepsilon_\sigma(t)\) depends measurably only on the past gate history and is shift covariant. Its joint law with the gate is therefore stationary and ergodic.

It is the only resident stationary law in E. Indeed any two such laws have the same symmetric stationary gate marginal. Couple equal initial gate states and all subsequent gate jumps; (5) contracts the remaining concentration difference to zero. Stationarity then forces the two joint laws to coincide. This does not require a bracket/accessibility condition for the resident subsystem.

## 3. The response expansion has a uniform remainder

For the resident with sigma=+1 write \(z^\varepsilon=x^\varepsilon_+-x_0\). The same Lyapunov estimate, with Young's inequality, yields

\[
\dot V\leq-\frac{V}{4p_{\max}}
              +4p^2Q_0^2\varepsilon^2.
\]

Taking the stationary pullback limit proves a bound \(\|z^\varepsilon\|_\infty\leq C_z|\varepsilon|\), with a finite constant independent of the gate history.

Define the bounded stationary linear response

\[
y(t)=\int_0^\infty e^{Js}q(x_0)\xi_{t-s}\,ds.
\]

Write \(N(z)=f(x_0+z)-Jz\). Smoothness gives \(|N(z)|\leq C_N|z|^2\) on the chosen neighborhood. The difference \(w=z^\varepsilon-\varepsilon y\) satisfies

\[
\dot w=Jw+N(z^\varepsilon)
  +\varepsilon\xi_t[q(x_0+z^\varepsilon)-q(x_0)].
\]

The last two terms are uniformly \(O(\varepsilon^2)\). Since J is Hurwitz, their convolution against \(e^{Js}\) is uniformly bounded by \(C_w\varepsilon^2\). The bounded entire solution therefore obeys

\[
\boxed{x^\varepsilon_+(t)=x_0+\varepsilon y(t)
                            +O(\varepsilon^2)}                  \tag{6}
\]

uniformly over all gate histories and all t. This proves the regular expansion needed for the invasion calculation; it is not an appeal to unverified linear-response behavior. The negative resident has the corresponding sign-reversed first-order response. The model is an ordinary differential equation between finite-rate jumps, so no white-noise or Ito correction enters this calculation.

## 4. Certified positive leading coefficient implies positive exact invasion

The resident has \(R_+=S\geq s_{\min}\) and bounded S. Thus its long-run per-capita growth is exactly zero. The rare opposite type therefore has invasion exponent

\[
\lambda(-\mid+)
=2\varepsilon\delta\,\mathbb E[h(L^\varepsilon_+)\xi].
\]

With \(h_0=h(L_0)\), \(J_c=\delta h_0S_0=aS_0D_0\), and \(\mathbb E\xi=0\), equation (6) gives

\[
\lambda(-\mid+)=-2\alpha\delta h_0\varepsilon^2
                    \mathbb E[y_L\xi]+O(\varepsilon^3).
\]

Because \(\mathbb E[\xi_t\xi_{t-s}]=e^{-2\kappa s}\),

\[
\mathbb E[y\xi]=(2\kappa I-J)^{-1}(-J_c,0)^T,
\qquad
\mathbb E[y_L\xi]
=-\frac{J_c b(D_0-S_0)}
 {(2\kappa)^2-2\kappa\operatorname{tr}J+\det J}.
\]

Consequently

\[
\boxed{\lambda(-\mid+)=K\varepsilon^2+O(\varepsilon^3),\quad
K=\frac{2\alpha\delta h_0J_c b(D_0-S_0)}
 {(2\kappa)^2-2\kappa\operatorname{tr}J+\det J}>0.}               \tag{7}
\]

The error is genuinely bounded: the uniform Taylor remainder above gives a finite C such that its absolute value is at most \(C\varepsilon^3\) after possibly reducing \(\varepsilon_0\). Taking \(0<\varepsilon<\min(\varepsilon_0,K/(2C))\) makes the exact exponent positive. If C vanishes, no additional restriction is needed. Gate symmetry gives the identical exponent \(\lambda(+\mid-)\). Kappa is fixed and positive; no uniform conclusion in a singular zero-switching limit is asserted.

For \(a=b=\delta=\gamma=1,\alpha=4,M=.95,\kappa=.25\), the independently outward-enclosed high state is Hurwitz and satisfies \(D_0>S_0\). [certify_diversity_coefficient.py](certify_diversity_coefficient.py) and [diversity-coefficient-certificate.json](diversity-coefficient-certificate.json) certify

\[
0.189215105435<K<0.189215105460.
\]

They establish the strict sign required by (7). This proof supplies the previously missing existence of a valid sufficiently small amplitude, but does not calculate a numerical value for that amplitude.

## 5. Persistence follows without averaging the type frequencies

Let

\[
\mathcal K=\{(R_+,R_-,L,\xi):R_\pm\geq0,
           (R_++R_-,L)\in E,\ \xi\in\{-1,+1\}\},
\]

and let \(\mathcal K_0\) be its two pure-type facets \(R_+R_-=0\). Section 1 makes this a compact forward-invariant state space. The all-extinct state is excluded by \(S\geq s_{\min}\); no trajectory inside this high region can approach the separate low basin. The smooth bounded vector fields and finite positive gate rate give a Feller process with invariant extinction facets and positive types remaining positive at every finite time.

There are exactly two ergodic stationary laws on \(\mathcal K_0\): the two unique resident laws established in Section 2. Under either law the resident's mean per-capita growth is zero and the absent type's mean growth is strictly positive by Section 4. Thus a positive weighted sum of species invasion rates is positive for every ergodic boundary law, with both weights equal to one.

Apply [Benaïm, Stochastic Persistence, §3.3 and Theorems 6.1 and 4.4](https://arxiv.org/pdf/1806.08450). Auxiliary feedback coordinates such as L are allowed. Concretely use

\[
V=-\log(R_+/M)-\log(R_-/M),\qquad
\mathcal LV=-(g_++g_-),
\]

where \(g_\pm=aD-\delta h(1\pm\varepsilon\xi)\). V is nonnegative on the interior, diverges at either extinction facet, and its generator extends continuously to \(\mathcal K\). Gate jumps leave V unchanged; its generator martingale is identically zero. On each boundary stationary law the mean generator is negative. Compactness supplies tightness. These facts verify the average-Lyapunov persistence criterion rather than simply naming it.

For every \(\eta>0\) there is a compact \(\mathcal K_\eta\subset\mathcal K\setminus\mathcal K_0\) such that, from any initially positive pair in \(\mathcal K\),

\[
\mathbb P\left(\liminf_{T\to\infty}
\frac1T\int_0^T\mathbf1_{\mathcal K_\eta}(X_t)dt
\geq1-\eta\right)=1.
\]

All occupation-measure weak limits are invariant probabilities supported on the interior. No approximation to the actual frequency dynamics was used in this step. The joint two-type process need not itself be contracting or have a unique mixing law.

## What has and has not been established

Within this declared concentration model, this is actual small-nonzero-forcing stochastic persistence of two physically different response classes in an actively renewing high organization. The high neighborhood has positive lower bounds on S,D,L. Paid productive rates \(aSD\) and \(bSD\), and L turnover \(\gamma L\), remain positive: the result is not preservation by an immortal inherited shield. Continuing chemical work is necessary. No external diversity reward or maintained population proportion appears.

The theorem conditions on both types being present in the high region. Paid copying error is a separate physical route into that interior; this proof does not implement a mutation-rate switch after a first birth or prove the full fixed-mutation model. It does not give a founder establishment probability. Source-free concentrations omit demographic loss, and exact finite-temperature births destroy the invariant zero-type facets used in this theorem. Finite-eta lineage retention requires a different statement, not a silent extension of this proof. Nor does the result prove spatial localization, collective inheritance, continuous-trait branching, increased available free energy, a biodiversity optimum or participant agency. Those limits do not undo the narrower existence result: reciprocal renewal can, under explicitly supplied chemistry and weak fluctuating responses, generate sustained functional retention without a diversity reward.
