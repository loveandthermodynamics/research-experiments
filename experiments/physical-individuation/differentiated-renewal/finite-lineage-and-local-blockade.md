# Finite lineages and a local renewal intervention

This note supplies finite-time criteria for the actual finite, driven spatial CTMC, with all reverse channels, fixed copying mutation and autonomous gates retained. It does not replace that process by an infinite-population invasion exponent. The results below follow directly from its finite generator and its recorded parent/target channels; no simulation is reported here.

## 1. Presence is not establishment

A reversible two-state site already supplies an exact warning. Let D become R1 at rate u and R1 become D at rate d, with \(u/d=e^{-\epsilon/\theta}>0\). Its stationary R1 probability is \(u/(u+d)>0\), although R1 has no copying ability. In the full model thermal births and R0-to-R1 copying errors can likewise maintain R1 occupancy without any particular R1 founder establishing a lineage. A stationary count, repeated appearance or a positive immigration-supported mean is therefore insufficient.

For a declared paid founder, annotate copying ancestry offline. A copy inherits the ancestry of the channel's recorded parent/helper. Mutation and back-mutation change genotype while preserving this ancestry. Material transport moves it. Passive loss and reverse copying remove the target instance and its ancestry. Unrelated thermal births and copies from outside the lineage do not acquire it; neither does recycled D merely because it contains the same material. Record L products attributed to a catalyst separately from R copying descent. These are readout annotations and never enter a rate.

The existing ledger distinguishes `instance_created`, `instance_removed`, `parent_instance`, `catalyst_instance` and helper/target sites. Reverse helper attribution can be reconstructed from the pre-event state even when the stored reverse-event catalyst field is zero. This is enough to reconstruct ancestry through the existing channels. If the final marked R instance disappears, its copying lineage is extinct even though chemically identical R1 can later reappear elsewhere. No chemical irreversibility or zero-temperature approximation is needed for this genealogical extinction.

Track both total lineage size \(N_\ell=N_{\ell,0}+N_{\ell,1}\) and its genotype-resolved counts. Growth of a lineage that has lost R1 is not evidence of continuing R1 function. A bounded establishment event can require a declared descendant threshold, R1 descendants at a declared later time, and turnover of the original founder. The thresholds, horizon and eligible founders must be specified before examining their later success. Retain founders that fail or appear late; do not select the largest surviving clan after the fact.

## 2. An exact lineage drift and an extinction bound

Augment the physical state only by the neutral ancestry colors. Let Q be its generator. All present chemical channels change global lineage R count by +1, -1 or zero; gate, reservoir and transport jumps do not change that global count. Define from the **actual propensities**

\[
B(x)=\sum_{e:\Delta N_\ell=1}q_e(x),\qquad
D(x)=\sum_{e:\Delta N_\ell=-1}q_e(x).
\]

B includes marked-parent copies of either genotype; D includes every marked-target passive or reverse productive loss. Fixed mutation, surrounding unmarked templates, finite fuel slots and protection remain in these rates. There is no external lineage immigration. Direct generator calculation gives

\[
QN_\ell=B-D,\qquad
N_\ell(t)-N_\ell(0)-\int_0^t(B_s-D_s)ds=M_t,
\]

with martingale quadratic variation compensator \(\int_0^t(B_s+D_s)ds\). Thus a realized birth count is not itself a per-capita growth coefficient. The full hazard and all failed lineages matter.

Here is a conservative, rigorous null. Suppose throughout a predeclared reachable domain, before the assay stops,

\[
B\leq bN_\ell,\qquad D\geq dN_\ell,\qquad 0<b<d.
\]

Set \(r=d/b>1\). Then

\[
Q r^{N_\ell}
=r^{N_\ell}(r-1)(B-D/r)\leq0.
\]

Stopping at extinction, first size m, domain exit or horizon T yields

\[
\boxed{\Pr(\text{reach m before extinction, exit and T})
\leq\frac{r^{N_\ell(0)}-1}{r^m-1}.}                \tag{1}
\]

For b=d the corresponding bound is \(N_\ell(0)/m\). If the subcritical rate bounds hold globally, also \(\Pr(N_\ell(t)>0)\leq N_\ell(0)e^{-(d-b)t}\). Finite population and finite horizon make the stopped-martingale argument bounded. A stronger establishment event contained in reaching m inherits its upper bound. Verified rate envelopes, not averages measured only on successful trajectories, are required. If an envelope is valid only in a domain, count its exit separately; do not silently extend the inequality afterward.

The exact finite-time alternative is the hitting probability for a specified establishment set A versus lineage extinction B. On the other augmented states it satisfies \(\partial_t h_t=Qh_t\), with boundary values 1 on A, 0 on B and initial value zero outside A. Its complete-state calculation may be large, but the same probability can be estimated from complete CTMC trajectories with all outcomes retained. Neither method substitutes stationary thermal occupancy for establishment.

## 3. A neutral ancestry check and a response counterfactual

An exact tag-blind bookkeeping null is available even when locations and genotypes differ. At a declared preparation time take a cohort C of k existing R instances. Independently choose one uniformly to mark. Conditional on the entire later physical event/parent history,

\[
\mathbb E_{J\in C}[N_J(t)\mid\text{history}]
=N_C(t)/k,
\]

where \(N_C\) counts descendants of any cohort member and excludes later unrelated sources. Every descendant has one recorded copying ancestor, so the identity is immediate. It detects ancestry/readout bias without assuming equal local reproductive opportunities. It is not itself a causal test of R1 advantage. Nor should one assume lineage fraction \(N_\ell/N_R\) is a martingale: position, genotype, reverse channels and newly entering founders can all matter. Common response alone does not guarantee aggregate lumpability at finite reverse rates.

For a phenotype-specific test, branch the same physical snapshot into differentiated-response and common-response generators. Keep energies, mutation, gate state and physical gate process fixed, changing only the declared paired passive mobility map for the chemical genotype. The same marked founder is a readout in both branches. Compare its finite-time establishment probability and genotype-resolved descendants, retaining every branch. This measures the total causal effect of the response map in that prepared context; it need not be an intrinsic, context-free fitness advantage. Never change rates according to the founder's ancestry mark.

Paid L production attributed to a lineage should be counted in both directions: forward production minus reverse consumption, using the actual helper. A positive net current establishes material production, not benefit. A lineage can establish by exploiting somebody else's L pool. Contribution therefore needs an additional physical intervention.

## 4. A reversible, tag-blind local-production blockade

Choose a candidate geometric patch P from the preparation state, freeze its coordinate mask, and predeclare the observer outcome and horizon. A patch selected using past morphology is an assay preparation, not a group recognized by the autonomous generator. Let \(\mathcal E_P\) be a fixed set of directed helper-target locations for paid L production, such as edges whose helper j lies in P. Apply a factor \(c\in[0,1]\) to **both directions** of every selected productive reaction pair.

In general the requirement is

\[
c_e(x)=c_{\bar e}(x^e).
\]

For c>0 this preserves the paired ratio

\[
\log\frac{q_e^c(x)}{q_{\bar e}^c(x^e)}
=\frac{W_e-\Delta H_e}{\theta}.
\]

At c=0 both edges are removed; the ratio is then undefined, but remaining channels still obey local detailed balance. A selector based on the unchanged helper location/genotype is safe for a target D/F-to-L/W reaction. A selector based on whether that target currently is L need not be edge symmetric. A moving cluster/centroid rule must be independently checked on every reverse pair. An ancestry-dependent producer selector would introduce observer feedback and is excluded.

The reaction energy/work ledger is unchanged per enabled jump. Multiplying paired rates is a barrier intervention, not a free energy source. This separation between rate ratios and kinetic coefficients is consistent with [Rao and Esposito, Physical Review X 6, 041064 (2016), equation 50 and the chemical-work balance in equation 78](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.6.041064). The finite-generator ratio above is also directly checkable edge by edge.

Compare near-patch blockade with remote and sham preparations from the same physical snapshot. Remote selections should match lattice degree, shape/edge count, blockade strength and pre-treatment substrate, fuel, producer and gate opportunities as closely as possible. Sham leaves c=1. **Equal degree alone matches neither blocked propensity nor lost productive current.** Report initial forward/reverse hazards, integrated hazard exposures along each arm and actual chemical work. An integrated removed-channel propensity along a treated path is an exposure, not an observed number of counterfactually prevented events. Do not choose a remote control to match a favorable post-treatment flux: that would condition on an intervention consequence. If suitable matches do not exist, report the mismatch.

Under the model's stated barrier-only Hamiltonian idealization, switching c does no direct state-energy work because the well energies do not change. That is not a proof that a physical blocker, measurement or switching apparatus is free. Record the intervention/preparation protocol; include apparatus work if making a full physical efficiency claim. Otherwise restrict the conclusion to a kinetic causal counterfactual of the modeled system. Continue counting actual chemical supply, fuel/waste exchange, heat and interaction energy in every arm. Matching initial H is not a construction-work certificate.

## 5. The causal claim and its exact generator test

For a declared bounded outcome Y, such as integrated local R or regeneration after founder turnover, define

\[
\Delta_{\rm near}=\mathbb E[Y_{\rm sham}-Y_{\rm near}],\qquad
\Delta_{\rm specificity}=\mathbb E[Y_{\rm remote}-Y_{\rm near}].
\]

Positive values support causal local production dependence over this horizon and preparation protocol. Separate local disappearance from transport by tracking the same descendants globally and their boundary crossings. A barrier intervention can redistribute templates without destroying their renewal.

For a finite-state terminal readout \(Y=f(X_H)\), write \(Q_\lambda=Q_0-\lambda B_P\), where \(B_P\) is the paired-channel generator being attenuated and \(0\leq\lambda\leq1\). With identical initial row distribution p, matrix differentiation gives the exact identity

\[
\frac{d}{d\lambda}\mathbb E_\lambda f(X_H)
=-\int_0^H p e^{tQ_\lambda}B_P e^{(H-t)Q_\lambda}f\,dt.       \tag{2}
\]

For an occupancy integral, integrate this expression over its terminal time. The sign is not fixed by detailed balance: production can protect, waste scarce fuel or displace copying substrate. The claimed positive return must be measured or bounded, not assumed from the word protection.

A local blockade affecting all nearby producers attributes dependence to that physical production region, not automatically to one variant or lineage. To join differentiation to renewal, repeat the same intervention contrast under differentiated and common response with the same mutation and gate protocol. A positive interaction

\[
[\mathbb E Y_{\rm diff,remote}-\mathbb E Y_{\rm diff,near}]
-[\mathbb E Y_{\rm common,remote}-\mathbb E Y_{\rm common,near}]>0
\]

supports an increased role of local production under differentiated response. Also report the absolute differentiated-versus-common outcome: dependence alone need not improve maintenance. Shared supply protocols do not imply equal realized work; a higher-activity result using more work is not automatically an efficiency result.

The conservative joint criterion is therefore paid appearance, immigration-free descendant establishment, physically different response, actual paid renewal current, and a favorable matched local causal return. Any can fail independently. Passing them would establish differentiated renewal within the supplied spatial chemistry. It would not by itself establish a self-bounded individual, unlimited biodiversity, available-free-energy gain or moral participant status. The Unitology distinction remains the independently consequential interests and capacities of participants, not a lineage's success in these assays.
