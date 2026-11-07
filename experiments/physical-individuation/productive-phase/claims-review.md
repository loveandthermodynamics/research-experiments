# Independent claims review

21 October 2026. **Pass for the reported scope; both prospective scientific claims fail.** This review checks the frozen decision rules, all saved seed-level analysis rows, the analytical constraints and the current `checkpoint.md`, `consolidation.md` and `RUN.txt`. It does not change a threshold, source, trajectory or prior result. No additional physical simulation was run.

## Decisions and evidence

The complete design is eight independent seed rows at each of two sizes and four arms: 64 histories. Each window contains 21 snapshots, averaged within seed before inference. Recomputing the saved window means and all 22 gate means, standard errors, Student-t7 intervals and decisions agrees exactly with `analysis.json`. All gates are estimable. Snapshots and events have not been treated as additional independent observations.

| Prospective claim | Side 16 | Side 32 | Decision |
|---|---:|---:|---|
| A: B mass in connected clusters of at least nine sites | 0.20288 [0.12929, 0.27646] | 0.25852 [0.23222, 0.28481] | Both required lower bounds >0.5 fail |
| B: specialist minus full-generalist template density near those clusters | −0.11782 [−0.14160, −0.09403] | −0.09558 [−0.10302, −0.08814] | Both required lower bounds >0 fail |

Intervals are two-sided 95% t7 intervals over the declared seed rows. Claim A passes 18 of its 20 gates but fails its conjunction. Its successful wavelength, abundance, current and size screens cannot substitute for the failed dense-mass screens. Claim B is contrary to the proposed specialist advantage, rather than merely imprecise: all eight paired differences are negative at each size, and its declared bootstrap sensitivity retains that direction. Neither cluster cutoff nor window should be changed to obtain a pass. The study does not establish finite maintained phase units in either architecture; no separate generalist phase claim was preregistered.

The positive result is narrower and substantial. Specialist late template densities are 0.5702/0.5734. Its whole-horizon exact time-average densities are 0.4745/0.4834, with positive lower confidence bounds for the contrast against the one-activity finite-horizon **expectation** bound of 0.05385. Both late productive currents have positive lower bounds at both sizes. These are simulation estimates compared with a mathematical expectation bound, not pathwise population ceilings or an exact determination of an unknown ensemble mean.

Original X instances clear in all 16 specialist histories. Paid copying by their descendants continues after clearance, and founder descendants remain at the horizon in all 16. The first qualifying paid Y appearance by time 20 occurs in 15/16, its copying subtree reaches four live instances in 11/16, and survives to the horizon in 9/16. The subtree includes both later phenotypes; thermal arrivals do not reconstitute extinct ancestry. These counts demonstrate paid genealogical succession, not causal mutant benefit, phenotype continuity, organizational memory, aging or the persistence of an individual.

## What the counterexamples now exclude

The exact balances and the physical comparisons agree. With no background A/B channel, one activity alone has work equal to a state gradient and an effective equilibrium Hamiltonian. Its high B aggregation does not establish a stationary paid renewal current. Both activities permit a physical closed cycle of affinity 2mu/theta = 24; irreducibility and positive paired rates imply positive stationary probability-flux entropy production. That theorem guarantees neither a large current nor a dense unit.

Complementary transformations do not require separate taxa. The exact one-genotype/two-inherited-phenotype relabeling remains observationally equivalent. The homogeneous half-strength generalist example already defeats necessity from reciprocal roles alone. The current unpenalized spatial generalist is a stronger comparative counterexample in this regime: near-cluster template density and total population are higher, while chemical work is also higher. Work per area is 72.6504/71.2852 for generalists versus 51.4922/53.0449 for specialists. This is neither a matched-function efficiency comparison nor a universal dominance theorem. It does rule out claiming that the supplied grammar has demonstrated a specialization advantage. Target-substrate exclusivity also means the existing per-parent hazard cap cannot justify imposing an arbitrary half-rate penalty on a promiscuous catalyst.

No new result closes the remaining obligations in `full-goal-obligations.md`. In particular, aggregate paid renewal is not local reciprocal integrity; B residence is not informational aging; ancestry is not functional memory; a finite alphabet is not open-ended biodiversity; and no additional level of individuation has been generated. That document is explicitly a pre-outcome snapshot and should remain unchanged as provenance. Thermodynamic performance cannot certify participant interests or moral legitimacy. Voluntary fusion, dependence and sacrifice remain compatible with the participant criterion; a separation requirement would not answer Unitology or Resident Evil. Hestia's social-prediction and moral-control commitments remain outside the chemical measurements.

## Next scientific question

The useful next step is an analytical physical constraint, not another outcome scan: **what realizable source–sink, transport or catalytic-capacity constraint can sustain dense local renewal, and can the best admissible full or adaptive generalist satisfy it?** Derive that constraint from explicit molecular states and the actual current/material balances at fixed useful renewal and resources. If an operating incompatibility or switching cost is needed, specify its physical mechanism and count its construction and repair cost; do not infer it from a maximum hazard or choose it to reverse the ranking. Any approximation must retain bypasses and boundary transport. A linear instability calculation alone would not exclude nucleation or prove maintained finite units. Only after a maintained candidate exists would one matched local interruption test its causal integrity. This keeps the failed phase mechanism and the generalist counterexample in the same question without treating aging or higher heredity as already solved.

## Verification and provenance

The separate independent physical audit passes all 64 histories, 11,551,005 accepted events, 100,263,216 proposals including rejections, and 5,184 samples, with maximum running first-law discrepancy 2.26e−9. That audit supplies complete trajectory replay; this claims review independently checks saved summaries and inference rather than claiming another raw replay. The checkpoint accurately distinguishes counted chemical work from excluded initial preparation and chemostat construction. `RUN.txt` matches the current CLI and correctly requires a disposable extraction for commands that regenerate outputs. The compact kit requires regeneration before absent raw histories can be replayed.

The prospective protocol binds `balance-and-cost.md` and `generalist-comparison.md`, but not their certificate scripts or JSON reports. The final certificate rerun and manifest additionally bind those artifacts; they cannot retroactively make them protocol-bound. This distinction is stated correctly in the packaging drafts. The review does not certify an archive that had not yet been finalized at review time.

Reviewed input SHA256 digests:

```text
protocol.json                 77809c13b9704fccf5976a955763fa5a4fdbe7dfc712e2befc303db68441f31f
study-design.md               af61c723f94c0d23b0c12ba46ee051ef3423bd3cbd274726888f4f4b84d76321
analysis.json                 f9f7e82dbf809f3daef9b94f2f6b2f75f9f164b2d86c32c07581137d1cd1c3da
full-goal-obligations.md       472ea4d34ca50b1c930b9514647f287d77995312078d12750c9847236204a253
checkpoint.md                 7bba225261518f4f4d7a5e0f20c21e21cd63c6a4993396d38edd971167f9f0e4
consolidation.md              ca73a4750c003cebe1aa810fa949070817f181ad51f6dd5d8cd71ad2570d2675
RUN.txt                       bc16ed0fb8e5a87704761a287f756586a5a69e748bc1bba3751e0bd6f4470ae0
independent-verification.json 32be097407b66c5aa94aa3efd624de2439b019328f458a6c16c9c03e19330b5f
certificate-reverification.json 6df6fadb0c49f4232d590787cc21feea7ac0f9e3a63da8ed4b757a2637f57e54
```

No mandatory claims, numerical or CLI correction remains in the reviewed drafts. The full research goal remains active and unproved.
