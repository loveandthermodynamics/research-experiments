# Independent support-readout and design review

The current `cohort_support_readouts.py` agrees with an independent reconstruction on eight manually specified synthetic histories. `verify_support_readouts.py` reconstructs per-founder ancestry roots, physical occupancy, original-instance death times and exact event-to-event interval summaries; it checks global and regional balances independently of the production readout. No model rates are called. These are adversarial observer fixtures, not natural chemical simulations or a thermodynamic validation.

The fixtures establish the following software distinctions:

- A cohort can leave the frozen region, lose its original founder and continue globally through a copied mutant. Local disappearance is not lineage extinction.
- Thermal replacement after the original cohort dies increases R occupancy without restoring ancestry.
- A helper can move out of the original square and continue making L as a cohort descendant. Cohort contribution and fixed-site contribution differ.
- Reverse L production is attributed using the helper's **pre-event** tag even when its event record contains `catalyst_instance=0`. A fixture has cohort forward/reverse/net counts 2/1/1, while its near-site forward/reverse/net counts are 0/1/-1.
- Two cohort instances can cross the shell in opposite directions during a swap. Gross crossings are retained rather than replaced by net displacement.
- Backmutation preserves copying ancestry. A new child born outside the region is a birth there, not an outward material crossing. Copies after original-instance clearance are counted separately.
- The periodic nearest-neighbor shell is respected.
- An empty initial cohort has no turnover event. Its presence/cleared-after-presence flags are false and its post-clearance integrals are zero despite subsequent immigrant R occupancy.
- Reverse copying can eliminate the cohort while later unrelated thermal and paid copies remain present.

All reconstructed global cohort, local cohort, global R, local R and local L balance errors are zero in the fixtures. The return values match the readout function, including final genotype counts, clearance time, exact integrals, forward/reverse L currents and crossings. The source uses the correct pre-event helper tags for reverse production. No readout correction is requested at this checkpoint.

The independent public API is `reconstruct_support(run, near_sites)`, returning `measurements` and `diagnostics`; `verify_run(run, near_sites)` additionally compares the production readout. Importing this module performs no I/O and does not import the readout under test. The command-line entry can verify complete JSON or gzip raw histories after they exist. The current `support-readout-verification.json` records eight synthetic fixtures and **zero physical raw runs**; it must not be cited as verification of an ensemble. Model feasibility, rates, work, heat and paired statistics belong to the separate causal-study verifier.

## Prospective interpretation of the assay

The prospective design appropriately freezes selection at preparation, retains empty/unmatched preparations, uses immutable helper-site factors on paired reaction directions, shares the autonomous gate stream, and keeps ancestry out of rates. The regional readout is fixed coordinates plus a nearest-neighbor shell. Its cohort annotation distinguishes material turnover and transport; it supplies no molecular ownership or autonomous group boundary.

Matching `q_near H_near = q_remote H_remote` matches initial **total removed propensity**. When the raw hazards differ, the removed fractions and retained mobilities differ. It therefore does not isolate location while holding blockade strength fixed. It also does not match forward production, reverse fuel-regenerating consumption, net current, future hazard, chemical work or apparatus cost. There is no established monotonicity that makes the unequal strengths automatically conservative. The legitimate contrast is a **fixed-location effect under the specified initial-hazard-matched, strength-varying intervention protocol**.

Disjoint helper squares and a center-distance criterion do not imply disjoint support pathways. Their possible L-target neighborhoods can overlap, and a remote helper can lie in the near square's expanded readout shell. The current metadata records target-neighborhood overlap; final reporting should preserve this interference qualification. Shared material and L transport can also mediate remote effects. A favorable remote-minus-near contrast would not prove that support belongs exclusively to the near square.

The two predeclared positive primary intervals would support dependence of near-region R occupancy on continuing production at the treated nearby locations over this horizon. They would not independently establish differentiated-lineage benefit, an isolated equal-strength location effect, a comoving individual, collective inheritance, memory, aging or subjecthood. The response interaction and absolute opposed-versus-aligned effect remain separately interpreted outcomes. Apparatus construction and measurement work are explicitly outside the model, so neither matching nor a positive occupancy contrast establishes full physical efficiency. These limitations are consistent with a useful external causal assay; they need not invalidate it or prompt a new protocol search.

## Completed continuation readout check

The separate `support-readout-physical-verification.json` now verifies all **384 continuation histories and 640,978 events** against the independent reconstruction. Every global/local balance error is zero and every production-readout field matches. This excludes the 64 preparation prefixes and does not replace their physics audit. The synthetic-only report is preserved and still correctly reports zero physical histories.

There are 64 preparations per arm; 42 contain an initial near-square R cohort and 22 do not. Descriptive continuation counts are:

| Response | Lesion | Original cohort cleared, of 42 present | Cleared with global R descendants at horizon | Cleared with new cohort copies afterward | Cohort outward/inward crossings | Mean global/local cohort R-time, over all 64 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Opposed | Sham | 26 | 10 | 13 | 140 / 77 | 178.38 / 96.74 |
| Opposed | Near | 33 | 15 | 17 | 161 / 94 | 165.77 / 80.45 |
| Opposed | Remote | 32 | 14 | 14 | 167 / 105 | 176.61 / 92.17 |
| Aligned | Sham | 27 | 16 | 16 | 207 / 128 | 224.23 / 112.70 |
| Aligned | Near | 31 | 11 | 14 | 167 / 105 | 168.40 / 80.53 |
| Aligned | Remote | 29 | 10 | 13 | 167 / 97 | 180.78 / 86.21 |

These distinguish three facts. Original chemical instances do turn over while descendants sometimes persist and copy further. Descendants also spend substantial time outside the frozen readout region. In two opposed/sham continuations, no R remains locally at the horizon while descendants of the prepared near cohort remain globally. Thus regional disappearance cannot be used as lineage extinction. The last count is the specific sufficient condition `final_local_R=0` with `final_cohort_size>0`, not a comprehensive count of every lineage located entirely outside a region that still contains unrelated R.

Conversely, 12 opposed/sham continuations end with global R present but zero descendants of the initial near cohort. Such R can come from other prepared cohorts or later sources; this readout alone does not identify it all as thermal immigration. Local occupancy and the marked cohort's continuity are separate measurements.

Attributed cohort net L currents are respectively 535, 428 and 535 in opposed sham/near/remote, and 695, 443 and 572 in aligned sham/near/remote. No reverse production was attributed to these cohorts or the near helper sites in these continuations, so their net counts happen to equal forward counts. Reverse channels are retained; the explicit synthetic fixture, rather than an observed reverse count in this dataset, validates the attribution branch. The counts identify catalytic provenance, not ownership or a guaranteed beneficial return.

Post-clearance integrals and cleared-with-descendants counts are descriptive. Lesions change whether and when clearance occurs, so these conditional windows cannot independently rank regeneration efficiency or establish a causal improvement after turnover. The predeclared paired occupancy contrasts, audited separately, remain the causal endpoints.

The verifier is frozen at SHA-256 `ec20c32a96eea92d9a581bbe4f214c0f5215850c101fe07d8632457be03c5cb5`; the physical readout report is `1501069fd2f3f0d7401145cda7b201329498a904ac4949f1632e38180dc3abc0`. No primary study source or natural history was changed.
