# Prospective finite-domain phase study

This first bounded study tests one supplied physical regime. It has no outcome-selected parameter pilot. A failure will remain a failure at this regime; the protocol does not authorize a rescue scan. Software smokes and independent physics checks precede execution. No physical outcome has been inspected when this design is first written.

## Chemistry and preparation

Use the separately audited `model-spec.md` and `active_phase_model.py`. Candidate values are theta=1, epsilon_B=2, J_BB=2.4, J_EB=2.5, mu=12, material mobility=1, enzyme mobility=.5, passive rate=.02 and catalytic rate=.2. These are supplied dimensionless interactions, not fitted material constants or a radius prediction. All trajectories start with zero B, randomized soluble A and randomized mobile E. Both reversible reaction channels remain available. No spatial group readout changes a physical rate.

Compare side16 (77 A/B carriers, 8 enzymes) and side32 (308 carriers, 32 enzymes), keeping inventories per lattice area exactly equal. At each size use three arms: full drive/affinity; no drive (mu=0, same Hamiltonian); and no enzyme-B affinity (J_EB=0, same drive). The last changes the Hamiltonian and is not a purely kinetic ablation. The enzyme inventory is conserved rather than generated or repaired, and ideal homogeneous fuel/waste reservoirs are supplied. Their chemical work is counted; reservoir construction is outside the model.

Use eight fresh integer seeds 720100–720107 at all size/arm combinations (48 physical records). A seed is one coupled replicate across arms/sizes; accepted events and correlated time samples are not independent replicates. Matched arms at one size have identical initial microscopic arrays. Run to T=400 with sampled arrays every5 time units, including endpoints. Early window is100–200 inclusive; late window is300–400 inclusive. Each window has21 snapshots, first reduced to one arithmetic mean per seed and metric. Full complete accepted-event ledgers remain available for independent physical replay. No trajectory or extinct/unformed state is removed.

## Analysis-only observables

`morphology.py` uses periodic nearest-neighbor B components, B-neighbor contact excess above the exact fixed-count random-placement expectation, enzyme-B excess association and a low-frequency spectral moment. The low-frequency length is 2*pi*sum(S)/sum(k*S), using0<|k|<=pi/2, with lattice spacing1. It is a descriptive wavelength with a declared cutoff, not a droplet radius or an independently derived length. Clusters of at least9 sites are an operational diagnostic; the cutoff supplies no natural dynamics. Undefined quantities are retained as null values.

Late dense clustering requires all three mean criteria at both sizes: B density>.05; fraction of B mass in components>=9 sites>.5; and B-neighbor excess>.15. Full late/early spectral growth is the ratio of the two within-seed window means. Compare this growth to the seed-matched no-drive growth separately at both sizes. Compare full late wavelengths and largest-B-component fractions between side32 and16. The affinity-free arm is a descriptive physical control, with the same complete reporting as other arms.

## Primary decision rule, fixed before outcomes

Use two-sided95percent Student t7 intervals over the eight seed-level values. Call the joint finite-horizon/domain arrest claim met only if every following gate holds:

1. At each size, lower bounds of the three full late dense-clustering metrics exceed .05, .5 and .15 respectively.
2. At each size, the upper bound of full late/early wavelength growth is <=1.25, and the lower bound of matched no-drive growth minus full growth is positive.
3. The upper bound of full late wavelength32/16 is <=1.5, and the lower bound of full largest-B-fraction16 minus32 is positive.

The tolerances1.25,1.5,9,.5,.15 and.05 are declared operational screens, not physical theorems. An undefined required metric in any seed makes the affected primary gate nonestimable and the joint claim not met. It is never a reason to drop that seed. Report all component gates and observed patterns even if the joint claim fails. A fixed10000-resample percentile bootstrap over complete seed rows is a sensitivity check and cannot override the primary rule. The statistic of ratios is explicitly the mean of seed-level ratios, not a ratio of pooled samples.

Formation and useful intrinsic size are separate. Even a met gate would be evidence only over these two domains and this finite time, leaving longer transient coarsening as an alternative. A failed gate does not disprove all driven-condensate mechanisms. No gene, copying variation, new functional role, heredity, aging, recursive individual or agency exists in this phase-only generator. Reaction currents/chemical state cycling and conserved carrier exchange are descriptive turnover diagnostics, not proof of those properties.

## Audit and stopping

Freeze model/spec/readout/protocol digests after independent rate, energy, conservation, detailed-balance, zero-drive and microscopic-ledger preflight. Run one timing-only software trace before the study if needed; do not use it to choose physical values. Source defects require a recorded repair and new preflight before outcome execution, not silent changes to a frozen generator. Report every physical record and all work/heat counters. No extra seeds, times, parameters or favorable windows are authorized by this study. A separately proposed study would require its own prospective protocol.
