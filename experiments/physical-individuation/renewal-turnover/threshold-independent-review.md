# Independent review of the open R,L threshold model

The open well-mixed model is a thermodynamically consistent conditional construction of reciprocal material maintenance. Its finite-temperature high state has certified local stability, consumes positive chemical work, and is unavailable at stationarity under either specified production/protection ablation. This is a substantive mathematical result under supplied chemistry. It does not validate the unsuccessful finite-fuel spatial renewal claim or establish organismhood, diversity, aging, recursive individuality or moral legitimacy.

## Checks performed

I inspected `verify_threshold.py`, `threshold-verification.json`, `interval_threshold_certificate.py`, `threshold-interval-certificate.json`, `finite-budget-bound.md` and the corrected `theory-bridge.md`. Both numerical programs reran successfully with file writes intercepted in memory; their generated JSON matched the saved reports. No root-owned source or output was changed.

Separately, a 75-digit Decimal reconstruction, using explicit mass-action currents rather than importing either verifier, found the low, saddle and high stationary R concentrations:

| State | R | Jacobian trace | Jacobian determinant |
| --- | ---: | ---: | ---: |
| Low | 3.91619708093245e-8 | -1.05000000412 | 0.049999862748 |
| Saddle | 0.0345629856573734 | -1.06912602203 | -0.038872197082 |
| High | 0.232847617046128 | -1.46569524131 | 0.179287877281 |

The stationary reduction, analytic Jacobian and stability classifications agree. Four independently reconstructed nonstationary free-energy balances had absolute discrepancies at most 1e-74 at this precision. This is a numerical arithmetic check of the identity, alongside its algebraic derivation below.

## Thermodynamics and the interventions

With theta=1, epsilon_R=epsilon_L=20, epsilon_D=epsilon_W=0, epsilon_F=40 and fuel/waste concentrations both 1, eta=exp(-20) is consistent for both productive reverse coefficients and both passive reverse coefficients. For example, the productive R-channel affinity is `20 + log(D/R)` and the passive R-to-D affinity is `20 + log(R/D)`. Their cycle affinity is 40. The L channel has the corresponding relation. The L-dependent passive mobility multiplies both directions and cancels from the affinity.

For ideal-mixture material free-energy density, the independent gradient components are `20 + log(R/D)` and `20 + log(L/D)`. Substituting the four oriented reaction currents gives exactly

    Gdot_material = 40*(J_copy + J_L) - sigma,

where sigma is the nonnegative dimensionless current-times-log-flux-ratio sum. At the high steady state the chemical power and sigma are both approximately 10.8358710840 per normalized volume and time. These are continuing chemostat costs, not extracted work or a minimum-dissipation principle. Reactor walls, mixing and reservoir-maintaining apparatus remain supplied idealizations.

Stopping L production removes both productive directions while retaining thermal D↔L conversion. Removing protection changes a paired kinetic factor while retaining the same species, material inventory and well energies. Starting these comparisons from the same high state is a matched intervention; no group score enters the equations. The stationary bounds are valid for all physical stationary states of the ablated systems:

    no protection: R <= eta*M/(1-M),
    no L production: R <= eta*M/(exp(-4*eta*M)-M).

For the second bound, stationarity implies L=eta*D and hence the protection factor is at least exp(-4*eta*M)>M. Both ceilings are approximately 3.9162e-8, many orders below the certified high R concentration. The finite-time ODE endpoints are approximations to this small baseline, not the bounds themselves.

## What the interval certificate proves

The outward interval calculation encloses a sign change in each of three positive physical R brackets. Along the exact L-nullcline, the scalar derivative equals `det(J)/g_L`. Its definite sign proves local uniqueness within each bracket. Negative trace and positive determinant certify the two locally asymptotically stable states; negative determinant certifies the saddle. These arguments are valid and do not require the preceding grid search to be globally exhaustive.

The certificate does not exclude additional equilibria outside those brackets or establish a global finite-temperature basin classification. Numerical integrations illustrate two basins. Their supplied initial L does not establish formation from zero L under a specified founder/preparation budget. The negligible-reverse threshold proof is stronger analytically within its stated approximation; it must remain distinct from the finite-temperature completion. For complete notation in the general threshold proposition, `min` over L≥0 or `inf` over L>0 avoids an unattained minimum when F is increasing; this does not affect the exhibited parameter regime.

The corrected distinction between linear stability and the nonhyperbolic equality case is appropriate. Likewise, passive founder-loss hazards are positive; finite-temperature productive reverse reactions add instance loss. Concentration renewal demonstrates chemical role turnover, not traced atom replacement, inherited organizational memory or a persistent subject. Copying ablation and controlled recovery are discussed analytically/prospectively; they are not additional interventions executed by `verify_threshold.py`.

## Separation from the spatial result

`finite-budget-bound.md` correctly derives `-dD(p||pi)/dt = sigma` and the finite integrated entropy-production budget for the autonomous reversible lattice, with communicating classes handled explicitly. Nonzero equilibrium R occupancy can give an unbounded occupancy-time integral without sustained nonequilibrium work. Neither occupancy-time nor that decay identity is an aging law.

The new construction changes to an open, well-mixed reactor with continuing fuel/waste control and a supplied exponential protection law. Its positive theorem complements the spatial negative result; it does not overturn it. In the spatial turnover confirmation, protected ongoing production remains inconclusive and the declared joint claim remains false. A finite molecular-count completion of the new model would also have noise and possible metastable switching rather than permanently separate deterministic attractors.

Verdict: the corrected mathematical bridge and certificate support conditional reciprocal maintenance within the supplied open chemistry. The full research objective remains open. No additional parameter experiment was run during this review.
