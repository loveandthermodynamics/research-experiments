# Independent review of the conditional radius bridge

Reviewed 2026-10-21. Verdict: **passes as a conditional three-dimensional radial calculation and accurately bounded literature bridge. No mandatory correction identified.** No new physical experiment or trajectory was run, and no reviewed source or result was edited.

Reviewed SHA256 digests:

- `active-droplet-bridge.md`: `7171682ea30161e4231e5abeeaa0c06d78de2f903996ea581c8da658108d476e`
- `certify_radius.py`: `a227ccb835ed6a4235ac6c244ae7d02749f65d29415b413e6652a73230e5f88c`
- `radius-certificate.json`: `a4b1625eccbebf3e53d1a784b8abe1d566a9dc13dda5d1e915af003f5553e9fe`

Differentiating the supplied excess material gives `4*pi*Delta_c*R^2*Rdot`. Dividing the stated diffusion influx minus volume sink by that factor gives exactly `Rdot=A/R-C/R^2-B*R`, including the factor `1/3` in B. With concentration measured in amount per volume, Gamma has units concentration times length; A, C and B consequently have units length squared/time, length cubed/time and inverse time. Both sides of the fold inequality have units length to the sixth/time cubed.

For positive A, B and C, `g(R)=A*R-C-B*R^3` starts negative, has its sole positive stationary maximum at `Rm=sqrt(A/(3B))`, and eventually becomes negative. The maximum is positive exactly when `4*A^3>27*B*C^2`, giving precisely two positive roots. At either root, differentiating the radius flow and substituting `C=A*R-B*R^3` gives `f'(R)=A/R^2-3*B`. Thus the smaller root is radially unstable and the larger radially asymptotically stable. Equality is a nonhyperbolic fold, not an ordinary stable radius; the draft does not claim otherwise.

I independently evaluated both saved brackets using exact rational arithmetic rather than the certificate's interval library. For A=3, B=C=1 the margin is exactly 81. The lower-root bracket `[.347296355333,.347296355335]` has endpoint g signs negative/positive and derivative in approximately `[21.87257810798,21.87257810827]`. The upper bracket `[1.532088886237,1.532088886239]` has signs positive/negative and derivative approximately `[-1.721933856721,-1.721933856717]`. These agree with the outward enclosures. Monotonicity on either side of Rm=1 gives uniqueness within each bracket. The certificate's embedded source digest matches the actual source.

The prior-work claims agree with the primary papers. [Zwicker et al.](https://arxiv.org/html/1603.01571v1) explicitly distinguish radial size stability from shape instability, and exhibit driven growth/division in a continuum model. [Kirschbaum and Zwicker](https://arxiv.org/html/2103.02921v1) supply a thermodynamically consistent reaction/phase-separation theory in which driven reactions and enzyme localization can control droplet size. The latter also explains why arbitrary fixed reaction laws cannot automatically be treated as thermodynamically consistent.

The present certificate validates only the stated reduced equation. It supplies neither a chemical-potential derivation of its flux assumptions nor a parameter map to the existing two-dimensional lattice. Fixed far-field supersaturation, approximately constant coexisting concentrations, a thin interface, spherical geometry and the assumed reaction/diffusion scale separation are conditional inputs. A physical application must check these approximations at its selected radius. Radial attraction does not prove nucleation, nonspherical stability, division, coexistence of several droplets, heredity or individuation. The bridge explicitly retains these limitations and does not rescue the failed response comparison or derive a universal biodiversity or moral principle.
