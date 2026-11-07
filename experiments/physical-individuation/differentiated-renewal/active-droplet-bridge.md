# A physical route from renewed flow to finite spatial units

The current lattice experiments test response differentiation and local support. They do not yet establish a dynamically selected individual size. Local attraction and copying alone do not guarantee that several finite organizations persist rather than disappear, merge or occupy the whole available domain. Which behavior our exact generator exhibits requires further analysis; no coarsening theorem is asserted for it here.

There is relevant prior work that makes the phase-transition intuition concrete. [Zwicker and colleagues, Growth and division of active droplets provides a model for protocells](https://www.nature.com/articles/nphys3984), with an [open preprint](https://arxiv.org/pdf/1603.01571), model phase-separated droplets maintained by chemical conversion and external fuel. Their analysis and simulations show regimes of finite droplets and shape instabilities leading to division. This supplies a physical model of spatial units; it does not establish complete organisms or open-ended biodiversity.

[Kirschbaum and Zwicker, Controlling biomolecular condensates via chemical reactions](https://arxiv.org/pdf/2103.02921), develop a thermodynamically consistent reaction-and-phase-separation model. Size control requires more than merely adding reactions: energy supply and different reactive behavior inside and outside can matter, with enzyme partitioning supplying a physical source of that difference. This is a candidate mechanism for emergent localization, not permission to put an observer's droplet label into our lattice rates.

## An explicit conditional size calculation

The following thin-interface calculation exposes the competing geometric terms. It is an illustrative approximation to this family of mechanisms, not a derivation of our current CTMC or a claim to a new droplet theory.

Assume molecular interactions have already supplied two coexisting phases, with concentration difference Delta c>0. For an approximately spherical droplet of radius R, its excess material is (4*pi/3)*Delta c*R^3. Assume a diffusion-limited influx 4*pi*D*(Delta c_sup*R-Gamma), where positive Gamma captures a curvature correction, and net volume-proportional conversion out of droplet material (4*pi/3)*k*c_in*R^3. This net conversion and outside replenishment must be supported by a stated driven chemical cycle. The assumed separation of diffusion, interface and reaction scales is essential.

Mass balance gives

\[
\dot R=\frac{A}{R}-\frac{C}{R^2}-BR,
\quad A=\frac{D\Delta c_{sup}}{\Delta c},\quad
C=\frac{D\Gamma}{\Delta c},\quad
B=\frac{k c_{in}}{3\Delta c},\qquad A,B,C>0.
\]

Stationary positive radii solve g(R)=A*R-C-B*R^3=0. Its unique positive maximum is at Rm=sqrt(A/(3B)). There are two positive radii exactly when

\[
\boxed{4A^3>27BC^2.}
\]

The smaller radius is unstable and the larger stable: at a stationary radius, the radial derivative is A/R^2-3B. Equality is the fold threshold. For A=3,B=1,C=1 in consistent dimensionless units, the two radii are 2*cos(4*pi/9) and 2*cos(2*pi/9), approximately .347296355 and 1.532088886, with radial derivatives about +21.87 and -1.72. [certify_radius.py](certify_radius.py) and [radius-certificate.json](radius-certificate.json) enclose the two roots and derivative signs with 50-decimal outward interval arithmetic. Material can turn over at the stable size while the size remains. The constants' dimensions differ (A is length squared per time, C length cubed per time, B inverse time), so this threshold is dimensionally consistent.

There is no preferred group size or biodiversity reward in that radius equation. Finite size follows from the competition between supply, curvature and consumption geometry. But neither spherical symmetry nor a preexisting phase interface was derived by this reduction. Stability to nonspherical modes, spontaneous nucleation, division, inheritance and interactions among droplets need the full field or particle model. A scalar radius cannot prove them.

## Consequence for this research program

This offers a concrete way to connect the threads if the current lattice establishes differentiated renewal but still lacks spatial integrity: use physically segregating, fuel-driven material with an intrinsic reactive length scale, then ask whether metabolic response classes help maintain or regenerate those units. Chemical potentials, reverse channels and input work must remain consistent. The difference between internal and external behavior should arise from molecular composition or partitioning, not a controller that recognizes a group after the fact.

The published droplet models are evidence that a phase transition and sustained flow can generate finite spatial organization under specified chemistry. They are not evidence that least action universally maximizes biodiversity, that every such droplet is an ethical subject, or that our own model has already joined those results. A successful Unitology defense must still distinguish preserved participant capacities from a whole's efficient propagation.
