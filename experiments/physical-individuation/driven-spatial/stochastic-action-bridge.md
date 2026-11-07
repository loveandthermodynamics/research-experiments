# A conditional stochastic-action bridge

A useful action principle follows from a supplied chemical generator. It does not independently choose that generator. [Lazarescu, Cossetto, Falasco and Esposito, Large deviations and dynamical phase transitions in stochastic chemical networks, §III, equations 29–30, 41–42 and 52–55](https://arxiv.org/pdf/1902.08416) derive the chemical-network Lagrangian/Hamiltonian structure from the large-volume master equation.

In dimensionless concentration units, let system size be Ω, directed reaction increments be \(\nu_\rho\), and propensities be \(\Omega a_\rho(x)\). For a fixed initial state the dynamical rate functional is

\[
\mathcal H(x,p)=\sum_\rho a_\rho(x)(e^{p\cdot\nu_\rho}-1),\qquad
\mathcal L(x,v)=\sup_p\{p\cdot v-\mathcal H(x,p)\},\qquad
I[x]=\int_0^T\mathcal L(x,\dot x)dt.
\]

Path probabilities have exponential scale \(e^{-\Omega I[x]}\), under the appropriate large-volume assumptions. Since \(\mathcal H(x,0)=0\), \(\mathcal L\geq0\). Convexity gives \(\mathcal H(x,p)\geq p\cdot\partial_p\mathcal H(x,0)\), hence

\[
\dot x=\partial_p\mathcal H(x,0)=\sum_\rho a_\rho(x)\nu_\rho
\quad\Longrightarrow\quad I[x]=0.
\]

Thus ordinary reaction-rate evolution minimizes this nonnegative fluctuation cost. This statement conditions on the reactions, mobilities, supply, initial state and scaling. It does not prove that their trajectories contain organisms. Different deterministic basins can all contain zero-action histories; the scalar cannot rank them without those physical conditions.

Our maintained high branch can be zero-action in exactly this sense: its copying, protective production and turnover solve its supplied open reaction equations. The earlier threshold proof identifies why that branch exists and requires continuing paid replacement. The stochastic formulation describes typical maintenance and the costs of atypical chemical departures; it adds no diversity reward or new selection law.

This Hamiltonian is a jump cumulant generator, not mechanical energy. Its conjugate p is a fluctuation variable, and the normalized action is a probability cost, not Hamilton's physical action measured in energy-times-time. Entropy production is another quantity: local detailed balance relates forward/reverse probabilities and chemical affinities, while symmetric mobilities still affect the entire generator and action. Zero fluctuation cost can coexist with positive continuing dissipation and supplied work.

Finite reverse chemistry must remain in the sum as separate directed channels. Omitting small uphill rates changes the generator and potentially its absorbing boundaries. Our reversible maintenance benchmark and source-free diversity theorem must therefore keep their distinct scopes.

The diversity theorem also retains a finite gate switching at order-one rate κ. Increasing chemical system size does not make that gate deterministic. Conditional on its realized history, the chemical Hamiltonian uses \(a_\rho(x,\xi_t)\), and the gate-driven rate equations have zero **chemical** action. Gate histories retain their own jump and waiting-time probabilities. This is not a single deterministic zero-noise limit for the joint system.

Finally, the small-forcing proof was checked against [Benaïm, Stochastic Persistence, version 3, §3.3 and Theorems 6.1 and 4.4](https://arxiv.org/pdf/1806.08450v3): compact common invariance, auxiliary L, a Feller switching process, invariant extinction facets, and positive weighted invasion at every ergodic boundary law match its hypotheses. Its conclusion is occupation-measure persistence for the source-free concentration process, not finite-population survival or finite-reverse lineage retention. That conditional functional-retention result can follow zero-chemical-action evolution without establishing biodiversity maximization, individuality or moral legitimacy.
