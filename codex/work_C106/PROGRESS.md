# C106 progress

2026-10-29.

[Proved work record] Read the project STATE.md, the complete bound background including both updates, REPORT.md sections 9.68 to 9.75, the current STATUS.md conclusions, and the relevant C69, C100, and free-fermion results. All writing stays in work_C106. The project STATE.md is outside the authorized folder and will not be edited.

[Assumed setting] The infrared sectors and the gauge-invariant scale M are those in the background. A local real-time realization is needed for the inherited velocity cap V = 6 e M.

[Proved preliminary reduction] A is the integral of N(u). Dyadic velocity windows improve the entropy estimate by replacing the largest velocity with the upper end of the window. A crossover law T(u) proportional to M (u/M)^p would give an integrable tail for p < 2, provided the selected sectors have a simultaneous one-sided entropy lower bound.

[Open investigation] A two-band nearest-neighbor free model with an avoided crossing appears to keep its Fermi velocities bounded away from zero while its thermal matching window collapses. The dispersion and the entropy limits will be derived before making a counterexample claim. This concerns the proposed crossover lemma, not the velocity-weighted bound, which is already proved for free chains.

[Open investigation] Operator dimension counts need both quantitative access to each independent infrared field and a relation between the number of fields and central charge. A stress tensor alone supplies neither relation. Entropy increments over disjoint temperature windows may close the borderline quadratic crossover case under a stronger, explicit matching assumption.

[Proved completed result] The avoided-crossing example is an explicit d = 4 nearest-neighbor chain. The bond norm is gauge-optimized exactly, and a smooth tuning holds v/M fixed. Direct entropy integration proves a fixed-fraction CFT matching endpoint of O(epsilon). Thus a positive onset bound depending only on d, v, and M fails, even within the free class. The example obeys A/M = 2r and does not refute the conjecture.

[Proved completed reductions] Wrote the integrable-tail and logarithmic-tail identities, the dyadic thermal theorem for p < 2, the stronger entropy-increment theorem at p = 2, and the quantitative operator-support lemma that would imply the target. Also proved that the specific universal thermal matching premises fail in the two-band example. The operator premises remain unproved.

[Numerical verification] The initial executable name python was unavailable. The first python3 run caught an incorrect intermediate expression for the velocity-squared root equation. Correcting it to (4t^2 x-epsilon^2 x^2)/(1+x) made the fixed-r assertion pass. Final checks include bond eigenvalues, equal extremal one-site marginals, entropy quadratures in both regimes, and entropy increments. All assertions pass; CHECKS.md records the final output.

[Proved literature record] Used only the browsing tool for two primary papers by Zou, Milsted, and Vidal on lattice-operator matching and numerical RG flows. They provide context, not a premise for the new analytic results. No claim of an exhaustive literature exclusion is made.

[Open final status] No general interacting N(u) estimate was proved. The concrete outcome is a proved obstruction to the proposed uniform crossover lemma and a precise statement of the missing operator-content/support hypotheses. RESULT.md, PROOF.md, and local STATE.md are complete. No writes outside this folder, no heavy numerics, no cluster jobs, no commits, and no pushes.
