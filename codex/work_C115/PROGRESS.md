# C115 progress

2026-10-29.

[Proved, work record] Read the project memory, the full bound background with all updates, REPORT sections 9.68 to 9.78, and every file in C106 and C109. Only this folder is writable for this task. The project memory will be represented by a local STATE.md handoff.

[Proved, spectral correction] A conformal level of dimension x in a sector of speed v costs 2 pi v x/L. An energy window 2 pi u x/L contains its tower through dimension u x/v, not through x when v >= u. Slow sectors contribute more states to a fixed low-energy window.

[Proved, obstruction under finite-level CFT matching] Two decoupled XX chains, one fast and one arbitrarily slow, have d=4 and can be normalized to M=1. At a fixed threshold u=1 their fast multiplicity is one, but their unfiltered ring excitation trace and low-energy counting function are unbounded across the family after the infrared limit. This is a failure of the proposed trace upper bound, not of the target bound on A.

[Open, current work] Develop a velocity-resolved spectral lemma with an explicit summation constant, distinguish full modular partition functions from vacuum descendant counting, and audit the dynamical channel-capacity premise against C106. Test the counting statements on free chains, XX stacks, and the recorded slow-mode constructions.

[Proved, completed spectral analysis] Wrote the strict infrared trace counterexample with the limits in the order L to infinity, retained oscillator count to infinity, then slow velocity to zero. Corrected the energy-window direction. Explained why a spatial quantum transfer matrix does not have a fixed d^2 column space and why a total E/|P| cutoff does not select fast constituents.

[Assumed missing input; proved reduction] For C_j in the dyadic velocity window below 4eM 2^(-j), C_j<=K(d)2^(j/2) gives A<=4e(2+sqrt(2))M K(d). A linear-j bound gives A<=16eM K(d). For full modular state content the dimension trace satisfies log Q_j(1/2)>=pi C_j/4. A new microscopic upper bound on a suitable velocity-resolved trace would close the proof. No such bound is asserted here.

[Proved, dynamical audit] Even a fixed-width packet premise plus SIE gives only a borderline 1/u multiplicity estimate per independent window. C106's shrinking pockets obstruct the proposed uniformly short IR packet construction. The required simultaneous capacity or multi-scale budget remains open.

[Proved and inherited numerical tests] Free N(u)<=min(log_2 d,2M/u); XX stacks have A=2M and obey the sufficient spectral inequality. Restored physical cell lengths in the slow-mode table. A vacuum hop changing both letters proves optimized M>=t for the spatial construction, so its leading payload velocities are indeed tiny relative to M on the recorded schedule. An all-n exponential envelope is stated conditionally, and rounded data are not treated as certified bounds.

[Numerical, completed verification] The first script launch stopped before any checks because the report initializer contained a stray unary plus on a string. Removed it and reran. All final assertions passed, including the subsequent gauge-witness check. The largest matrix was 16 by 16. CHECKS.md contains the passing output.

[Proved, work record] Browsed primary sources for the quantum transfer matrix, modular partition-function convention and SIE. Corrected the transfer-matrix authors against the paper metadata. No external source is used as evidence for a general counting theorem. No heavy computation or cluster access was needed.

[Open, final status] RESULT.md, PROOF.md, CHECKS.md, checks.py and local STATE.md are complete. The general interacting bound and a universal velocity-resolved spectral budget remain open. No files outside work_C115 were modified.
