# C118 handoff

2026-10-29.

[Open] No universal K(d), f(d), k proved for the requested ratio bound. No counterexample constructed. Parent STATE.md is intentionally not edited because the task forbids writes outside work_C118.

[Proved] Periodic XX, even L >= 4, all beta >= 0: R <= 2 log 2 + (5 pi/3) beta/L. For h = J h_XX this is k = 10 pi/3, f = 1, M = |J|/2. Proof retains both Jordan-Wigner parity sectors and bounds the momentum quadrature error by integral |f''| <= 4 beta. Independent XX products have K(d) = 2 log d and the same k, f = 1. This does not improve the already known bound on free-chain A.

[Proved] Diagonal classical bonds, including oriented ones: R <= log d using entrywise nonnegative transfer matrices, without matrix positive semidefiniteness.

[Proved] All-beta target implies L(E_(2L)-2E_L) <= k M f(d) at each compatible L. Averaged seam has norm >= 2M on XX, Heisenberg and SU(3), at any L. Partial-swap interpolation lacks continuous symmetry, with explicit commutator witnesses. R is neither convex nor concave already on XX L = 4. A bound at a single tau plus continuum scaling data cannot account for slow-sector multiplicity without additional microscopic information.

[Proved] R = beta<V>_rho + I + two marginal relative entropies. The relative entropies are not capacity bounded. The reverse-direction identity is R = beta<V>_sigma-D(sigma||rho). No beta/L bound emerges from positivity alone.

[Numerical] Local checks complete in about two seconds. Full ED through XX/Heisenberg length 8 and SU(3) length 6; linear-cost XX paired lengths through 8192. All proved comparisons pass. Convexity, concavity, R <= MI, and averaged seam norm <= 4M/L fail. A sampled trial bound with K = 2 log d, k = 10 pi/3, f = log_2 d passes but is unproved for interacting chains.

[Unperformed] SLURM script prepared with 64 GB, one CPU, 24 hours for full spectra through length 12, including SU(3)'s block of dimension 34650. Runtime is not benchmarked. No Typhon connection, submission, or external write. No independent review.

Files: RESULT.md summary; PROOF.md derivations; CHECKS.md numerical findings; PROGRESS.md work record; checks.py and ring_tools.py implementation; checks_local.json full output; run_summary.json saved console summary; run_typhon.slurm unsubmitted larger run.
