# C120 progress

2026-10-29.

[Work record] Read the project memory, all background updates, REPORT sections 9.68 through 9.80, and the C116/C118 notes and source files. Parsed every saved JSON field in those folders. Their settled results are inputs. All research changes and Python execution are confined to work_C120. The parent STATE.md will not be modified under the task restriction.

[Proved] Dimension does not bound arbitrary-reference max-relative entropy, Petz or sandwiched Renyi divergences, or negative log fidelity. Support min-relative entropy vanishes between finite-temperature Gibbs states. Exact mixed-trace identities retain an unbounded normalization contribution.

[In progress] Derive conditional min-entropy bounds and the missing reference comparison, explicit positive-map trace counterexamples, and a signed-history transfer comparison. Test each proposed replacement on full ED of XX, Heisenberg and SU(3), with full-trace free-fermion checks at large XX length.

[Open] No universal ratio bound or counterexample to that bound has been obtained.

[Proved] Wrote complete derivations in PROOF.md. Conditional min-entropy gives R = Hmin(sigma)-Hmin(rho)+log[mu(Y)/mu(X)], with an entropy difference bounded by 2 log q. The remaining log ratio retains beta times the doubling-energy defect. Petz, sandwiched and fidelity identities likewise retain unnormalized mixed traces.

[Proved] Constructed a strictly positive completely positive qubit map for which the transfer trace ratio is unbounded at fixed superoperator dimension four. Derived a real-matrix comparison with an explicit negative-pair penalty. Map positivity, matrix positive semidefiniteness and entrywise nonnegativity must be distinguished.

[Numerical] Initial suite passed in about six seconds. Full ED covers 27 temperature/model pairs, Petz orders 0.05 through 0.9, conditioned sandwiched and max-divergence calculations, support limits, fidelity, trace distance and conditional-min-entropy witnesses. Also checked 24 finite Trotter history networks and full XX traces through paired length 16384. Ground-state Slater determinants extend state-overlap checks through paired length 1024.

[In progress] Final verification adds the doubled-ring temporal contraction as an independent check of the spatial transfer denominator, checks the pure-state XX overlap against ED, and records explicit failed candidates and numerical limitations.

[Numerical] Final suite passed in 10.35 seconds. Maximum Petz identity residual 1.92e-13. Spatial and temporal contractions agree to 6.22e-15 on L and 6.11e-14 on 2L. The ground-state XX determinant agrees with ED below 1e-12. The stable large-ring thermal recursion agrees with the inherited formula below 1.31e-11 on the calibration sizes. Failed inequalities are reported explicitly in CHECKS.md.

[Work record] RESULT.md, PROOF.md, CHECKS.md and local STATE.md are complete. Finite-temperature large-ring overlap calculations and full conditional-min-entropy optimizations were not performed. The large-ring numerical tests distinguish full thermal ratios from zero-temperature state-overlap checks. No independent review was performed. The universal ratio estimate remains open.
