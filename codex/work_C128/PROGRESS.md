# C128 progress

2026-11-06.

[Assumed scope] All writes and Python execution remain in work_C128. No network shell command, cluster connection, or parallel agent is used. Project-level STATE.md cannot be modified under the folder-only rule; a local handoff will be written.

[Proved, inherited] Read the shared target, C127 specification, full bound background and reconnection note, REPORT sections 9.76 through 9.88, C124 RESULT, and C125's proof, results, checks, handoff and source scripts. Existing numerical artifacts are not used for orientation. The positive transfer is naturally a two-column construction, and XX already excludes an all-frequency linear crossed cone.

[Heuristic under audit] A nonlinear fermionic product from the exact XX crossed spectrum should yield a direct normalized-trace bound through convex midpoint quadrature. Independent Poisson domination of neighboring column histories has a more basic difficulty: paired jumps have identical continuous times, making the joint history weight singular relative to independent column clocks.

[Proved] Completed the explicit nonnegative decorated-history contraction and the separate positive-operator construction. Proved a Poisson generating-function bound for the count of incident events for all requested XXZ models. Independent column-clock domination fails because shared bond jumps are singular in the product of independent event-time measures. Uniformly bounded conditional intensities also fail on explicit bridges.

[Proved] Derived the XX normalized-trace bound by convex midpoint quadrature, including its exact integral, an improved ratio coefficient and a literal scalar Fock product. The lower denominator bound retains the crossed vacuum normalization. The unchanged XX majorant fails for Heisenberg already through its exact four-site zero-temperature slope.

[Assumed and open] Stated a sufficient interacting many-state comparison using the valid nonlinear lattice frequency rate. Proved that its logarithmic ultraviolet tail still gives order beta M/L. The many-state comparison itself remains open; no counterexample to all scalar Fock majorants or to finite C(d) was found.

[Numerical] Completed 192 finite-Trotter trace tests, 35 exact XX trace tests with cutoff error bounds, 36 local event-count generating-function checks, explicit gate and paired-jump tests, and frequency/integral tests. Recorded 24 failures of the unchanged XX bound on interacting finite-Trotter spectra. Corrected an overflow in the initial tail calculation; the final suite passes its stated assertions.

[Work record] RESULT.md, PROOF.md, CHECKS.md, STATE.md, checks.py, checks.json and run_typhon.slurm are complete. The extended cluster scan is prepared but unsubmitted. No outside file was modified.
