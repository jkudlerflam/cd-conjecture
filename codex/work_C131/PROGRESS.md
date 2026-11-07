# C131 progress

2026-11-07. Run finished.

[Proved] The compactness reduction, explicit even-ring/all-length/admissibility-graph distinctions, and conditional finite-cover assembly are in PROOF.md. The requested table gives precise local lemmas and distinguishes proved inputs from assumptions. The gapped route is split into ring-energy convergence, stability, and coexistence versions. No analytical attempt was made on (d) or (e).

[Source audit] Hastings-Koma clustering does not compare ring energies. Michalakis-Zwolak stability requires Local-TQO and Local-Gap in addition to frustration-freeness. Borgs-Kotecky is a classical cubic-torus theorem; Borgs-Kotecky-Ueltschi explicitly assumes spatial dimension >=2. A zero-temperature quantum-ring corollary was not established. C130 remains a parity-preserving draft, with no ambient radii from which to estimate a coverage fraction.

[Numerical] Completed 12 starts each at L=2,3,4 with 1000-evaluation caps, then reran the same starts with 4000-evaluation caps. Best values found are 4, 2.0369383345, 2.7912527285 respectively. Reported convergence, residuals, and numerical gauge brackets are in NUMERICS.md. The final gauge check retains the best feasible dual iterate; this repaired loose bounds from failed final SLSQP steps. It did not change the best values. All ring-construction, gauge-invariance and XX baseline checks passed. Heisenberg benchmarks were evaluated through N=16.

[Deferred as requested for large computation] Full multistart at L=2,4,6,8 and odd L=3,5,7 is in submit_c131.slurm, 32 starts per length. It was not submitted, and no Typhon connection was attempted. No claim about suprema at L=6,8 or large-L drift is made.

RESULT.md, PROOF.md, NUMERICS.md, RUNNING.md and the local STATE.md are complete. C(2)<infinity remains unproved. Only work_C131 was modified; no shared files, commits, pushes, or external messages.
