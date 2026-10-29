# C116 progress

2026-10-29.

[Assumed] The setting and claim tags are those of specs/_bound_background.md, including every update through round six. Read REPORT.md sections 9.68 through 9.78 and the project memory. Work is confined to work_C116. No symmetry assumption is added to the target.

[Proved] A same-aspect-ratio comparison needs coefficients that cancel the extensive term without canceling the CFT torus function. The unweighted ratio Z(L,beta)/Z(2L,2 beta) does not do this.

[In progress] Derive exact partition-function comparisons and their scaling limits, audit spatial transfer dimension and reflection positivity, and run small full-spectrum checks on XX, Heisenberg and SU(3). Heavy numerical work, if needed, will be submitted to Typhon through SLURM.

[Open] No crossover-free bound on A has yet been obtained.

[Proved] For R(L,beta) = 2 log Z(L,beta) - log Z(2L,beta), the CFT scaling limit has large-aspect-ratio slope pi A/4. An all-scale bound R <= K(d) + k M f(d) beta/L would imply A <= (4k/pi) M f(d). The ordinary gluing estimate gives only |R| <= 4 beta M.

[Proved] A positive spatial transfer matrix of temperature-independent dimension D would give R <= log D, which fails already in XX. A two-site gate's operator positivity does not imply positivity after rotating its indices into a spatial transfer matrix.

[Numerical] Small full-spectrum checks completed in under one second, with largest charge block dimension 90. XX and Heisenberg rings through L = 8, SU(3) through L = 6. The XX full-trace product formula agrees with exact diagonalization to 5.7e-14 and permits O(L) checks through L = 512 without large diagonalizations. No heavy numerics ran locally.

[Blocked computation] First Typhon SSH attempt timed out during banner exchange. A retry is being attempted after continuing analytic work. A SLURM-only larger-check script is prepared, but no job has been submitted.

[Proved] The transfer-rank obstruction can be strengthened without positivity. Correct-pressure finite transfer representations of XX require log D_beta >= beta (sqrt(2)-4/pi) - 4 log 2. A microscopic d^2 gate rank cannot be promoted to a temperature-independent spatial transfer dimension.

[Proved] The corrected same-shape combination Q = 4 log Z(L,beta) - log Z(2L,2 beta) is contaminated by reciprocal velocities. A two-XX-sector family fixes d = 4, M = 1/2 and A = 1 while Q in the scaling limit diverges as one velocity tends to zero. R cancels that term. Full traces at L = 2048 illustrate the distinction with O(L) work.

[Proved] A gapped spectator with finite vacuum multiplicity and the stated ring gap and vacuum-energy conditions contributes no slope. An exact massive Ising spectator obeys 0 <= R <= log 2. A twisted XX example shows R < 0 on some even rings, so all-length positivity would add an unjustified assumption.

[Numerical] check_comparisons.py verifies a nonsymmetric complex-bond comparison, scalar and coboundary invariance, Golden-Thompson, a twisted XX sign counterexample, pressure subtraction, and the two-velocity example. All checks passed. Largest additional matrix dimension is 256; the script completes in less than one second.

[Blocked computation] The later Typhon retry also timed out during banner exchange. No job was submitted and no remote files were created. Larger diagonalizations remain unperformed; they were not run locally.

[Work record] RESULT.md, PROOF.md, CHECKS.md and local STATE.md written. The universal target remains open. The precise missing estimate and its implied constant C = 4k/pi are stated. No independent referee review was performed.
