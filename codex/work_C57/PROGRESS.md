# C57 progress

[read] Background, C51 scalar criterion and collaborator thermal notes read. The task is estimate (b), with ordinary off-half-filled weak XXZ first, then the retarded and multichannel extensions. All writes stay in this folder.

[read] The collaborator already supplies a fixed-history single-scale determinant/winding lemma, a conditional thermal-history recursion, and a thermal summation theorem. These are inputs to audit and cite, not new results to reproduce. C55, C56 and C58 PROGRESS.md were absent at the first check.

[in progress] Primary-source audit of Benfatto–Mastropietro and related free-energy results, and of rigorous XXZ quantum-transfer-matrix results. No interacting thermal coefficient or estimate P is yet claimed.

[read, new literature lead] Rutgers 127th Statistical Mechanics Conference, short-talk B12, 17 December 2024, Bruno Renzi with Alessandro Giuliani, Vieri Mastropietro and Venanzio Alessio Santini, announces "Free energy universal finite temperature correction in non-solvable quantum chains". Primary abstract: https://www.math.rutgers.edu/images/127_SMC_short_talks.pdf . No full manuscript or precise model/estimate located yet. Do not claim the question has no prior treatment; equally, this abstract cannot be invoked as the missing theorem.

[read] Giuliani–Mastropietro 1211.3030 equations (3.84), (3.90), (3.91) give a directly analogous scalar/history finite-size estimate for weakly perturbed Ising, whose quartic localization vanishes. BM 0907.2837 Theorem 1.1 is an exponent/susceptibility result, not a thermal-pressure result. MP 2108.13708 Theorem 8.5, (8.61), is reference beta cancellation; Section 8.2 constructs an analytic free energy, but does not state estimate P.

[in progress] Deriving a four-history bootstrap that keeps both the irrelevant defect and thermal powers. Bounding lattice and reference thermal histories separately loses the decaying defect factor. Also calculating the first perturbative coefficient in the off-half-filled XXZ pressure as a normalization and ensemble check. These are not an all-order interacting theorem.

[proved partial, completed write-up] PROOF.md now gives the four-history
bootstrap, including its Hessian term, a weighted-history version, terminal
scalars, and the explicit unverified GN hypotheses. It also derives the
first XXZ interaction coefficient at fixed field and fixed density. The
fixed-density answer is T^2/(12t^2)+O(T^4), agreeing with the first-order
dressed velocity. Higher orders are not controlled after division by T^2.

[read] SOURCE_AUDIT.md records the verified primary estimates and the
Rutgers announcement. No published estimate P for XXZ was located in the
inspected sources. C55 and C56 still leave the all-scale retarded theorem
open in their progress reports; C58's ensemble work is independent.

[numerical] Final checks.py run passes. Direct quadrature converges to both
ensemble coefficients. A nonlinear four-history identity is checked using
90-digit arithmetic. The preliminary finite-temperature tolerance and
binary-float cancellation failures are recorded in PROOF Section 9.

[status] RESULT.md and STATE.md deliver the partial result. The requested
fixed-interaction XXZ theorem, estimate (b), and its retarded/multichannel
extensions remain unproved. No counterexample or obstruction to the
physical thermal prediction was found.
