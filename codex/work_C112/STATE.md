# C112 session memory

2026-10-29

[Proved within CFT] Pure stress-tensor modulation has kappa = 1/(24 pi), with A = sum c_a v_a. Ring adaptation deficit is A G/(24pi)+(pi A/6)(1/ell_w-W/L^2). Fundamental cosine n = 1 has zero CFT deficit. For n fixed, D/G = (A/(24pi))(1-1/n^2). Gaussian line response is A epsilon^2/(48 sqrt(pi) sigma). Details and normalization source are in PROOF.md Sections 1-3.

[Proved] Literal all-representatives smooth-modulation lemma is false. For a = diag(-delta,delta,...), X = |0><1|+h.c., and h = (a tensor I+I tensor a)/2+t(X tensor I-I tensor X), the uniform chain is independent onsite spins, M = delta, gap = 2delta, A = 0, but the smooth deficit ratio tends to t^2/(2delta^2). This works for every d >= 2. Its minimizing representative has t = 0 and zero deficit, so it does not refute a specified-minimizer variant.

[Proved] At XXZ Delta = -1, the standard representative has D = 0 for every positive w. Adding t(Sz_i-Sz_(i+1)) leaves the uniform chain and M = 1/2 unchanged, but the fundamental cosine has small-amplitude ratio t^2 L/[4(L-1)sin(pi/L)^2]. A separate spin-coherent trial bound proves divergence at fixed epsilon as L grows, so an exchange of limits is unnecessary. The endpoint has quadratic magnons and no nonzero 1/L Casimir coefficient; approaching from the critical phase gives A -> 0. See PROOF.md Section 5.

[Numerical] Completed lightweight checks only, main suite 3.43 seconds plus a subsecond scalar check. XX through L = 128 matches the circle formula at the level of 0.0005 to 0.0008 relative error for n = 2 and a Gaussian. XXZ standard sectors through L = 12 have finite-size deviations, shrinking for XX, Delta = 1/2 and Heisenberg. SU(3) through L = 9 and SU(4) through L = 8 were also checked. Gauged ferromagnetic endpoint through L = 12 matches the analytic second-order coefficient.

[Numerical] Three full-space random bonds per d at L = 8,5,4 respectively were checked both raw and gauge-optimized. Largest raw ratios are 0.0773776745, 0.2207099063, 0.0766413 to 0.0766465. Largest optimized ratios are 0.0446512756, 0.1987820769, 0.0533731 to 0.0533768. These are smoke tests, not a moderate-size survey or thermodynamic phase assignments. M has explicit numerical primal and dual bounds. Largest eigenvector residual is 4.40e-14.

[Numerical] Exact onsite scalar calculation gives ratio 499999.999980 at d = 2,3,4 with t = 1, delta = 0.001. The supremum is analytically infinite. Scalar ferromagnetic trial-state bounds reach ratio >= 6495.97 at L = 8192, epsilon = 0.1. Neither scalar calculation is a large many-body diagonalization.

[Proved, infrastructure observation] Two Typhon SSH attempts timed out during banner exchange. No job submitted. cluster_job.sbatch and the guarded --cluster suites are prepared but unrun. No heavy work ran locally or on the login node. The guard was tested successfully.

[Assumed, unresolved] A revised lemma in a specified energy-density representative, its relation to the gauge-optimized M, and suppression or subtraction of additional microscopic gradient response remain open. The all-sectors bound on A is not proved or refuted. The requested moderate-size unrestricted survey remains blocked by Typhon access.

[Proved, file map] RESULT.md is the summary, PROOF.md contains derivations, NUMERICS.md describes reproduction and limitations, PROGRESS.md records the work. smooth_checks.py generates small_results.json; scalar_checks.py generates scalar_results.json. additional_results.json is an earlier subset. cluster_job.sbatch is unsubmitted. Parent STATE.md and every file outside work_C112 were left untouched under the task-specific instruction.
