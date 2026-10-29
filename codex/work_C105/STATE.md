# C105 session memory

2026-10-29

[Assumed, unresolved] General target A <= C M f(d) remains unproved. The intended stress-energy statement must specify the energy-density gauge and distinguish smooth energy from raw lattice oscillations.

[Proved] New XX certificate. For w_i = a_i a_(i+1), a of finite support,

    Delta(w) <= (1/(2 pi)) sum_i (a_(i+1)-a_i)^2.

It follows from -D_a|T|D_a <= D_a T D_a <= D_a|T|D_a and the trace-norm bound. With r_n = 2[Gamma((n+2)/2)/Gamma((n+1)/2)]^2 and a_i^2 = r_i r_(N-1-i)/N, the weights are x(N-x)/N and the allowance is (1/(4 pi))log N+O(1). This is weaker than the existing free theorem A <= 2M, but supplies a weighted lower-energy certificate.

[Proved] Positive uniform-interval decomposition obstruction. For XX, g(n) >= 1/sqrt(2)-2/pi > 0 for all n. Every positive decomposition has total coefficient at least max w. Its allowance is therefore at least g_* max w, even if all overlaps are allowed.

[Proved] Fixed-range local-marginal SDP obstruction. Full-rank rho_r permits sigma_r = (1+t)rho_r-tI/d^r with a lower bond energy and the same local overlap constraints. It supplies a feasible weighted energy below e0 W by a positive constant times W for every fixed r. XX marginals are full rank at every finite r. Growing range or nonlocal constraints are not excluded.

[Proved] Exact XX half-line formula is Equation (2.1) of PROOF.md. The first-bond deficit is 4G/pi^2-1/pi. [Numerical] At x = 8191 the raw scaled deficit is 111.42982189; the centered (1,2,1)/4 filter gives 0.01326197011, close to 1/(24pi). [Heuristic] The raw staggered contribution is predicted to be (-1)^(x-1)/(2pi x log x). A uniform spectral remainder bound is missing, so its asymptotic divergence is not labeled proved.

[Numerical] Completed small tests on XX, Heisenberg, XXZ(1/2), SU(3), all interval inequalities pass. Naive extension Delta <= -e0_centered D(a)/2 fails for Heisenberg at N = 12. [Proved] Reduced the failure to an exact N = 4 singlet witness. With a = (3/4,1,1,3/4), the trial energy is -sqrt(7)/4, below the proposed lower bound (25/8)(1/2-log(2)). Keep this failure visible.

[Proved, infrastructure observation] Two Typhon SSH attempts timed out during banner exchange. No remote jobs submitted. Only subsecond local numerical checks ran. Large finite Gaussian script and SLURM file are prepared but unrun. No large interacting job is prepared.

[Proved, file map] RESULT.md is the summary. PROOF.md contains the derivations and explicit qualifications. NUMERICS.md documents reproduction. xx_halfline.py, checks.py and certificate_checks.py produce their matching JSON results. cluster_job.sbatch is unsubmitted. PROGRESS.md records the work. Parent STATE.md and every other file outside work_C105 were left untouched under the task-specific instruction.
