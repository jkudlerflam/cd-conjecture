# C129 handoff

2026-11-07. Scope is the U(1) spin-1/2 doubling defect, all lengths, fields, and DM angles. The full explicit C_U1 theorem is unresolved. This local memory replaces a project-level update because this task prohibits outside writes.

Read RESULT.md first, then PROOF.md for derivations and the literature audit. CHECKS.md gives commands and the actual numerical scope.

[Proved] Field reduction works sectorwise. The DM twist is Ltheta. A minimizing flux is Lpi modulo 2pi for every sector, proved by the configuration-basis diamagnetic inequality. The exact twist upper bound is JL[1-cos(pi/L)]/2. The physical doubling reduction costs at most pi^2 M(B)/(4L). All these statements include odd L and L=2.

[Proved] M(B) is half the range of J_z/4 +/- B/2 and -J_z/4 +/- J/2. At B=0 it is J/2 for |Delta|<=1 and J(1+|Delta|)/4 otherwise. M(B)>=M_0. Do not erase the field from M when using the field reduction.

[Proved] Original defect is zero for J_z<=-J, nonpositive for J=0, and zero in the saturated region |B|>=J+J_z with J_z>=-J. The free class J_z=0 has explicit original constant 8-4sqrt(2)+pi^2/4.

[Proved] New elementary partial sector bounds, for Delta>=-1, are Q<=2pi^2 n^3/L and Q<=(Delta+1)L^2, where Q=L(F_(2L)(2n)-2F_L(n))/M_0. They cover n^3<=L and the ferromagnetic crossover band (Delta+1)L^2<=A. They do not cover the whole low-density or ferromagnetic neighborhood.

[Proved, literature audit] TOOLBOX overstates Kozlowski 1508.05741. Proposition 3.4 has the finite-boundary expansion. Theorem 3.5 is only half-filled convergence for -1<Delta<=1. Corollary 4.3 already contains the 1/L energy term below half filling. Fixed-Delta density uniformity of the leading estimates is stated, Delta-uniform constants are not. The paper is formulated on even lengths. PROOF Section 12 identifies constants and a printed strip-width issue for negative Delta. The original Yang-Yang PDF was inaccessible, so its Delta range is checked through Kozlowski's restatement, not a direct original-text reading.

[Numerical] 250 Bethe-versus-ED sector checks through L=10 pass, error at most 7.24e-14, paired Bethe rings through 20. Max sampled Q=4 at two-site Heisenberg. Twist tests through L=12 and 144 field/DM checks pass. 72 new variational checks pass. JSON files retain outputs. No independent referee pass.

[Assumed, pending] run_typhon.slurm and bethe_scan.py prepare fixed-Delta and crossover scans through L=4096, paired length 8192. Not submitted or run. Large-solver convergence and cancellation must be inspected from output diagnostics. No cluster connection attempted.
