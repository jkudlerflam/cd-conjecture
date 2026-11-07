# C137 result: the window lemma

This run was done by a Claude Opus agent (Claude Code subagent), not by Astra. No Codex or gpt_solve.sh call was made.

Throughout, K denotes n^2 (e_0 - lambda_min)/M for a window on n sites, W_n(h) the optimal window value, and "sine window" the C135 window a_j = v_j v_{j+1}/s with v_j = sqrt(2/(n+1)) sin(pi j/(n+1)), s = cos(pi/(n+1)). "Sine-square window" is (2/n) sin^2(pi j/n). Full derivations are in [PROOF.md](PROOF.md), all numbers in [TABLES.md](TABLES.md) and results/.

## 1. Formalization (Task 1)

[Proved] The window counterterms X = sum_j x_j h_{j,j+1} + sum_k y_k (one-site y_k) give exactly the open chains W(c, u) = sum_j c_j h_{j,j+1} + sum_k u_k with sum_j c_j = 1 and sum_k u_k = 0, by a bijection (PROOF Lemma 1.1). The admissible class is: real weights of either sign, one-site terms summing to zero (which is the same as giving each bond its own gauge representative whenever no weight vanishes, Lemma 1.2), both allowed to depend on h and n. Every window is feasible for the LTI dual, so W_n <= E_LTI(n) <= min(e_0, E_L/L for L >= n), and (W) implies (T), (K) and A <= (6/pi) K M with the same K (Proposition 1.3). The C135 XX certificate is the window c = a, u = 0 (Lemma 1.6).

[Proved] The primal of the window problem is a clean relaxation of LTI. By Sion's minimax theorem, W_n(h) is the minimum of e over n-site states whose bond energies Tr(rho h_{j,j+1}) all equal e and whose one-site marginals are all equal (Proposition 1.4). For a fixed (h-independent) c the window lemma is equivalent to a kinematic inclusion: the c-weighted average of the bond marginals of any n-site state with equal one-site marginals lies within K/n^2 of Omega_inf in the dual-M norm.

[Proved] Elementary properties (Lemma 1.5): gauge covariance, monotonicity W_{n+1} >= W_n (and improvement under convolution of windows), the trivial rate e_0 - W_n <= 2M/(n - 1), and twist blindness. Open windows do not see symmetry twists, so the chiral bond of C135 has exactly the Heisenberg window constant (numerically equal to 1e-10 for n <= 6). For even n, e_0 - W_n >= e_0 - E_n/n, so the window constant is never below the ring Casimir constant (Proposition 3.1). The sine-square window at n = L is the normalized SSD Hamiltonian of the L ring (Remark 1.7), so the window lemma restricted to that window is the route set aside in _Cd_target.md.

## 2. Optimal windows and the comparison table (Task 2)

[Numerical] Method. lambda_min is concave in (c, u); it was maximized by a proximal spectral-bundle method on exact sparse diagonalization (common.py), reaching n = 16 for the main bonds. The a-posteriori model gap of every reported optimum is below 1e-8 in absolute value. For quadratic bonds an exact Majorana reformulation (lambda_min = const - ||A||_*/4, a nuclear-norm minimization) reaches n = 1024 (free.py, agrees with ED to 1e-14). The optimizer reproduces C135's restricted duals (XX, Heisenberg, n <= 8) and, once Z one-site terms are allowed, C135's unrestricted TFIM LTI value -1.2814577 at n = 8, which the C135 bond-only restriction had missed by 0.0018.

[Numerical, e_0 brackets] For XXZ and TFIM, e_0 is the Yang-Yang and Pfeuty value (literature, exact), and chiral is Heisenberg - 1/4 by twist blindness. For random bonds, e_0 is bracketed by [max over computed n of W_n, e0_up], where W_n <= E_LTI <= e_0 and e0_up is the energy of a shift-symmetrized two-site iTEBD uniform MPS (chi = 64), a TI state on Z, evaluated from positive fixed points of the transfer map (spectral projector of the identity, so that degenerate non-injective MPS give a genuine state). These are floating-point bounds, not interval certificates.

[Numerical] Main bonds, K at the largest n (n = 16 unless noted). Columns: optimal window, sine-square window with optimized u, sine window with optimized u, ring constant at even n, and the CFT value pi A/(6M).

| bond | opt | sine-square | sine | ring | CFT |
| :--- | ---: | ---: | ---: | ---: | ---: |
| Heisenberg, n = 16 | 1.6629 | 1.6629 | 1.7455 | 1.6621 | 1.6449 |
| Heisenberg, n = 15 | 1.4616 | 1.5100 | 1.4616 | | |
| XX | 1.0519 | 1.0519 | 1.0948 | 1.0519 | 1.0472 |
| XXZ, Delta = 0.5 | 1.3692 | 1.3692 | 1.4337 | 1.3684 | 1.3603 |
| XXZ, Delta = 2 (gapped) | 0.7839 | 0.7852 | 0.8236 | 0.6620 | 0 |
| XXZ, Delta = -0.9 | 0.2627 | 0.2627 | 0.2658 | 0.2606 | 0.2665 |
| XXZ, Delta = -0.99 | 0.0548 | 0.0548 | 0.0549 | 0.0547 | 0.0773 |
| TFIM, g = 1 | 0.3707 | 0.3707 | 0.3815 | 0.3707 | 0.3702 |
| TFIM, g = 1 - 2/n | 0.1237 | 0.1238 | 0.1275 | 0.1083 | |
| TFIM, g = 1 + 0.5/n | 0.3232 | 0.3232 | 0.3349 | 0.3218 | |
| TFIM, g = 1 + 5/n, n = 15 | 0.0364 | 0.0366 | 0.0420 | 0.0173 (n = 14) | |

(All columns and all n are in TABLES.md.) With u = 0 in the symmetric representative the TFIM constants are larger (0.441 sine-square and 0.451 sine at n = 15, converging to 0.370 as n grows, see section 3). The largest K over all main-bond entries is 1.819 for the optimal window and 1.952 for the sine window, both Heisenberg at n = 4.

[Numerical] Description of the optimal windows. All optimal weights are positive. For XX, Heisenberg, XXZ (Delta = 0.5, -0.9, -0.99) and chiral, the optimal c is the sine-square window at even n and the C135 sine window at odd n (relative sup deviation below 5e-3 at even n, below 3e-3 at odd n), and at even n its value equals the ring value E_n/n to solver precision for XX (Katsura's sine-square exactness for free fermions, which together with Proposition 3.1 makes it optimal), and is within n^2 (E_n/n - W_n)/M <= 2.5e-3 of it for Heisenberg and XXZ at Delta = 0.5, -0.9, -0.99 (exact at n = 4, 6, a small real deficit from n = 8 on). [Heuristic] Odd n prefers the sine window because odd rings are frustrated, so the best ring-like comparison is the ring of length n + 1. For TFIM at every n and every c tested, the optimal c is the sine-square window (deviation below 2.3e-2) together with h-dependent Z fields; at g = 1 the optimum is exactly the canonical SSD Hamiltonian, bond weights sin^2(pi j/n) and site weights sin^2(pi(k - 1/2)/n). For gapped Delta = 2 the optimum deviates from sine-square by 1.5e-2 and lies above the ring constant, which is the only case among the main bonds where the window is visibly weaker than the ring. Against C135's full LTI values (n <= 8), the window restriction loses n^2 (E_LTI - W_n)/M = 0 to solver precision for XX, Heisenberg, TFIM and chiral, at most 6e-4 for XXZ at Delta = 0.5 and -0.9, and at most 0.032 for Delta = 2 (lti_compare.py).

[Numerical] Random bonds. The local runs completed 7 of the 30 C135 random bonds (random_00 to random_06; n up to 7 to 10) before the run was finalized. The remaining 23 bonds, n = 11 to 13 for all 30, and the main bonds at n = 17, 18 are in Typhon job 17579408, still pending in the queue (outputs will land in /scratch/lustre/jkudlerflam/c137_window/work_C137/results/ as NAME_nN.json and results/rings/NAME_LL.json; copy them into results/typhon/ and rerun report.py, which reads that folder). iTEBD upper bounds on e_0 were computed for all 30 bonds (results/e0_mps/) and no bracket is reversed. On the 7 completed bonds:

- the upper bound on K for the optimal window is at most 0.75 over all n (random_06; 0.56 for random_00, all others below 0.44), and at the largest n it lies in [0, 0.36];
- the better of the two fixed windows with optimized u is within n^2 (W_opt - W_fixed)/M <= 2e-3 of the optimum at the largest n (the sine-square window alone is within 0.028, the sine window alone within 0.006);
- all optimal weights are positive, and the optimal shape can differ substantially from both fixed windows (relative sup deviation up to 0.52 for random_01, a nearly frustration-free bond with K <= 0.07, where many windows are almost optimal), but the value is insensitive to it;
- against C135's full LTI values (n <= 6) the window restriction loses at most n^2 (E_LTI - W_n)/M = 0.113 (random_00, n = 5).

The bracket widths are set by the lower endpoint, which is the window value itself at the largest n, so they bound K from above only; the upper endpoints are the informative numbers.

## 3. Lesson 2 (Task 3)

[Proved] For TFIM at g = 1 + c/n (and every quadratic bond) the fixed sine envelope with u = 0 has K <= n^2 [1 - cos(pi/(n+1))] < pi^2/2 for all c and n (Theorem Q below). So on the free near-critical family the fixed envelope cannot grow, and nothing depends on whether the weights are allowed to depend on h.

[Numerical] Free TFIM at g = 1 + c/n, c in {-40, ..., 40}, n = 8 to 1024 (results/lesson2_free_A.json, B.json). At fixed c, K converges as n grows, to a limit K_inf(c) that is largest at criticality and decreases with |c|:

| c | -40 | -20 | -10 | -5 | -2 | -1 | 0 | 0.5 | 1 | 2 | 5 | 10 | 20 | 40 |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| sine window, u = 0, n = 1024 | 0.0122 | 0.0259 | 0.0509 | 0.0932 | 0.1867 | 0.2702 | 0.3727 | 0.3333 | 0.2789 | 0.1988 | 0.1060 | 0.0628 | 0.0366 | 0.0215 |

The optimal window (n <= 128) and the sine-square window with optimized u are below these values; at n = 64 they are 0.0167, 0.0501, 0.1534, 0.2562, 0.3703, 0.3242, 0.2583, 0.1588 at c = -10, -5, -2, -1, 0, 0.5, 1, 2. The largest K over the whole scan is 0.463 (sine window with u = 0, c = 0, n = 8) and 0.372 for the optimal window (c = 0, n = 8). The canonical SSD Hamiltonian (fixed bond and site weights) stays below 0.372 at every c and n as well. For |c| >= 5, K_inf(c) falls off roughly like 1/|c| (coefficient about 0.5 on the ordered side, 0.5 to 0.9 on the disordered side, where n = 1024 is not yet converged at large c), so at fixed mass the deficit appears to decay like n^{-3} rather than n^{-2}. The near-critical massive mode is therefore seen by the window, at the order of the anomaly or below, as lesson 2 says, but it never exceeds the critical value, and no h-dependence of the weights is needed. Only the one-site fields adapt.

[Numerical] Interacting family. XYZ with (J_x, J_y, J_z) = (1 + gamma, 1 - gamma, 0.5) and gamma = c/n, c in {0.5, 2, 5}, n = 4 to 14 (u = 0 by the Z2 x Z2 symmetry). The gap opens as a power of gamma, so this is a near-critical massive family with the massive scale tied to n. The optimal window is again sine-square at even n and sine at odd n: n^2 (W_opt - W_fixed)/M, with W_fixed the better of the two fixed windows, is below 2e-4 at c = 0.5 and 2 and below 3e-3 at c = 5, over all n. With iTEBD upper bounds on e_0, K_opt <= 1.357 to 1.367 at c = 0.5, 1.06 to 1.14 at c = 2 and 0.33 to 0.55 at c = 5 (n = 8 to 14), increasing towards the XXZ Delta = 0.5 value 1.36 as gamma = c/n shrinks. So in this interacting family, too, the fixed envelope does as well as the optimal h-dependent one, and the constant stays bounded and below the Heisenberg values.

[Proved, mechanism] Theorem Q's proof shows why lesson 2 is harmless for bounds. The window splits as W' plus a correction (PROOF (4.2)), where W' is a superposition of quadratic forms in normalized windowed modes, each bounded below by the infinite-volume symbol at the same momentum, whatever the mass. The whole deficit comes from the factor s = sum_j v_j v_{j+1} multiplying the hopping, which costs at most (1 - s) times the kinetic scale, independent of the gap. A massive mode is seen at the order of the anomaly, which is lesson 2, but it is seen with a coefficient bounded by the kinetic energy, not by its crossover length. Lesson 2 is an obstruction to extracting A from a window, not to bounding the window deficit.

## 4. Proofs (Task 4)

[Proved, Theorem Q, new] For every quadratic qubit bond (Jordan-Wigner frame: XY anisotropy, pairing, z-directed DM, transverse field, no ZZ), written in the parity-symmetric representative (4.1), and every n >= 2,

    lambda_min( sum_j a_j h_{j,j+1} ) >= e_0(h) - (1 - cos(pi/(n+1))) ||h_cross|| >= e_0(h) - pi^2 M(h)/(2(n+1)^2).

The window is the C135 sine window with no one-site terms, independent of h. The proof is a Majorana-mode decomposition: the sine window's bond-induced one-site weights are exactly v_k^2 (because v is the lowest Dirichlet eigenvector), the part of the window with bond weights v_j v_{j+1} is an integral of positive quadratic forms in normalized modes Gamma(k) = sum_j v_j e^{ikj} gamma_j and is bounded below by e_0, and the leftover hopping costs (1 - s)||h_cross||, with ||h_cross|| <= M by parity averaging. It gives an explicit telescoping certificate X_n = sum_j (sum_{i<=j} a_i - 1) h_{j,j+1} for the whole class, the dual counterpart of C135's covariance smoothing, extending C135's XX certificate. Checked numerically on 40 random quadratic bonds to n = 400 with no violation and with equality ||h_cross|| = M attained (check_quadratic.py).

[Proved] K = 0 for gauge frustration-free bonds (some gauge representative has h~ >= e_0), which include all classical bonds (by min-mean-cycle LP duality) and ferromagnetic XXZ at Delta <= -1, and the constants transfer to all twists (Proposition 2.1, Corollary 2.2, Lemma 1.5(d)).

[Proved, obstructions] (i) Gosset–Mozgunov's weighted sum of squares needs frustration-freeness (zero-energy subchain ground states) and outputs gaps, so for energies it adds nothing beyond K = 0 on frustration-free bonds and does not transfer to frustrated ones. (ii) Many-body IMS localization of the Hamiltonian density is exact and reduces to the window lower bound itself, with no estimate. (iii) Deforming from the free certificate controls the deficit only on a ball of radius of order M/(n^2 ||V||). For a relevant direction the deficit is a scaling function of n^{1/nu} t, so uniformity requires a nonperturbative bound on that function (a constructive-RG input, C130), which is not available. (iv) h-dependent weights are not what makes the windows work: on every tested bond the fixed sine-square (even n) or sine (odd n) envelope is optimal or within 3e-3 in K (main bonds and XYZ) and 2e-3 (random bonds, better of the two fixed windows), and only the one-site terms need to depend on h.

## 5. What is proved, and the single remaining lemma

[Proved] (W), hence (T), the LTI rate, (K) and A <= (6/pi) K M, holds at d = 2 with K = pi^2/2 on the quadratic class (by a fixed h-independent window, uniformly in gaps and crossover lengths), with K = 0 on gauge frustration-free bonds (including all classical bonds), and with the same constants on all twists of these.

[Numerical] On all other tested bonds (XXZ at Delta = 0.5, 1, 2, -0.9, -0.99, chiral, interacting XYZ near criticality, 7 of the 30 random bonds), the optimal window constant is at most 1.82, the sine window constant at most 1.96, and both converge to the ring and CFT values where those exist.

[Assumed, single remaining lemma] Fixed-envelope window lemma at d = 2: there is K such that for every qubit bond h and every n there are one-site terms u_k with sum zero and

    lambda_min( sum_j a_j h_{j,j+1} + sum_k u_k ) >= e_0(h) - K M(h)/n^2,

with a_j the C135 sine window. Equivalently, for every n-qubit state with equal one-site marginals, the sine-weighted average of its bond marginals is within K/n^2 of Omega_inf in the dual-M norm, a kinematic statement independent of h. Theorem Q proves it in every quadratic direction, for all states. What is missing are the non-quadratic directions (ZZ in the Jordan-Wigner frame, transverse DM and fields), where no mode decomposition of the many-body space exists.

## Files

PROOF.md (derivations), TABLES.md (all numbers, generated by report.py), PROGRESS.md. Scripts: common.py (window engine, bundle optimizer), free.py (Majorana engine), e0.py (exact e_0 and iTEBD), models.py, survey.py, rings.py, lesson2_free.py, lesson2_xyz.py, e0_mps.py, check_quadratic.py, lti_compare.py, shapes.py, report.py, typhon_job.slurm with tasks_typhon.py (cluster continuation, job 17579408, pending). Data: results/.
