# C149 progress log

Run by a Claude Opus agent (2026-11-07), not by Codex. J=1 throughout. Notation of C145 PROOF Section 1 and C146.

## Step 0, reading
Read specs/C149_comparison_principle.md, specs/_strip_context.md, specs/_Cd_target.md, specs/_bound_background.md, work_C146 (RESULT, PROOF), work_C145 (RESULT, PROOF), work_C142 PROOF Part II, work_C138 (RESULT, PROOF), work_C147 REFEREE, work_C136 RESULT and PROOF Section 1. Copied bethe_ferro.py, continuum.py, check_attractive.py from work_C146.

## Step 1, analysis (first pass)
- Comparison principle. With the Bethe map in theta form, F_j(x) = L p(x_j) - sum_k theta(x_j - x_k) - 2 pi I_j, a maximum-principle argument at the index maximizing x_j - y_j gives F(y) >= F(x) => y >= x for all real configurations, with no ordering hypothesis (theta nonincreasing, p strictly increasing). Uniqueness of the real solution is a corollary.
- Residual identity for quantile configurations: if Z~(y_j) = I_j/L for an increasing Z~ with Z~(A) = -nu/2, Z~(B) = nu/2, then F_j(y)/(2 pi L) = D(y_j) - Phi_j, D the continuum counting residual of rho~ = Z~' on [A,B], Phi_j the quadrature error of theta(y_j - .) against the nodes. The sign part of theta is integrated exactly at a node. Writing Z~ = Z_0 - Psi with Z_0 the continuum counting function of r_Q, D = (I - T)Psi up to end terms.
- Heuristic sizes: Phi_j is about eta^2 rho'/(e^(2 pi lambda) - 1), at most (eta/(2 pi L))|rho'/rho| at first order, so a uniform twist Psi = delta Z_Q with delta of order eta/L gives a sandwich of width of order one root, with no logarithm. The Fermi edge needs separate treatment (rho'/rho is of order 1/eta in the edge layer).
- Conversion: a sandwich of width O(1) roots on the unit neighbourhood of the shell gives comparability once M_0 is large compared to the width, and an integrated derivative bound eta int_S |sigma' - r'|/r of order (width)/M_0, through int |T'(. - x_k) - T'(. - y_k)| <= ||T''||_1 |x_k - y_k|.

## Step 2, numerics on the residual (explore_half.py, explore_disp.py, scan_phi.py)
- Half filling, explicit trial y_j = (2/b) atanh(tan(pi I_j/L)) (quantiles of r_0, D = 0 exactly): Phi_j = -F_j(y)/(2 pi L) agrees with eta^2 rho'/(e^(2 pi lambda) - 1) to 2e-3 eta/L at 1+Delta = 1e-2 and 2e-5 eta/L at 1e-6, and max L|Phi|/eta -> b/(2 pi). True displacement L|Z_0(x_j) - I_j/L| <= 0.263 roots, about b/4.
- The old graded Nystrom solver is mesh-limited near Delta = -1 (r_Q(0) off by 5 to 10 percent at 1e-6, Q = 8), because interpolation errors are amplified by B. Replaced by fftsolver.py: uniform grid h = eta/10, Gregory end weights of order 7, FFT Toeplitz matvec, CG on the symmetrized system. m(Q) converges to 1e-12, r_Q(0) -> 1/(2 gamma) at large Q.

## Step 3, the simplest barrier
- Uniform barriers (constant shift of the quantum numbers, Psi = kappa) suffice numerically. The general identity F_j(y)/(2 pi L) = [Z_0(y_j) - I_j/L] - T*(N_0 - N_y)(y_j) holds for every configuration, so no room condition beyond existence of the top node (kappa0 < 1/2 + eta(L - 2n)/(2 pi)) is needed.
- Hypothesis (R_kappa): quantiles of Z_0 at levels I_j/L +- kappa0/L are super/subsolutions. rk_test.py: holds with kappa0 = 0.4 at 1+Delta = 0.5, 0.1, 1e-2, 1e-4 for L to 4096, all fillings; worst case half filling, margin 0.093 eta/L at 1e-4 (asymptotically (2/pi)(kappa0 - b/4)).

## Step 4, conversion (the main gain)
- The sandwich gives the cumulative bound ||N_x - N_0|| <= (kappa0 + 1/2)/L relative to the continuum at m(Q) = nu (not at Q_L). With the factorization T = T_A * T_B (both positive unimodal, Fourier identity), this gives log-free (D) and (D1) on all of R: eta L|sigma - s_nu| <= (2 kappa0 + 1) eta cot(eta)/pi, eta^2 L |sigma' - s_nu'| <= (2 kappa0+1) 2 eta^2 T_A(0) T_B(0) <= (2 kappa0 + 1) 8/pi^2. Checked on Bethe roots (check_density.py).
- Energy identity with Q_nu in place of Q_L (no mismatch term, a boundary layer of size O(a/L) instead). C138 Theorem B then runs with the log-free constants: K = 136.4, C_strip = 685 at kappa0 = 0.4 (constants_C149.py).
- Coordinator message: C150 needs only sigma_L >= c_* r_(Q_L) on S'. Delivered when Q_L <= Q_nu (c_* >= 0.714); Q_L < Q_nu in every sampled case; the transfer for Q_L > Q_nu is proved only for a dense Fermi edge.

## Step 5, final numerics and write-up
- rk_test.py at 1 + Delta = 1e-6 (L = 16 to 8192, 49 cases): (R_kappa) at kappa0 = 0.4 holds with margin >= 0.088 eta/L (half filling); kappa0 = 0.25 fails, kappa0 = 0.3 holds with 0.020. In all, 173 cases from 1 + Delta = 0.5 to 1e-6.
- check_density.py (30 cases at 1e-2, 1e-4, 1e-6): cumulative error <= 0.758/L, eta L|sigma - s_nu| <= 0.317, eta^2 L|sigma' - s_nu'| <= 0.217, all inside the conditional bounds; Q_L < Q_nu always.
- PROOF.md and RESULT.md written. Strip reduced to (R_kappa); conditional constant 685 (kappa0 = 0.4) on the whole attractive side. Typhon not used.
