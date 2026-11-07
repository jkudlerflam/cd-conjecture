This run was done by a Claude Opus agent (2026-11-07), not by Codex.

# C145 result: the ferromagnetic strip, the shell, and C138's continuum bounds

Proofs are in [PROOF.md](PROOF.md), the log in [PROGRESS.md](PROGRESS.md). No independent referee pass of the new material. J=1, eta = pi - gamma, Delta = -cos eta, lambda(x) = eta L r_Q(x).

## Bottom line

[Assumed, unresolved] The logarithm is not removed, and the strip -1 < Delta < -1 + delta_0 is not closed unconditionally. What is left is smaller and sharper than before: one shell inequality (G) on the Bethe roots and two one-line continuum statements (W1'), (W3'). All three hold numerically with wide margin down to 1 + Delta = 1e-6.

[Proved, conditional on (G), (W1'), (W3')] On the whole attractive side -1 < J_z/J <= 0 (delta_0 = 1),

    |F_L(n) - L e_bulk(Delta, n/L)| <= 136.61/L,     L(E_(2L) - 2E_L) <= 686 M,

and so on the whole U(1) class J_z/J in [-1, 1], with every field, DM angle and L >= 2, the single constant C = 4078 (from the repulsive side, C142).

## Task 1: the logarithm

[Proved] Route (c), carried out (PROOF Part A). With C142's log-weak (D), (D1) inserted into C138 Theorem B, the core term is already uniform: the variation of g = r/sigma on the deep core K_1 = {lambda >= 4C_D(eta)} satisfies eta Var_(K_1)(g) <= C_G1 <= 13.19 (analytic, via C_D1/C_D <= 4.88; 10.11 numerically). The log comes only from the first-order charge of the shell S = {2 <= lambda < 4C_D(eta)}, which carries about 2(M_1 - M_0)/eta, that is order log(1/eta)/eta, roots. So the shell cannot be counted away or handed to Theorem A, and only second-order control on S is needed.

[Assumed] The inequality left, (G): on S, with M_0 = 2,

    sigma_L >= (3/4) r_Q     and     eta Var_S(r_Q/sigma_L) <= C_GS,

uniformly in eta, L, n. It is weaker than C138's (D), (D1) on S (an integrated derivative bound, nothing pointwise), and asks nothing on K_1 or on the sparse tail.

[Numerical] check_shell*.json, 36 cases with K_0 = {lambda >= 2} nonempty, 1+Delta in {1e-2, 1e-4, 1e-5, 1e-6}, L to 16384: min sigma/r on K_0 is 0.995 (needed 0.75; >= 1.0000 except one mesh-limited near-half-filled case), and eta Var_(K_0)(r/sigma) <= 0.0018 (Theorem A.2 uses C_GS = 1). At 1 + Delta <= 1e-4 every sampled K_0 lies inside the shell, so these test (G) on the whole shell. The comparability half has the margin of (D) itself, eta L sup|sigma - r| <= 0.32 against 0.5.

[Proved, obstruction] Routes (a) and (b) both fail at the same place (PROOF Part D). The second-order charge of the cells at distance u in (eta, 1) has a derivative term whose du/u integral is C_D1 log(1/eta)/lambda^2, even with a log-free C_D1. A continuity bootstrap in eta with hypothesis "sigma >= r/2 on {lambda >= M}" returns only (1/pi) log(1/eta) at the boundary of that set. Closing (G) needs a derivative bound of the form |sigma' - r'| <= (C/(eta L))(1 + 1/dist(y, sparse set)), which is a new statement.

## Task 2: (W1)-(W3)

[Proved] (W2), uniformly in eta and Q (PROOF Theorem B.1). A reflected copy of r_0 is a supersolution: r_0 - r_Q <= r_0(2Q - x) + r_0(2Q + x) on [-Q, Q], hence

    r_Q(x) >= (1 - 4 e^(-2d)) r_0(x),   d = Q - |x|,

so kappa = 0.459 at d = 1 (C138 sampled 0.762). No Wiener-Hopf factorization is needed. Numerics (check_W2*.json) agree, with the only apparent violations shrinking under mesh refinement.

[Proved] |eps_Q'| <= 2a tanh(2Q) r_Q (sharpens C136 (5.2)). Hence (W1) follows from (W1'): v_Q nondecreasing on [0, Q], with C_W = 4pi; and (W3) follows from (W3'): eta Z_Q(Q)^2 >= c_Z tanh(2Q), with C_e = pi/c_Z. The edge-singular parts of r'/r and Z'/Z cancel exactly in v_Q(Q), so both reduced statements are about the edge size of Z_Q and the shape of v_Q, not about a singular profile.

[Assumed] (W1'), (W3'). [Numerical, check_W13.json, 56 cases, 1+Delta 0.5 to 1e-6, Q 0.003 to 3] v_Q monotone in 55 cases (one relative dip of 3.4e-6 in an underresolved case); min eta Z_Q(Q)^2/tanh(2Q) = 0.808 at 1+Delta = 1e-6, with values that depend only on Q/eta in the dilute window, consistent with a limit near 0.8. Used: c_Z = 3/4. Sampled Var(v_Q)/eta <= 7.81, e_bulk''/eta^2 <= 1.24.

The Fourier symbol 1 - T^(w) = 2cosh(gamma w/2) sinh(eta w/2)/sinh(pi w/2) and its Gamma-function factorization would give the edge behaviour of Z_Q at Q = infinity (half-line problem), but a version uniform in Q with explicit constants was not carried out.

## Task 3: referee of C138 Theorem B

[Proved] Correct as stated (PROOF Part C): Lemma 1's Peano constant, the core variation, the tail and straddling-cell charges, the mismatch and the arithmetic of (4.3) (124.87 <= 125, 627.5 <= 628) all check. Remarks: C_W = 8 is sharp at the free point at half filling, not just sampled; (D1) is needed up to the Fermi point whenever K_c reaches it, where C138's numerics stopped short. With the proved kappa(d) and optimized split distance the C138 constant becomes 128.99 (C = 647.4).

## Task 4: assembly

| region | C | source | status |
| --- | --- | --- | --- |
| 0 <= J_z/J <= 1 | 4078 | C142 Corollary I.10 | proved, not refereed |
| -1 + delta <= J_z/J <= 0 | 6 * 10^15 delta^(-25/3) | C136 | proved, refereed |
| -1 | 0 | C129 Section 8 | proved |
| -1 < J_z/J <= 0 | 686 (661.6 with sampled inputs) | Theorem A.2, PROOF E.1 | conditional on (G), (W1'), (W3') |

The seams at 0 and 1 glue as in C142 (same reductions, M >= M_0, finite-matrix continuity at Heisenberg). The strip constant 686 = 5 * 136.61 + pi^2/4, with terms 11.15 (core), 3.29 (Var v), 100.25 (first-order tail, M_0 = 2), 16.76 (straddling cells), 5.17 (mismatch).

## Task 5 and files

check_W2.py, check_W2_refine.py (W2), check_shell.py, check_shell_big.py (shell, outputs check_shell_*.json, check_shell_big_*.json, logs log_shell_*.txt), check_W13.py (W1', W3'), constants_C145.py (all constants, output constants_C145.json), summarize.py (summary_C145.json). Solvers copied from work_C142: bethe_ferro.py, continuum.py, check_attractive.py. All runs local, each under about two minutes; Typhon was not needed.
