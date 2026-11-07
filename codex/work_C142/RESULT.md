This run was done by a Claude Opus agent (2026-11-07), not by Codex.

# C142 result: the two U(1) strips

Proofs are in [PROOF.md](PROOF.md), the log in [PROGRESS.md](PROGRESS.md). No referee pass yet. J=1.

## Strip (L*), 63/64 < Delta <= 1: closed

[Proved] (L*) holds along the branch at every 0 < gamma <= pi/2, with margin. The 6 percent shortfall in C140 came from one tail integral. C140 Lemma 3.3 bounded int_Q^inf r0 by (2/b) r0(Q), and the sharp constant is pi/(2b), since atan(t)(1+t^2)/t is increasing on (0,1] with value pi/2 at t=1 (PROOF Lemma I.1). So

    L r0(x_n) >= b/(2(pi - gamma)) >= b/(2pi) = 0.1592 b   (needed 0.1324 b, ratio >= 1.202),

and C140 Theorem 3.4 holds with the uniform margin sigma_L >= theta_1 r0 on [x_1, x_n], theta_1 = 1 - 2(2-ln2)(pi-gamma)/pi^2 >= 0.168 (PROOF Theorem I.3). The continuation argument then runs on all of 0 <= Delta < 1, written out in full (the closedness step and the identification facts that C140's referee asked for, PROOF Theorem I.4). Also proved there: the negative mass of sigma_L beyond the outermost root is at most (pi - 2gamma)/(2pi L) per side, half of C140's two-sided 1/L.

[Proved] C139 Section 8.2 (Theorem C) redone line by line with the proved sawtooth sup|Cum R| <= (1 - gamma/pi)/L <= 1/L (PROOF I.3, Lemmas I.5 to I.8). The constants are L||G|| <= 45/4, |sigma - r_Q| <= (2g/(pi gamma))/L, |sigma' - r_Q'| <= 4M_1 g/L, core <= 46.96/L, ends <= 765.50/L, curvature <= 2.47/L. The empty-core case is handled by the first-order bound, so only L >= 29 is needed. Result (Theorem I.9 and Corollary I.10):

    |F_L(n) - L e_bulk(Delta, n/L)| <= 841/L       (0 <= Delta <= 1, every L >= 2, 0 <= n <= L/2),
    L(E_(2L) - 2E_L) <= 4078 M                      (0 <= J_z/J <= 1, every field, DM angle, L >= 2).

4078 = 5 * 815 + pi^2/4 uses analytic majorants of each term. The exact gamma-dependent expression evaluated on a 200001-point grid gives 811.30 and 4058.96 (constants_rep.json). Heisenberg is included by finite-matrix continuity at fixed L. This replaces C139's 345600 on [0, 63/64] and C140's unchecked 19935 + pi^2/4. The constant is above C139's conditional 2215 because the proved sawtooth is up to 1/L, not 1/(2L).

[Numerical] check_repulsive.py, 160 cases (Delta in {0.5, 0.9, 0.99, 0.999}, L from 4 to 256, seven fillings each). Every new inequality holds. The extreme values against the proved bounds are as follows: L r0(x_n)/b >= 0.442 (proved >= 0.159); min sigma/r0 on [x_1, x_n] >= 0.974 (proved >= 0.168); negative mass zero to grid accuracy; L sup|Cum R| <= 0.5000 (proved <= 1); L||G|| <= 0.532 (proved <= 11.25); L sup|sigma - r_Q| <= 0.92; L|F - L e_bulk| <= 0.91 (proved <= 815). ED (ed_check.py) matches Bethe to 2.3e-14 at Delta = 0.99, 0.999 for L <= 12, and the ED defect L(F_2L(2n) - 2F_L(n))/M_0 for L <= 6 peaks at 4.000 (Delta = 1, L = 2).

## Strip (D), -1 < Delta < -1 + delta_0: not closed, reduced to a logarithm

[Proved] (D) and (D1) hold up to a factor log(1/eta), pointwise on the whole occupied interval, uniformly in filling and Fermi boundary (PROOF Theorem II.6):

    eta L |sigma_L - r_(Q_L)|     <= C_D(eta)  <= (4/pi)(log(1/eta) + pi/2),
    eta^2 L |sigma_L' - r_(Q_L)'| <= C_D1(eta) <= 3.61 + (12/pi^2)(log(1/eta) + pi/2)   (0.83 in place of 3.61 as eta -> 0).

The mechanism has three steps. First, sigma_L - r_Q = int calR_Q(x, y) dR_L(y), with calR_Q = (I - T_Q)^(-1) T_Q the restricted dressing kernel. Second, the midpoint-cell oscillation bound |int f dR_L| <= Var(f)/L holds because z_L is monotone on the attractive side. Third, a comparison through the outside channel gives Var_y calR_Q(x, .) <= 2(1+q) calR(0) uniformly in Q, where calR(0) <= (log(1/eta) + pi/2)/(pi eta) is the peak of the full-line kernel. This improves C136's 1/(eta^2 L) (which is 1061/(eta L) at 1+Delta = 1e-6) to 9.97/(eta L) there.

[Assumed, the remaining inequalities] Two inequalities are left.

1. (D) without the logarithm, i.e. eta L sup|sigma_L - r_Q| <= C_D with C_D independent of eta, on the core where eta L r_Q >= 4C_D. The log comes from cells at distance between eta and 1 from x, where calR(x - y) is about (1/(pi eta)) log(1/|x-y|) and the oscillation bound charges first order. A second-order bound there needs sigma_L comparable to r_Q, and a controlled derivative error, on the unit neighbourhood of x. Theorem II.6 gives that only once eta L r >= C log(1/eta), and the bootstrap does not close on the shell C_0 <= eta L r <= C log(1/eta). Numerical margin: the true eta L sup|sigma - r_Q| is at most 0.32 (C138, L to 8192, 1+Delta to 1e-6; 0.317 in the 95 cases of check_attractive_*.json here), against the proved 3.87, 5.48, 7.00, 9.97 at 1+Delta = 1e-2, 1e-3, 1e-4, 1e-6.
2. (W1)-(W3) of C138, which are continuum statements, not addressed here except that the positivity method gives (W2) only at distance of order log(1/eta) from the Fermi point. Numerical margins (C138): Var(v_Q)/eta <= 8.0, e_bulk''/eta^2 <= 4/pi, min r_Q/r_0 at distance >= 1 is >= 0.762.

[Proved, conditional on (W1)-(W3) only] With Theorem II.6 inserted into C138 Theorem B (whose proof uses (D), (D1) only on its core), |F_L(n) - L e_bulk| <= K_*(eta)/L with K_*(eta) = 193.9 C_D(eta) + O(1), which is about 247 log(1/eta). This is length-uniform at each eta but still diverges, and its minimum with C138 Theorem A ((5pi/4) eta L) is also of order log(1/eta). C138 Theorem B itself has not been refereed.

[Numerical] check_attractive_*.json: see the table below.

95 cases: 1+Delta in {0.5, 0.1, 1e-2, 1e-4, 1e-6}, L in {16, 64, 256, 1024}, n from 2 to L/2 (L/2 is even half filling, Q_L = infinity). All pass.

| 1+Delta | eta | cases | max eta L sup|sigma - r_Q| | proved C_D(eta) | max eta^2 L sup|sigma' - r_Q'| | proved C_D1(eta) |
| --- | --- | --- | --- | --- | --- | --- |
| 0.5 | 1.047 | 19 | 0.037 | 0.67 | 0.085 | 1.46 |
| 0.1 | 0.451 | 19 | 0.123 | 2.08 | 0.183 | 3.02 |
| 1e-2 | 0.1415 | 19 | 0.233 | 3.87 | 0.216 | 4.64 |
| 1e-4 | 0.01414 | 19 | 0.308 | 7.00 | 0.213 | 7.53 |
| 1e-6 | 0.001414 | 19 | 0.317 | 9.97 | 0.207 | 10.35 |

The true values saturate near 1/pi while the proved ones grow like (4/pi) log(1/eta), which is the logarithm that remains. The derivative column is evaluated away from the Fermi point (within eta/4 of Q_L the mesh derivative of r_Q is not resolved, as in C138). calR(0) by quadrature matches the inverse Fourier transform to 1e-8, and lies below its closed-form bound (1762.6 against 1830.3 at 1+Delta = 1e-6).

## Assembly on J_z/J in [-1, 1]

| region | constant C in L(E_2L - 2E_L) <= C M | source | status |
| --- | --- | --- | --- |
| [0, 1] | 4078 | PROOF Corollary I.10 (C139 Lemmas 3.1-3.5, C140 Lemmas 3.1-3.4 refereed; Lemmas I.1-I.8 new) | proved, not refereed |
| [-1+delta, 0] | 6e15 delta^(-25/3) | C136 | proved, refereed |
| -1 | 0 | C129 Section 8 | proved |
| (-1, -1+delta) | 628 conditional on (D), (D1), (W1)-(W3); about 1235 log(1/eta) conditional on (W1)-(W3) | C138 Theorem B, PROOF II.4 | open |

The seams are as follows. At J_z/J = 0, both compact theorems hold on closed intervals containing Delta = 0 with the same C129 reductions and the same M >= M_0, so on [-1+delta, 1] the single constant is max(4078, 6e15 delta^(-25/3)) = 6e15 delta^(-25/3). At J_z/J = 1, the bound follows by continuity at fixed L. If (D) and (W1)-(W3) are granted, the whole class has the single constant 4078 (C138's 628 is smaller).

## Files

constants_rep.py (all Part I constants as functions of gamma, output constants_rep.json), check_repulsive.py (output check_repulsive.json, check_repulsive_summary.json), check_attractive.py (outputs check_attractive_<eps>.json), ed_check.py (output ed_check.json). Solvers were copied from C140: bethe_scan.py, check_c129.py, continuum_rep.py, bethe_ferro.py, continuum.py. Typhon was not used, since each run took minutes locally.
