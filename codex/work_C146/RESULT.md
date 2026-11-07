This run was done by a Claude Opus agent (2026-11-07), not by Codex.

# C146 result: the edge bound, the velocity variation, and the shell

Proofs are in [PROOF.md](PROOF.md), the log in [PROGRESS.md](PROGRESS.md). No independent referee pass of the new material. Notation of C145: J = 1, eta = pi - gamma, Delta = -cos eta, lambda = eta L r_Q.

## Bottom line

[Proved] Two of the three statements left by C145 are closed, (W3) through (W3') with an explicit constant and (W1) through a replacement of (W1'). [Assumed] The shell inequality (G) is not proved. Its derivative half is reduced to its comparability half up to a log log(1/eta) factor, so the one statement left for a uniform strip constant is still (G), and the one statement left for a log log constant is the one-sided bound sigma_L >= (3/4) r_Q on the shell.

[Proved, conditional on (G) with M_0 = 2, C_GS = 1] On -1 < J_z/J <= 0, every field, DM angle and L >= 2: |F_L(n) - L e_bulk| <= 154.0/L and L(E_(2L) - 2E_L) <= 773 M. On the whole U(1) class J_z/J in [-1, 1] the constant is 4078 (repulsive side, C142).

## Task 1, (W3'): proved

[Proved] Theorem 1. For every eta in (0, pi/2] and Q > 0,

    eta Z_Q(Q)^2 >= pi Q / (4(1 + 16Q/(3 pi^2))).

So (W3') holds with c_Z = 0.373 (clean value 0.309), and e_bulk'' <= C_e eta^2 with C_e = 8.43 (clean 10.16). No Wiener-Hopf factorization and no limit equation are needed. The proof has three steps: W(Q) = int Z_Q has W'(Q) = 2Z_Q(Q)^2 (domain derivative); Z_Q(Q) is nondecreasing in Q (translate the larger interval so the right edges match, then positivity), so W is convex and Z_Q(Q)^2 >= W(Q)/(2Q); and W = <1, (I - T_Q)^(-1) 1> is bounded below variationally with the semicircle trial function and the symbol bound 1 - T^(w) <= eta(|w| + 2/pi).

[Proved, checks] As Q -> infinity the bound tends to 1.453 against the exact pi/2. In the dilute window Q = s eta the limit equation is the Lieb-Liniger equation with c = 1, and the bound Z^2 >= pi s/4 sits a factor 2 below the weak-coupling value pi s/2 [Heuristic for the exact value], which also identifies C145's limit of eta Z^2/tanh 2Q as pi/4 = 0.785. [Numerical, check_W3.json, 161 cases, 1 + Delta from 0.5 to 1e-6, Q from 0.003 to 4] Z_Q(Q) is monotone in Q, dW/dQ = 2Z^2 to 1.5 percent, and the true eta Z^2 exceeds the bound by at least a factor 1.578.

## Task 2, (W1): replaced and proved

[Proved] Theorem 2. For every eta and Q,

    Var_[-Q,Q](v_Q) <= 4a tanh 2Q + 4a I(b) + 12 a tanh(2Q)(Z_Q(Q) - 1) <= V_W := 92.9,

with I(b) <= 5.998, and V_W -> 6 pi^2 = 59.2 as eta -> 0. Theorem A.2 of C145 uses Var(v_Q) only through the core term Var(v_Q)/(6L), so a bound independent of eta suffices and the factor eta of C138's (W1) is not needed. The mechanism is an exact first-order identity obtained by freezing v(x_0) in eps'' - v(x_0) r':

    v' = F_1 - (r(Q)/r)[(v(Q) - v) G + (v(Q) + v) H],   F_1 = a[(R a1')^2 - (R a1)(R a1'')]/(R a1)^2,

where G, H are the restricted resolvent kernels at the two Fermi points. The edge terms are never resolved, only integrated: int_0^Q (G + H) = Z(Q) - 1 exactly. F_1 >= -4a everywhere, and away from the edge F_1 is nonnegative up to an error exponentially small in the distance, from the reflection bound of C145 Theorem B.1 applied to derivatives.

[Assumed, not needed] (W1'), monotonicity of v_Q, is still unproved. [Numerical, check_W1_*.json, 30 cases, 1 + Delta 0.5 to 1e-6] The identity holds to 1e-3 (1.3e-2 in mesh-limited cases), F_1 >= 0 in every resolved case (a Turan-type inequality, not proved), Var(v)/a <= 3.0.

## Task 3, (G): not proved

[Proved] Proposition 3.1. If sigma_L >= c r_Q on the shell S = {2 <= lambda < 4C_D(eta)}, then eta Var_S(r/sigma) <= (1/c)[2 eta log(M_1/M_0) + eta kappa_eta |S|], with eta kappa_eta <= 4.45 (a local bound |sigma'| <= kappa_eta sigma from |T'| <= kappa_eta T) and |S| of order log(C_D(eta)), that is log log(1/eta). With c = 3/4 this gives C_GS(eta) = 17.7 at eta = 0.1 and 25.1 at eta = 1e-12.

[Assumed] What is left, (G_c): sigma_L >= (3/4) r_Q wherever 2 <= eta L r_Q < 4C_D(eta). [Numerical, C145 check_shell*.json] Margin min sigma/r = 0.995 against 0.75 (36 cases, L to 16384, 1 + Delta to 1e-6). For a uniform constant the full (G) is still needed (the log log comes from |S|). [Proved, obstruction] The logarithm in the first-order charge is carried by the sparse cells' first moments mu_j about their roots (each root is exactly the sigma-median of its cell). A local bound on mu_j turns the log into a log log; removing it needs the spacing regularity of the sparse roots (gradient of the background density of order sigma, not sigma/eta), which was not proved.

## Task 4, assembly

| region | C in L(E_2L - 2E_L) <= C M | status |
| --- | --- | --- |
| 0 <= J_z/J <= 1 | 4078 | proved (C142), not refereed |
| -1 + delta <= J_z/J <= 0 | 6 * 10^15 delta^(-25/3) | proved (C136), refereed |
| J_z/J = -1 | 0 | proved (C129) |
| -1 < J_z/J <= 0 | 773 | conditional on (G) only |
| -1 < J_z/J <= 0 | 838 at eta = 0.1 to 867 at eta = 1e-12, order log log(1/eta) | conditional on (G_c) only |

Strip terms (K = 154.0): core 11.15, velocity variation V_W/6 = 15.48 (was 3.29 under (W1')), first-order tail 100.25, straddling cells 16.76, mismatch pi^2 C_e/8 = 10.39 (was 5.17 under c_Z = 3/4). C_strip = 5K + pi^2/4 = 772.5; with C_GS = 5 it is 788.3.

## Task 5 and files

check_W3.py (Theorem 1 and its ingredients, check_W3.json), check_W1.py (Lemma 2.1 and Theorem 2, check_W1_*.json, logs log_W1_*.txt), constants_C146.py (every constant above, constants_C146.json). Solvers copied from work_C145: continuum.py, bethe_ferro.py, check_attractive.py. All runs local, each under fifteen minutes; Typhon was not needed. No Bethe-root computation was rerun; the shell margins are C145's.
