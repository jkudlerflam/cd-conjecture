# C138 result

This run was done by a Claude Opus agent, not by Astra. 2026-11-07. Details and proofs are in [PROOF.md](PROOF.md); the log is [PROGRESS.md](PROGRESS.md). No independent referee pass has been made.

## Bottom line

[Assumed, unresolved] (X) is not proved in the window 1 << eta L << 1/eta, where eta = pi - gamma and epsilon = 1+Delta = 2 sin^2(eta/2). It is reduced, with an explicit constant, to one discrete estimate (D) on the Bethe roots and four continuum estimates (D1), (W1)-(W3). All five hold numerically with margin down to epsilon = 10^-6. (X) itself holds numerically with K_* of order eta.

[Proved] Two unconditional results.

1. Task 1. The nondegenerate limit near gamma=pi is reached by dividing the equations by eta, not by rescaling rapidities. The normalized kernel symbol [1-T^(w)]/eta is at least 2/pi for all w and tends to w coth(pi w/2). In the step-free form L p(x_j) = sum_(k!=j) phi(x_j-x_k), phi(y)=2 atan(tan eta/tanh y), the eta->0 limit at fixed L is the strictly convex log-gas L tanh x_j = 2 sum coth(x_j-x_k). Its exact sum rule sum sech^2 x_j = n(L-n)/(L-1) reproduces the degenerate first-order energy -epsilon n(L-n)/(L-1) of C129 Section 10. The crossover parameter is eta L, the ratio of kernel width to root spacing.
2. Theorem A. For every -1<Delta<=0, L>=2 and 0<=n<=L/2,

       |F_L(n) - J L e_bulk(Delta,n/L)| <= (5pi/4) eta J.

   This gives (X) on eta L <= A_0 with K_* = (5pi/4)A_0, and the sector doubling constant (15pi/2)A_0 in units of M_0, a Bethe-ansatz analogue of the C129 band that grows linearly instead of quadratically in eta L. It does not reach the window.

[Proved, conditional] Theorem B. If (D), (D1), (W1)-(W3) hold with constants C_D, C_D1, C_W, kappa, C_e, then for every -1<Delta<=0, L>=2 and 0<=n<=L/2,

    |F_L(n) - J L e_bulk(Delta,n/L)| <= K_* J/L,
    K_* <= (2pi/(9C_D))[C_D + C_D1(1+1/kappa)] + pi C_W/12
           + 16 pi C_D (1+2/kappa) + 16pi/3 + pi^2 C_e/8.

No length threshold is needed (A=0). With C_D=1/2, C_D1=2, kappa=7/10, C_W=8, C_e=13/10, each above the largest sampled value ((D1) was not sampled within a quarter kernel width of the Fermi point), K_* <= 125, and the attractive U(1) class J_z/J in [-1,0], every field, DM angle and L>=2, would close with C = 5K_* + pi^2/4 <= 628. This number is conditional.

## The exact missing estimate

[Assumed] (D): with r_Q the restricted continuum density at the counting boundary Q=Q_L (r_0 = 1/(2 gamma cosh(pi x/gamma)) at even half filling),

    |sigma_L(x) - r_(Q_L)(x)| <= C_D/(eta L)   wherever L eta r_(Q_L)(x) >= 4 C_D,

uniformly in Delta in (-1,0], L and n. Here sigma_L = a1 + L^(-1) sum_k T(x-x_k) is the counting density of the Bethe roots. The derivative version (D1) is |sigma_L' - r_Q'| <= C_D1/(eta^2 L) on the same set. C136's argument gives these only with eta L replaced by eta^2 L, which is why its constant diverges at every length, not only in the window: with the C136 density bound, the first-order cost of the tail beyond the clipping point is of order 1/(eta L).

[Numerical] Status of (D). At half filling, eta L sup over roots |sigma_L - r_0| is between 0.230 and 0.317 for epsilon=10^-2,...,10^-6 and L=16,...,8192; the maximum sits on isolated tail roots, where sigma_L has a spike of height near 1/(pi eta L). C136's bound for the same quantity is 10.8 to 1061. At a finite Fermi boundary the maximum over [-Q_L,Q_L] is at most 0.28 in seven cases with epsilon=10^-3 to 10^-6 and eta L from 2.9 to 29. dens_check2.py extends this to epsilon in {0.9, 0.5, 0.1, 10^-2, 10^-4, 10^-6}, L in {16, 64, 256, 1024} and n from 2 to L/2-2 (120 cases): the maximum of eta L sup|sigma_L - r_Q| over [-Q_L,Q_L] is 0.317, 0.302, 0.185, 0.054, 0.009, 0.002 at epsilon = 10^-6, 10^-4, 10^-2, 0.1, 0.5, 0.9. For (D1), eta^2 L sup|sigma_L' - r_Q'| on [0, Q_L - min(eta,Q_L)/4] is at most 0.22 for epsilon <= 0.1. A first pass at epsilon in {0.5, 0.9} gave values up to 3.6, growing linearly in L; this was a mesh artifact (a one-sided difference for r_Q' at x=0, where r_Q'(0)=0 by evenness: the error halves when the mesh is doubled and does not depend on L). With r_Q'(0)=0 imposed the maximum is 0.019 (dens_check2_0.9_0.5.json). The small-epsilon values are upper estimates that may contain the same artifact.

[Numerical] Status of the continuum hypotheses. Var(v_Q)/eta <= 8.0 (equality at the free point); e_bulk''/eta^2 <= 4/pi (equality at the free point, about 1.0 near Delta=-1); r_Q/r_0 >= 0.762 at distance >= 1 from the Fermi point, for Q from 0.02 to 8 and epsilon from 1 to 10^-6. Two cases at epsilon=10^-6, Q>=4 are underresolved (they show r_Q/r_0 up to 1.05, violating the proved r_Q<=r_0).

[Numerical, negative] The C136 currency does not repair. The cumulative root error L sup|G| is 0.67 to 0.76 at half filling at every epsilon, but at a finite Fermi boundary it is concentrated at the outermost roots and grows as eta->0: the mismatch L(nu - m(Q_L)) is 0.32, 0.77, 1.6 at fixed eta L between 18 and 29 for epsilon=10^-2, 10^-3, 10^-4, and 4.6 at epsilon=10^-6, L=2048, n=512. A uniform bound L||G|| <= C therefore appears false near a Fermi boundary, while (D) stays bounded. This is why Theorem B is phrased with the pointwise density error.

## Numerics for (X) and the doubling defect

[Numerical] Doubling defect, scan_defect.py, epsilon=10^-k with k=2,...,6, L=4,...,4096 paired with 2L, n in {1,2,3,4,8, sqrt L, L/2-sqrt L, L/2-2, L/2-1, L/2} and dyadic densities 1/128 to 1/2:

| epsilon | eta | pi v/2 | max L Lambda/M_0 | max over pi v/2 | max abs L(F_L - L e_bulk) |
| --- | --- | --- | --- | --- | --- |
| 1e-2 | 0.1415 | 0.1160 | 0.1160 | 1.0000 | 0.0386 |
| 1e-3 | 0.0447 | 0.0356 | 0.0356 | 0.9999 | 0.0118 |
| 1e-4 | 0.0141 | 0.0112 | 0.0112 | 0.9993 | 0.0036 |
| 1e-5 | 0.00447 | 0.00352 | 0.00349 | 0.9926 | 0.0010 |
| 1e-6 | 0.00141 | 0.00111 | 0.00101 | 0.9132 | 0.00018 |

Here L Lambda = L[F_2L(2n) - 2F_L(n)], M_0=J/2, and v=(pi J/2)sin(gamma)/gamma. In every case 0 < L Lambda/M_0 <= pi v/2, the half-filled Casimir limit; the ratio increases with L at fixed density, passes one half near eta L between 3 and 10, and is at most 0.004 at L<=7 for epsilon=10^-6. The last column uses the exact ladder identity L(F_L - L e_bulk) = -(1/2) sum_j Lambda_j/4^j with Lambda_j = 2^j L[F_(2^(j+1)L)-2F_(2^j L)], truncated after at least three terms (tail below 1e-6), and is at most pi v/6 in magnitude. No C129 band or dilute violation occurs. So (X) and its one-sided version hold numerically with K_* of order eta J, three orders of magnitude below the conditional 125.

[Numerical, consistency checks] Delta=-1 limit: at eta=1e-5 the Bethe energy equals -epsilon n(L-n)/(L-1) through the proved sum rule to 1e-8. C129 band: L Lambda/M_0 <= epsilon L^2 in all cases. C136 constant at fixed delta: the finite-boundary identity C136 (4.3), evaluated with the new continuum solver, agrees with the ladder to 4e-7 at epsilon=10^-3, L=64. Casimir: at half filling L(F_L - L e_bulk)/(-pi v/6) at L=8192 is 1.0000, 1.0000, 0.9999, 0.9985, 0.9847 for epsilon=10^-2,...,10^-6 (diag_half.py, closed-form e_bulk(1/2)), and L Lambda/M_0 tends to pi v/2. ED: the solver matches exact diagonalization to 5.8e-15 at L<=12.

[Numerical, scope] Continuum diagnostics near half filling with Q_L of order log L at epsilon <= 10^-5 are mesh-limited and are excluded (diag_identity.py rows with n=L/2-1 at k>=4 have mesh changes larger than the quantity). Typhon was not used: all Bethe solves through L=8192 took minutes locally, and no computation was on the critical path.

## Files

bethe_ferro.py (solver), scan_defect.py and scan_K.json, energies_K.json (defect scan), continuum.py (restricted dressing), diag_half.py and diag_half.json, diag_identity.py and diag_K.json, conv_check.py and conv_check.json, dens_check.py and dens_check.json, dens_check2.py and dens_check2.json, cont_diag.py, cont_diag_wide.py and cont_diag.json, cont_diag_wide.json, check_c129.py (ED, copied from C136).
