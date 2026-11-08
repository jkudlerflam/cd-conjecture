This run was done by a Claude Opus agent (2026-11-07), not by Codex.

# C149 result: the ferromagnetic strip by a discrete comparison principle

Proofs are in [PROOF.md](PROOF.md), the log in [PROGRESS.md](PROGRESS.md). No independent referee pass. J = 1, eta = pi - gamma, Delta = -cos(eta), Q_nu the continuum Fermi point (m(Q_nu) = nu = n/L), Z_0 the continuum counting function of r_(Q_nu), s_nu = Z_0'.

## Bottom line

[Assumed, not closed] The strip is not closed unconditionally. On this route it is reduced to one inequality about explicit trial configurations, with no Bethe root in it:

    (R_kappa)  the quantiles y^+ of Z_0 at the shifted levels I_j/L + kappa0/L satisfy F(y^+) >= 0
               (equivalently, by symmetry, the quantiles y^- at I_j/L - kappa0/L satisfy F(y^-) <= 0).

[Proved, conditional on (R_kappa) with kappa0 = 0.4] For every -1 < Delta <= 0 (eta_0 = pi/2, the whole attractive side), every L >= 2 and 0 <= n <= L/2,

    |F_L(n) - L e_bulk(Delta, n/L)| <= 136.4/L,     L(E_(2L) - 2E_L) <= 685 M,

for every bond of the U(1) class with -1 < J_z/J <= 0, every field and DM angle. With C142 on the repulsive side the whole U(1) class has C = 4078. Neither (G), (W1') nor (W3') is used. For C150's route the same hypothesis gives sigma_L >= 0.714 r_(Q_L) on its shell S' whenever Q_L <= Q_nu (the case met in every sampled configuration), so L|D| <= 134 there.

[Numerical] (R_kappa) at kappa0 = 0.4 holds in all 173 sampled cases (1 + Delta from 0.5 to 1e-6, L from 16 to 8192, all fillings), with the smallest margin 0.088 eta/L at half filling and 1 + Delta = 1e-6, close to its predicted asymptote (2/pi)(kappa0 - b/4) eta/L = 0.096 eta/L as eta -> 0.

## The key lemmas

[Proved] Comparison principle (Theorem 1.1). F_j(x) = L p(x_j) - sum_k theta(x_j - x_k) - 2 pi I_j satisfies: F(y) >= F(x) componentwise implies y >= x, for all real configurations, all eta in (0, pi/2] and all real quantum numbers, with no ordering hypothesis. Proof by a maximum principle at the index maximizing x_j - y_j, using only that theta is nonincreasing and p strictly increasing (Rheinboldt's M-functions). Uniqueness of the real Bethe solution is a corollary.

[Proved] Residual identity (Lemma 2.1). For every configuration y, F_j(y)/(2pi L) = [Z_0(y_j) - I_j/L] - T*(N_0 - N_y)(y_j), with N the cumulative counting functions. For the shifted quantiles this is +kappa minus a T-average of kappa plus a sawtooth, so (R_kappa) is a first-order quadrature inequality for an explicit regular configuration against the gain (1 - q) kappa = (2 eta/pi) kappa of a twist.

[Proved, conditional on (R_kappa)] The conversion (Part 3). The sandwich y^- <= x <= y^+ gives the cumulative bound L sup |N_x - N_0| <= kappa0 + 1/2 relative to the continuum at Q_nu (Lemma 3.1). The kernel factorizes, T = T_A * T_B, into two positive unimodal kernels of the same family (Lemma 3.2, a Fourier identity). One integration by parts against T', or against T_A' after T_B, then gives log-free (D) and (D1) at every point of the real line (Theorem 3.3):

    eta L |sigma_L - s_nu| <= (2kappa0 + 1) eta cot(eta)/pi <= 0.573,
    eta^2 L |sigma_L' - s_nu'| <= (2kappa0 + 1) 2 eta^2 T_A(0) T_B(0) <= 1.459      (kappa0 = 0.4),

against C142 Theorem II.6's (4/pi)(log(1/eta) + pi/2) and its derivative analogue. The logarithm of C142 came from charging each cell at first order through the full dressing kernel; here the roots are controlled at the level of counting, which costs no logarithm.

[Proved] Energy identity at Q_nu (Lemma 4.1). F_L(n) - L e_bulk(nu) = sum_j eps^(x_j) - L int_(-Q_nu)^(Q_nu) eps_Q sigma_L, with eps^ the natural extension of eps_(Q_nu). There is no density mismatch term (the term that needed (W3)), only a boundary layer between Q_nu and Q_L, controlled by Corollary 3.4. C138 Theorem B then runs with the constants above (Theorem 4.3).

## Strip constant (kappa0 = 0.4, constants_C149.py)

| term | sup over eta | eta -> 0 |
| --- | --- | --- |
| core (2a/9)(1 + C_D1 A_core/(eta C_D)) | 5.53 | 4.83 |
| velocity variation V_W(eta)/6 | 13.43 | 9.87 |
| first-order tail 16 (a/eta) C_D A_tail | 114.88 | 114.88 |
| straddling cells 16a/3 | 16.76 | 0 |
| cells beyond Q_nu, 4aq(kappa0 + 1/2) | 4.04 | 0 |
| boundary layer 4aq^2(kappa0 + 1/2)^2 | 2.27 | 0 |
| K (pointwise sup of the sum, at eta = 0.289) | 136.39 | 129.58 |

C_strip = 5K + pi^2/4 = 684.4. For kappa0 = 0.3, 0.35, 0.45: 619.9, 652.1, 716.9.

## Numerics for (R_kappa)

rk_test.py evaluates F(y^+) exactly for the explicit trials. The continuum is solved by a new solver (fftsolver.py: uniform grid h = eta/10, Gregory end weights of order 7, FFT Toeplitz products, CG; m(Q) converged to 1e-12), because the graded Nystrom solver of C138 is mesh-limited near Delta = -1 (errors of 5 to 10 percent in r_Q(0) at 1 + Delta = 1e-6, amplified by B). Fillings n in {L/32, L/8, L/4, 3L/8, L/2 - L/16, L/2 - 8, L/2 - 2, L/2 - 1, L/2}, 173 cases in all. Entries: min over j and n of F_j(y^+)/(2 pi L) in units of eta/L.

| 1 + Delta | L | kappa0 = 0.25 | 0.3 | 0.4 | 0.45 | worst case | max true displacement L|Z_0(x_j) - I_j/L| |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.5 | 16 to 1024 | 0.154 | 0.187 | 0.246 | 0.275 | n = L/2 | 0.027 |
| 0.1 | 16 to 1024 | 0.102 | 0.141 | 0.220 | 0.260 | n = L/2 | 0.095 |
| 1e-2 | 64 to 4096 | 0.030 | 0.066 | 0.141 | 0.179 | n = L/2 | 0.176 |
| 1e-4 | 256 to 4096 | -0.008 | 0.025 | 0.093 | 0.129 | n = L/2 | 0.249 |
| 1e-6 | 16 to 8192 | -0.035 | 0.020 | 0.088 | 0.124 | n = L/2 (n = 7L/16 for 0.25) | 0.263 |

Consequences checked on Bethe roots (check_density.py, 30 cases, 1 + Delta = 1e-2, 1e-4, 1e-6, L to 4096, n from L/8 to L/2 - 1): L sup|N_x - N_0| <= 0.758 (bound kappa0 + 1/2 = 0.9), eta L sup|sigma_L - s_nu| <= 0.317 (bound 0.573), eta^2 L sup|sigma_L' - s_nu'| <= 0.217 (bound 1.459), the factorization T = T_A * T_B to 1.5e-11 at 1 + Delta = 1e-2 (2e-5 at 1e-6, limited by the check's quadrature grid), and Q_L < Q_nu in every case. The true displacement L|Z_0(x_j) - I_j/L| never exceeds 0.263, inside the sandwich width 0.4 that (R_kappa) predicts.

Consistency checks (PROOF 5.3). Free point: (R_kappa) is exact, F_j(y^+) = 2 pi L kappa. eta -> 0 at fixed L: (R_kappa) divided by eta tends to the same statement for the convex log-gas of C138 Section 2.3. C142 Theorem II.6: the logarithmic C_D(eta) = 9.97 at 1 + Delta = 1e-6 is replaced by 0.573. C138's observed growth of the cumulative error near a Fermi boundary is an artifact of normalizing at Q_L; at Q_nu the cumulative error is bounded.

## What remains

(R_kappa) itself. It needs (i) a sharp quadrature lemma for a T-average of the sawtooth at the nodes of an explicit quantile configuration, with the sparse-limit constant b/(2 pi) (a crude first-order bound loses a factor 2, which is the whole margin between kappa0 > b/4 and kappa0 < 1/2 at half filling), and (ii) control of r_(Q_nu) near the Fermi edge on the kernel scale, where only the reflection bound (distance >= 0.7) and the log-Lipschitz bound |r'| <= kappa_eta r are available. Numerical margin at kappa0 = 0.4: at least 0.088 eta/L in units of the residual, against the twist gain (2/pi) 0.4 = 0.255 eta/L, so the quadrature term uses at most 65 percent of the gain.

Files: rk_test.py (rk_*.json, log_rk_*.txt), check_density.py (check_density_*.json), constants_C149.py (constants_C149.json), fftsolver.py, trial.py, explore_half.py, explore_disp.py, scan_phi.py, barrier_test.py; bethe_ferro.py, continuum.py, check_attractive.py copied from work_C146. All runs local; Typhon was not used.
