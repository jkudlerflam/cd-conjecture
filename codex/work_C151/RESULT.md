This run was done by a Claude Opus agent (2026-11-07), not by Codex.

# C151 result: the ferromagnetic strip through the Destri-de Vega representation

Proofs are in [PROOF.md](PROOF.md), the log in [PROGRESS.md](PROGRESS.md). No independent referee pass of the new material. J = 1, eta = pi - gamma, Delta = -cos eta, lambda = eta L r_Q.

## Bottom line

[Assumed] The strip is not closed unconditionally. The one statement still missing is a positive lower bound on the density on the shell, sigma_L >= c r_Q on S = {M_0 <= lambda < 4C_D(eta)} for some fixed c > 0. The DdV representation does not supply it (Part 5 of PROOF), for the same underlying reason C145 and C146 found: it needs regular spacing of the shell roots.

[Proved, conditional on that statement only] For every eta in (0, pi/4], L >= 2, 0 <= n <= L/2,

    |F_L(n) - L e_bulk| <= K(c, M_0)/L,     L(E_(2L) - 2E_L) <= (5K + pi^2/4) M,

with K = 76.72 and C = 386.1 at c = 3/4, M_0 = 2. That hypothesis is C146's (G_c) word for word, from which C146 got only a log log(1/eta) constant; with the DdV representation it gives a uniform one. Any c > 0 works: K = 106.7, 202.9, 491.7 at c = 1/2, 1/4, 1/10 (C = 536, 1017, 2461). No derivative bound ((D1)), no velocity variation ((W1)) and no variation of r/sigma on the shell are used.

[Proved, unconditional] |F_L(n) - L e_bulk| <= K_0(eta)/L on (0, pi/4] with K_0 = 202, 278, 426, 570, 994 at eta = 1e-1, 1e-2, 1e-4, 1e-6, 1e-12, about 31.4 log(1/eta) + 140. Still a logarithm, but about one eighth of the slope of C142 (247 log(1/eta)) and with no unproved continuum input.

## The key estimates

1. [Proved] Positivity in the strip. With z = x + is, Re T(z) = (sin 2eta/pi)(cosh 2x cos 2s - cos 2eta)/|cosh 2z - cos 2eta|^2 and Re a1(z) is positive in the same way, so Re sigma_L > 0 on |Im z| < min(eta, pi/4). Hence Im theta(x + is) > 0 for s > 0, the only zeros of 1 + (-1)^n e^(2 pi i theta) near [-Q, Q] are the Bethe roots, and at the counting boundary theta(+-Q) = +-n/2 the function equals 2.

2. [Proved] Exact identity at every filling. For f analytic near the rectangle and real on R,

       sum_j f(x_j) - L int_(-Q)^Q f sigma_L = -(1/pi) Im int_Gamma f'(z) Log(1 + (-1)^n e^(2 pi i theta(z))) dz,

   for any path Gamma from Q to -Q in the closed upper half of the rectangle. With f = eps_Q this is F_L - L[E + h(nu - m)] (C136 (4.3)). On real stretches the integrand is the sawtooth (half the C138/C142 oscillation charge); on lifted stretches it is bounded by |f'| G(Psi), Psi = 2 pi Im theta, G(t) = -log(1 - e^(-t)) + arcsin(e^(-t)).

3. [Proved] Harnack. Re sigma_L is a positive harmonic function on the disc of radius eta about any real x, so

       2 pi Im theta(x + i eta/2) >= alpha_H eta L sigma_L(x),   alpha_H = 2 pi (2 log 1.5 - 1/2) = 1.954.

   The decay of the nonlinear term above x is controlled by the real density at x alone. This is what removes every derivative and variation bound: the lifted integral over the dense set is at most (4 kappa_eps ell_B/L) sum_k 2^(k+1) M_0 G(alpha_H c 2^k M_0), a geometric series.

4. [Proved] Strip bound on the dressed energy, |eps_ext'(x + is)| <= 2a kappa_eps r_Q(x) for s <= eta/2, kappa_eps <= 2 (4/3 as eta -> 0), from |T(x + is)| <= [sin^2 eta/(sin^2 eta - sin^2 s)] T(x) and |eps_Q'| <= 2a r_Q.

Contour: real axis on the sparse tail {lambda < M_0} (cost 2 pi M_0 A_tail, A_tail = 3.461), lifted at height eta/2 over K_0 = {lambda >= M_0}, vertical connectors at cell boundaries (or the Fermi-edge pieces at +-Q, which carry the CFT term), plus the mismatch pi^2 C_e/8. Terms at c = 3/4, M_0 = 2: tail 43.50, straddling cells 5.92, connectors 5.28, lifted core 11.62, mismatch 10.40.

## Why the shell bound does not follow (PROOF Part 5)

[Proved, computations] (i) Bootstrapping sigma >= c r through the dressed DdV equation leaves a Hilbert-transform sum over scales eta < u < 1 of size G(alpha_H c lambda(x)) log(1/eta), so it closes only for lambda >= (1/(alpha_H c)) log log(1/eta). (ii) The sparse tail next to the shell must be charged first order, which reproduces C145's geometric obstruction at the inner edge of any hypothesis region. (iii) Using the oscillation of the nonlinear term to remove the log needs |sigma'(y + i eps)| of order sigma on the unit scale, that is regular spacing of the shell roots; Harnack only gives (2/eta) sigma. [Heuristic] For a locally periodic root lattice the mean of the nonlinear term over a period is exactly zero, which is why the true answer has no logarithm; excluding locally irregular configurations needs the Bethe equations beyond positivity (for instance Yang-Yang convexity).

Relation to C150's target (H_I), int over {2 <= lambda < 2C_D} of dy/(eta L sigma) <= C_I. [Proved] sigma >= c r on the shell implies (H_I) with C_I <= 6.07/c, and implies the hypothesis of Theorem 2.1 here. Neither (H_I) nor the logarithmic variant of PROOF 2.3 (eta L sigma >= 0.52(1 + delta) log lambda, enough for the lifted integral) implies the other. Nothing here proves any of them.

## Numerical checks

| check | result |
| --- | --- |
| theta(x_j) = I_j from the continued counting function | 2e-14 |
| DdV identity, half filling, eta in {pi/2, 0.5, 0.1, 0.01}, L in {16, 64, 256}, heights eta/8, eta/4, eta/2 | equals L(F - L e) to 1e-7 or better, height independent |
| free point Delta = 0 | L(F - L e) = -0.52596, -0.52375, -0.52361 at L = 16, 64, 256; -pi/6 = -0.52360 |
| CFT, eta = 0.5, half filling | -0.149262 at L = 256 against -pi v/6 = -0.149275 |
| DdV identity, finite filling (8 cases, eta in {0.3, 0.1}) | equal to nine digits at heights eta/4, eta/2 |
| Fermi-edge pieces against -pi v_F/(6L), lambda(Q) = 4.3 | -0.0719834 (edge), -0.0719833 (total), -0.0719943 (CFT) |
| Harnack bound, 46 Bethe cases, 1 + Delta to 1e-6, L to 16384 | Psi/(alpha_H eta L sigma) >= 1.532; on the shell Psi/(pi lambda) >= 0.991 |
| shell comparability, 34 cases with S nonempty | min sigma/r = 0.9916 (two near-half-filled cases at 0.996 and 0.992 are mesh-limited; all others >= 1.000) |
| max eta L |sigma - r| | 0.317 |

Shell table (check_shell_1e-*.json):

| 1 + Delta | cases | L max | min Harnack ratio | min sigma/r on S | min Psi/(pi lambda) on S | min eta L sigma/log lambda on S |
| --- | --- | --- | --- | --- | --- | --- |
| 1e-2 | 15 | 2048 | 1.535 | 1.0000 | 0.991 | 2.72 |
| 1e-4 | 10 | 4096 | 1.532 | 0.9964 | 0.996 | 2.71 |
| 1e-5 | 5 | 8192 | 1.532 | 0.9916 | 0.991 | 2.72 |
| 1e-6 | 5 | 16384 | 1.532 | 1.0000 | 0.9997 | 2.72 |

DdV nonlinear integral equation at half filling (nlie_main.json), solved without Bethe input (Picard at eta >= 0.5, Newton-Krylov with continuation in eta from 0.5 below), compared with Bethe roots:

| eta | L | max abs(theta_NLIE - theta_Bethe) on Im z = eta/4 | max root difference | L(F - L e), NLIE | L(F - L e), Bethe |
| --- | --- | --- | --- | --- | --- | --- |
| pi/2 | 256 | 6e-12 | 2.6e-7 | -0.5236079663 | -0.5236079710 |
| 0.5 | 256 | 5e-14 | 4.2e-7 | -0.1492621124 | -0.1492621129 |
| 0.1 | 256 | 1.5e-13 | 5.1e-7 | -0.0269213968 | -0.0269213969 |
| 0.03 | 256 | 4.9e-10 | 4.9e-7 | -0.0076415793 | -0.0076415793 |
| 0.03 | 64 | 2.7e-5 | 3.9e-6 | -0.0049105805 | -0.0049105733 |
| 0.01 | 16 | 1.8e-4 | 5.8e-6 | -0.0002066422 | -0.0002066415 |
| 0.01 | 64 | 6.5e-5 | 2.0e-6 | -0.0007107523 | -0.0007107539 |

(Root differences are limited by linear interpolation of the counting function on the grid; the small-L differences at eta <= 0.03 come from grid truncation, visible as a residual 2e-7 to 5e-7 of the Bethe counting function in the equation. Plain Picard iteration diverges for eta <= 0.1, reflecting the norm pi/(2 eta) of the dressing. At eta = 0.01 the continuation needed steps of ratio 1.1 (nlie_eta001.json); L = 256 at eta = 0.01 was not run.)

## Files

ddv.py (counting function, Log(1 + w), contour quadrature), check_half.py and check_half.txt (identity at half filling, free point), check_finite.py and check_finite.json (finite filling, Fermi-edge pieces, C136 identity), check_shell.py and check_shell_1e-*.json with logs (Harnack bound, shell margins), nlie.py, nlie_main.json and nlie_eta001.json with logs log_nlie*.txt (DdV NLIE solve), constants_C151.py and constants_C151.json (every constant above). Solvers bethe_ferro.py and continuum.py copied from work_C146. All runs local, each under about ten minutes; Typhon was not needed.
