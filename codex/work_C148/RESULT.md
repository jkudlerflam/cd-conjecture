This run was done by a Claude Opus agent (2026-11-07), not by Codex.

# C148 result: the easy-axis side of the U(1) class, J_z/J > 1

Proofs are in [PROOF.md](PROOF.md), the log in [PROGRESS.md](PROGRESS.md). No referee pass yet. J = 1, Delta = cosh(lambda), M_0 = (1 + Delta)/4.

## Bottom line

[Proved, one step computer-checked] For every bond of the U(1) class with J_z/J > 1, every field, DM angle and L >= 2,

    L(E_(2L) - 2E_L) <= 10781 M.

This comes from the sector theorem (PROOF Theorem A), uniform in lambda > 0, L >= 2 and 0 <= n <= L/2,

    |F_L(n) - L e_bulk(Delta, n/L)| <= K/L,   K <= 218.4 ab <= 4311 M_0,

since ab/M_0 = 4pi^2 tanh(lambda/2)/lambda <= 2pi^2. Matched doubling gives (5/2)K <= 10778 M_0, and the C129 twist cost adds pi^2/4. The exact lambda-dependent constant evaluated on a grid gives 4153 M_0 and C = 10386 (constants.json). The sup is approached as Delta -> 1+. For Delta >= cosh(1.36) = 2.08, K/M_0 <= 462, and K/M_0 -> 0 like 1/lambda as Delta -> infinity.

The computer-checked step is one inequality between monotone elementary functions of lambda (PROOF Theorem 4.4), bracketed on 36 cells of [1, 1.36] with largest value 0.939 < 1 (pos_rigorous.py). Interval arithmetic would make it fully formal.

Gluing at Delta = 1. Both C142 (4078 on 0 <= J_z/J <= 1) and Theorem A are uniform up to Delta = 1 from their sides, with the same C129 reductions and M_0 = J/2 there, and the identification is a single analytic Bethe branch through Delta = 1 in XXX rapidities. So C = 10781 holds on all of J_z/J >= 0, and with C136 the constant on J_z/J >= -1 + delta is max(10781, 6e15 delta^(-25/3)).

## How it works

The Bethe roots live on the circle T = R/(pi Z), with a1 and K = a2 Poisson kernels. This is the repulsive structure of C139/C140/C142 in the limit gamma -> 0 (||K||_1 = 1, hole kernel ||J||_1 = 1/2), on a compact rapidity domain. The C139/C140/C142 tools carry over as follows.

- [Proved] The hole kernel factorizes as J = a1 (*) r0, with r0 = (K_e/pi^2) dn(2K_e x/pi). It is positive and circle-unimodal, by a circle version of the convolution lemma (PROOF Lemma 2.1). Also sum_(m>=1) J^((*)m) = K exactly.
- [Proved] Duality, the density bounds r0 <= r_Q <= 2r0, 1/2 <= Z <= 1, Var(Z) <= 1/2, the restricted inverse <= 3, eps_Q' >= 0, |v| <= ab, Var(v) <= ab(4 + 12 kappa_w) and e_bulk'' <= 2ab all hold uniformly in the Fermi point Q in (0, pi/2], half filling (Q = pi/2) included (PROOF Section 3).
- [Proved] The counting-density lower bound. The C140 pairing argument works on the circle with two changes: the antipode case, and the two-sided out arc, so both outer roots are left over. The outermost-root location is improved by elliptic functions to L r0(x_n) >= 4 sin(pi/8) k K_e/pi^2. Together they give sigma_L >= (1 - Phi) r0(x_n) on [x_1, x_n] with Phi <= 0.804 for lambda <= 1 (analytic) and Phi <= 0.939 for 1 <= lambda <= 1.36 (cells). For lambda >= 1.36, sigma_L > 0 on all of T for every configuration, directly from min a1 > (1/2) max K. The negative mass beyond the outermost roots is at most 1/L per side.
- [Proved] Identification at both parities. C142's branch is continued through Delta = 1 in XXX rapidities u, where the equations are analytic in s = gamma^2 = -lambda^2. It then continues along lambda in (0, infinity) by the open-closed argument. E = F_L(n) follows by the identity theorem, since F_L(n) is analytic (Perron-Frobenius at flux L pi).
- [Proved] Sawtooth s = 3/2 (s = 1/2 for lambda >= 1.36), ||G|| <= 13.5 s/L, density errors D/L and D1/L, quadrature on the core, first-order ends and curvature, all as in C142 Lemmas I.5 to I.8 with circle tail bounds (PROOF Sections 5 to 6). No lower threshold on L is needed.

The spec's four regions needed no separate theorem. (a) Delta near 1: the constants have finite lambda -> 0 limits. (b) Near h_c1: dilute holes are a short out arc, uniform in Q. C129's dilute bound is not used. (c) Half filling, handled directly; parity is discussed below. (d) Large Delta: K/M_0 -> 0.

## Numerics

- [Numerical] Bethe = ED in all 324 sectors, L = 2..12, Delta in {1.001, 1.01, 1.1, 1.5, 2, 3, 5, 10, 20}, max error 9.2e-14 (ed_checks.json).
- [Numerical] End-to-end ED L(E_2L - 2E_L)/M with 43 fields in [0, 1.05 h_c2] (including just above h_c1), DM in {0, 0.3, 1}, L = 2..6, 5805 samples: max 3.9993 (Delta = 1.001, L = 2, B = 0), continuous with Heisenberg's 4.
- [Numerical] Bethe doubling defect L(F_2L(2n) - 2F_L(n))/M_0, 621 cases to 2L = 1024: max 2.79 (L = 4), 2.47 at L = 128 (Delta = 1.001), 0.17 at Delta = 20. No solver failures.
- [Numerical] Intermediate inequalities (261 cases, L to 64). min sigma_L/r0 on T >= 0.969 (proved >= 0.06 on [x_1, x_n]). L r0(x_n) >= 0.64. Negative mass is 0. L sup|Cum R| = 0.500 (proved 3/2). L||G|| <= 0.57 (proved 20.25). L|F - L e_bulk|/M_0 <= 1.82 (proved 4311).
- [Numerical] h_c1 = 2a r0(pi/2) = (2/pi) sinh(lambda) K(k) k' matches F_L(L/2 - 1) - F_L(L/2) at L = 1024 to 4e-6 (Delta = 2, 5).
- Half-filling parity (task (c)) [Numerical]. Even L: the defect is exponentially small beyond the correlation length (7e-5 at Delta = 2, L = 130). Odd L: the spinon sector gives a positive defect that decays with L (L times defect over M_0 is 1.82, 1.38, 0.51 at Delta = 2 and L = 9, 33, 129).

## Open

Nothing on J_z/J > 1 at the level of the target. Remaining items: interval-arithmetic certification of the 36-cell check, and a referee pass. The constant is crude, mainly because s = 3/2 (the out arc is two-sided on the circle) and because the ends are first order.

## Files

circle.py (circle kernels, restricted continuum, finite-L diagnostics), bethe_scan.py and check_c129.py (copied from C142), scan_defect.py (scan_small.jsonl), scan_large.py (scan_large.jsonl), ed_checks.py (ed_checks.json), pos_window.py and pos_check.py (pos_check.json), pos_rigorous.py (pos_rigorous.json), constants.py (constants.json), quick_test.py. Typhon was not used, since each run took minutes locally.
