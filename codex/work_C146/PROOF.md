# C146 proof notes: the edge bound, the velocity variation, and the shell

2026-11-07. Written by a Claude Opus agent, not by Codex. J=1 throughout. Tags: [Proved] is a mathematical deduction, [Proved, cited] uses an earlier run (named, with its referee status), [Numerical] is floating point without enclosures, [Heuristic] has an unproved step, [Assumed] is a hypothesis.

Notation of C145 PROOF Section 1. eta = pi - gamma in (0, pi/2], Delta = -cos eta, q = 1 - 2eta/pi, B = 1/(1-q) = pi/(2eta), a = pi sin eta, b = pi/gamma in [1, 2), T(x) = sin(eta)cos(eta)/(pi(sinh^2 x + sin^2 eta)) >= 0 even and decreasing in |x|, with Fourier symbol T^(w) = sinh((pi-2eta)w/2)/sinh(pi w/2) in [0, q]. a1(x) = sin(eta)/(pi(cosh 2x + cos eta)), e = -a a1. On [-Q, Q], T_Q f(x) = int_(-Q)^Q T(x-y) f(y) dy, R = (I - T_Q)^(-1) = sum_k T_Q^k (positive, norm at most B on L^infinity), r = R a1, Z = R 1, eps_Q = R(e - h) with h = -a r(Q)/Z(Q), v = eps_Q'/r. r_0 = 1/(2 gamma cosh bx) is the full-line density.

Facts used from earlier runs, all [Proved, cited]: r and Z are even and nonincreasing on [0, Q], a1 <= r <= r_0, 1 <= Z <= B (C136 Section 2, refereed); eps_Q' > 0 on (0, Q] and (I - T_Q) eps_Q' = e' = -a a1' (C136 (5.1), refereed); 0 <= v <= 2a tanh(2Q) on [0, Q] (C145 Lemma B.2, not refereed, re-derived below where used); e_bulk''(nu) = v(Q)/(2 Z(Q)^2) at nu = m(Q) (C136 (5.3), refereed); 0 <= r_0 - r_Q <= r_0(2Q - x) + r_0(2Q + x) on [-Q, Q] (C145 Theorem B.1, re-derived in C145 and checked here line by line); Theorem A.2 of C145 (the conditional strip theorem; its ingredients from C138 Theorem B were refereed in C145 Part C).

Summary.

- [Proved] Part 1, (W3). eta Z_Q(Q)^2 >= pi Q/(4(1 + 16Q/(3 pi^2))) for every eta and Q. Hence (W3') holds with c_Z = 0.373 (and with the clean value 0.309), and e_bulk'' <= C_e eta^2 with C_e = 8.43.
- [Proved] Part 2, (W1). The monotonicity (W1') is not proved, but the bound it was used for is: Var_[-Q,Q](v_Q) <= V_W(eta) <= 92.9 for every eta and Q (59.2 as eta -> 0). Theorem A.2 only needs Var(v_Q) bounded, not of order eta.
- [Assumed] Part 3, (G). Not proved. The derivative half of (G) is reduced to the comparability half, at the price of a factor of order log log(1/eta) (Proposition 3.1). What is left is one one-sided pointwise statement, sigma_L >= (3/4) r_Q on the shell.
- Part 4. Assembly: C_strip = 773 conditional on (G) alone; C_strip(eta) <= 868 for eta >= 1e-12 (growing like log log(1/eta)) conditional on the comparability half of (G) alone.

# Part 1. The edge bound (W3)

## 1.1 Two exact facts about Z_Q(Q)

[Proved] Lemma 1.1 (domain derivative). Let W(Q) = int_(-Q)^Q Z_Q(x) dx. Then W is C^1 on (0, infinity), W(0+) = 0, and

    W'(Q) = 2 Z_Q(Q)^2.

Proof. Differentiate Z_Q(x) = 1 + int_(-Q)^Q T(x-y) Z_Q(y) dy in Q (T is smooth, so Z_Q is smooth in (x, Q)): with Z_Q(-Q) = Z_Q(Q) by evenness,

    (I - T_Q) d_Q Z_Q = Z_Q(Q) [T(. - Q) + T(. + Q)],   so  d_Q Z_Q = Z_Q(Q) R[T(. - Q) + T(. + Q)].

Then W'(Q) = 2 Z_Q(Q) + int d_Q Z_Q = 2Z_Q(Q) + Z_Q(Q) <1, R[T(.-Q) + T(.+Q)]>. R is symmetric, so <1, R f> = <Z_Q, f>, and <Z_Q, T(. - Q)> = (T_Q Z_Q)(Q) = Z_Q(Q) - 1, the same at -Q. Hence W' = 2Z(Q) + 2Z(Q)(Z(Q) - 1) = 2Z(Q)^2. W(Q) <= 2QB gives W(0+) = 0. End.

This is the Z-analogue of m'(Q) = 2 r_Q(Q) Z_Q(Q) (C136 (2.3)).

[Proved] Lemma 1.2 (edge monotonicity). Q -> Z_Q(Q) is nondecreasing on (0, infinity).

Proof. Let Q' > Q, delta = Q' - Q and u(x) = Z_(Q')(x + delta) on I' = [-Q - 2delta, Q]. Shifting variables in the equation for Z_(Q'), u = 1 + T_(I') u on I'. For x in I = [-Q, Q], contained in I', the part of T_(I') u coming from I' minus I is nonnegative (T >= 0, u >= 1), so (I - T_Q)(u - Z_Q) >= 0 on I. Positivity of R gives u >= Z_Q on I, and at x = Q this reads Z_(Q')(Q') >= Z_Q(Q). End.

[Proved] Corollary 1.3. W is convex on (0, infinity), and

    Z_Q(Q)^2 >= W(Q)/(2Q).

Proof. W' = 2Z_t(t)^2 is nondecreasing in t by Lemma 1.2. So W(Q) = int_0^Q 2Z_t(t)^2 dt <= 2Q Z_Q(Q)^2. End.

## 1.2 A variational lower bound on W

[Proved] Lemma 1.4. For every eta and w, 0 <= 1 - T^(w) <= eta(|w| + 2/pi).

Proof. For w > 0, sinh(pi w/2) - sinh((pi - 2eta)w/2) = int over [(pi-2eta)w/2, pi w/2] of cosh t dt <= eta w cosh(pi w/2), so 1 - T^(w) <= eta w coth(pi w/2). With coth y <= 1 + 1/y (equivalent to e^(2y) - 1 >= 2y), w coth(pi w/2) <= w + 2/pi. The symbol is even. Positivity is T^ <= q < 1. End.

[Proved] Lemma 1.5. For every eta in (0, pi/2] and Q > 0,

    W(Q) >= pi Q^2 / (2 eta (1 + kQ)),     k = 16/(3 pi^2) = 0.5404.

Proof. I - T_Q is a positive definite symmetric operator on L^2[-Q, Q] (norm of T_Q at most q < 1), so W(Q) = <1, (I - T_Q)^(-1) 1> = sup_f [2<f, 1> - <f, (I - T_Q) f>]. For f supported in [-Q, Q] and extended by zero, <f, T_Q f> = <f, T*f>, so by Plancherel and Lemma 1.4

    <f, (I - T_Q) f> = (1/2pi) int |f^(w)|^2 (1 - T^(w)) dw <= eta [<f, |D| f> + (2/pi) ||f||_2^2],

with |D| the operator of symbol |w|. Take f = A s, s(x) = (Q^2 - x^2)^(1/2) on [-Q, Q]. The classical identity |D| s = 1 on (-Q, Q) (equivalently: the Fourier transform of s is pi Q J_1(Qw)/w, and int_0^infinity J_1(t)^2 dt/t = 1/2) gives <s, |D| s> = int s = pi Q^2/2, and ||s||^2 = 4Q^3/3. Hence

    W(Q) >= sup_A [A pi Q^2 - A^2 eta (pi Q^2/2 + 8Q^3/(3pi))] = (pi Q^2/2)^2 / (eta (pi Q^2/2)(1 + 16Q/(3pi^2))). End.

## 1.3 The edge theorem

[Proved] Theorem 1. For every -1 < Delta <= 0 (eta in (0, pi/2]) and every Q > 0,

    eta Z_Q(Q)^2 >= pi Q / (4(1 + kQ)),   k = 16/(3 pi^2).

Consequently:

(i) (W3') holds with c_Z = inf_Q pi Q/(4(1 + kQ) tanh 2Q) >= 0.373; the clean split (tanh 2Q <= 2Q for Q <= 1/2, tanh 2Q <= 1 for Q >= 1/2) gives c_Z >= (pi/8)/(1 + k/2) = 0.309.

(ii) e_bulk''(nu) <= C_e eta^2 for every nu in [0, 1/2], with C_e = 4 sup_Q h(Q) <= 8.43, h(Q) = tanh(2Q)(1 + kQ)/Q (clean value 4 * 2(1 + k/2) = 10.16).

Proof. Corollary 1.3 and Lemma 1.5 give Z_Q(Q)^2 >= W(Q)/(2Q) >= pi Q/(4 eta (1 + kQ)). For (ii), e_bulk'' = v(Q)/(2Z(Q)^2) with v(Q) <= 2a tanh(2Q) (C145 Lemma B.2, re-derived in 2.1 below) gives e_bulk'' <= a tanh(2Q)/Z(Q)^2 <= (4a eta/pi) h(Q) <= 4 eta^2 h(Q) since a = pi sin eta <= pi eta. At nu = 1/2 (Q = infinity) e_bulk'' is the Q -> infinity limit and the same bound holds. The two numerical sups are over a geometric grid of 4 * 10^5 points on [1e-6, 5], each piece bounded above by monotonicity of its factors (tanh(2Q) and (1+kQ) increasing, 1/Q decreasing), with Q < 1e-6 and Q > 5 bounded in closed form (constants_C146.py). End.

So Wiener-Hopf is not needed. The square-root edge appears here as convexity of W plus the semicircle trial function, which is the minimizer of <f, |D| f> - 2<f, 1> on an interval.

## 1.4 Consistency checks

* Q -> infinity [Proved]: the bound tends to eta Z^2 >= pi/(4k) = 3pi^3/64 = 1.453, below the exact half-filled value pi/2 = 1.571 (Z(infinity)^2 = pi/(2(pi - gamma)), the classical dressed charge at zero field).
* Dilute window [Heuristic for the exact value, Proved for the bound]: with Q = s eta, eta -> 0 at fixed s, the kernel T(eta u) eta du tends to du/(pi(1 + u^2)), and Z tends to 2 pi times the Lieb-Liniger density with c = 1 and Fermi rapidity s. The bound gives Z(s)^2 >= pi s/4. The Lieb-Liniger Luttinger parameter at weak coupling is K = Z(s)^2 ~ pi/sqrt(gamma_LL) with gamma_LL = 4/s^2 (semicircle density), so Z(s)^2 ~ pi s/2, a factor 2 above the bound. This also identifies the limit function of Q/eta seen in C145's numerics: eta Z^2/tanh(2Q) -> pi/4 = 0.785 from above, consistent with the sampled minimum 0.808.
* Q << eta [Proved]: the bound is pi Q/4, which is below the trivial eta Z^2 >= eta; nothing is lost where it matters, since tanh(2Q) is also of order Q.
* Free point eta = pi/2 [Proved]: T = 0, Z = 1, and the bound reads pi/2 >= pi Q/(4(1+kQ)), true since pi Q/(4(1+kQ)) < pi/(4k) = 1.45 < pi/2.

[Numerical, check_W3.json, 7 values of 1 + Delta from 0.5 to 1e-6, 23 values of Q from 0.003 to 4 each, mesh N = 900] Z_Q(Q) is nondecreasing in Q in every case; dW/dQ agrees with 2Z(Q)^2 to 1.4 percent (finite differences on the Q grid); W(Q)/(2Q Z(Q)^2) <= 0.999 (Corollary 1.3); and eta Z(Q)^2 divided by the bound of Theorem 1 is at least 1.578 in every case. The minimum of eta Z(Q)^2/tanh(2Q) is 0.808 at 1 + Delta = 1e-6, against the proved c_Z = 0.373.

# Part 2. The velocity variation (W1)

How (W1) is used. In C145 Theorem A.2 the only appearance of Var(v_Q) is the core term (1/(8L)) sup_(K_0) g Var(v_Q) <= Var(v_Q)/(6L) (g = r/sigma <= 4/3). So a bound Var(v_Q) <= V_W, with V_W independent of eta and Q, is all that is needed; the factor eta in C138's (W1) is not used anywhere else. [Proved, bookkeeping on C145 PROOF Section A.2.]

## 2.1 An exact first-order identity for v

Write phi = -a1' and let G = R[T(. - Q)], H = R[T(. + Q)] (nonnegative functions on [-Q, Q], the restricted resolvent kernel at the two Fermi points).

[Proved] Lemma 2.1. On [-Q, Q],

    v' = F_1 - (r(Q)/r) [ (v(Q) - v) G + (v(Q) + v) H ],
    F_1 = a [ (R a1')^2 - (R a1)(R a1'') ] / (R a1)^2.

Proof. u := eps_Q' satisfies (I - T_Q) u = -a a1' (C136 (5.1)), so u = -a R a1' and v = u/r = -a (R a1')/r. Differentiating f = s + T_Q f and integrating by parts in y (d_x T(x-y) = -d_y T(x-y)) gives, for any f on [-Q, Q],

    (I - T_Q) f' = s' - f(Q) T(. - Q) + f(-Q) T(. + Q).

For f = u (odd, s = -a a1') this is (I - T_Q) u' = -a a1'' - u(Q) J_+, and for f = r (even) it is (I - T_Q) r' = a1' - r(Q) J_0, with J_+ = T(. - Q) + T(. + Q), J_0 = T(. - Q) - T(. + Q). Fix x_0 and the number c = v(x_0). Then D = u' - c r' solves (I - T_Q) D = -a a1'' - c a1' - r(Q)[v(Q) J_+ - c J_0] (u(Q) = v(Q) r(Q)), so

    D(x_0) = -a (R a1'')(x_0) - c (R a1')(x_0) - r(Q)[(v(Q) - c) G + (v(Q) + c) H](x_0),

using R J_+ = G + H and R J_0 = G - H. At x_0, D = u' - v r' = r v', and -c (R a1')(x_0) = c u(x_0)/a = v^2 r/a at x_0. Divide by r(x_0): v' = -a(R a1'')/r + v^2/a - (r(Q)/r)[...], and v^2/a = a (R a1')^2/r^2. End.

[Proved] Lemma 2.2 (properties of G and H). (a) H(x) = G(-x). (b) G >= H on [0, Q]. (c) int_(-Q)^Q G = Z(Q) - 1, hence int_0^Q (G + H) = Z(Q) - 1 and int_0^Q H <= (Z(Q) - 1)/2.

Proof. (a) T_Q commutes with the reflection P f(x) = f(-x), and T(. + Q) = P T(. - Q). (b) G - H = R J_0 with J_0 odd and J_0 >= 0 on [0, Q] (|x - Q| <= |x + Q| there). On odd functions R acts through the reduced kernel T(x - y) - T(x + y) on [0, Q], which is nonnegative for x, y >= 0 (|x - y| <= x + y and T decreasing in |.|); its Neumann series is positive, so R J_0 >= 0 on [0, Q]. (c) int G = <1, R T(. - Q)> = <Z, T(. - Q)> = (T_Q Z)(Q) = Z(Q) - 1; then (a) gives int_0^Q (G + H) = int_(-Q)^Q G, and (b) gives the last claim. End.

[Proved] Lemma 2.3 (lower bound on F_1). For x in [0, Q], with d = Q - x and e(d) = 4 e^(-2bd),

    (F_1)_- <= a min(4, b^2 e(d)(4 + e(d))/(1 - e(d))^2)   (the second entry only when e(d) < 1).

Hence int_0^Q (F_1)_- <= a I(b), I(b) = int_0^infinity min(4, b^2 e(4+e)/(1-e)^2) dd <= I(1) = 5.998 for b in [1, 2].

Proof. First entry. a1 = c_1/(C + kappa) with C = cosh 2x >= 1, kappa = cos eta in [0, 1), and a direct computation gives a1''/a1 = 4(C^2 - kappa C - 2)/(C + kappa)^2. Since -(C + kappa)^2 <= C^2 - kappa C - 2 <= (C + kappa)^2 (the left inequality is 2C^2 + kappa C + kappa^2 >= 2, true for C >= 1), |a1''| <= 4 a1. So F_1 >= -a (R a1'')/r >= -a R(a1'')_+/r >= -4a, using positivity of R.

Second entry. Let delta_k = r_0^(k) - R a1^(k) on [-Q, Q], k = 0, 1, 2. Restricting the differentiated full-line equation r_0^(k) = a1^(k) + T*r_0^(k) to [-Q, Q] gives (I - T_Q) delta_k = g_k, g_k(x) = int_(|u| > Q) T(x - u) r_0^(k)(u) du. From r_0 = c sech(bx), r_0' = -b r_0 tanh(bx) and r_0'' = b^2 r_0 (tanh^2 - sech^2)(bx), |r_0^(k)| <= b^k r_0, so |g_k| <= b^k g_0 and |delta_k| <= R|g_k| <= b^k R g_0 = b^k delta_0. By C145 Theorem B.1, 0 <= delta_0 <= e r_0 with e = e(x) = [r_0(2Q - x) + r_0(2Q + x)]/r_0(x) <= e(d). Now N = (R a1')^2 - (R a1)(R a1'') = (r_0' - delta_1)^2 - (r_0 - delta_0)(r_0'' - delta_2), and N_0 = r_0'^2 - r_0 r_0'' = b^2 r_0^2 sech^2(bx) >= 0. Expanding and bounding each cross term by the estimates above,

    N >= N_0 - 2|r_0'||delta_1| - r_0 |delta_2| - delta_0 |r_0''| - delta_0 |delta_2| >= -b^2 delta_0 (4 r_0 + delta_0) >= -b^2 r_0^2 e(4 + e),

and r = r_0 - delta_0 >= (1 - e) r_0. So F_1 = a N/r^2 >= -a b^2 e(4+e)/(1-e)^2, which is increasing in e, so e may be replaced by e(d).

The integral. Substituting s = bd, I(b) = (1/b)[4 s_* + (b^2/2) int_0^(e_*) (4 + e)/(1 - e)^2 de], where e_* solves b^2 e(4 + e) = 4(1 - e)^2 and s_* = (1/2) log(4/e_*), and int_0^(e_*) (4+e)/(1-e)^2 de = 5/(1 - e_*) - 5 + log(1 - e_*). On b in [1, 2] this is maximal at b = 1, I(1) = 5.9982, and I(2) = 3.996 (constants_C146.json; the closed form agrees with direct quadrature to 1e-9). End.

## 2.2 The variation bound

[Proved] Theorem 2. For every eta in (0, pi/2] and Q in (0, infinity],

    Var_[-Q,Q](v_Q) <= 4a tanh(2Q) + 4a I(b) + 12 a tanh(2Q) (Z_Q(Q) - 1)
                   <= V_W(eta) := a(4 I(b) - 8) + 12 a B <= 92.9,

and V_W(eta) -> 6 pi^2 = 59.2 as eta -> 0.

Proof. Q = infinity: v = ab tanh(bx), Var = 2ab <= 4a. Finite Q: v is odd with v(0) = 0, so Var_[-Q,Q] v = 2 Var_[0,Q] v = 2[v(Q) + 2 int_0^Q (v')_-] (because int |v'| = int v'_+ + int v'_- and int v' = v(Q)). By Lemma 2.1, (v')_- <= (F_1)_- + E_+ with E = (r(Q)/r)[(v(Q) - v) G + (v(Q) + v) H]. On [0, Q]: r(Q)/r <= 1 (r nonincreasing), G, H >= 0, 0 <= v <= M := 2a tanh(2Q) (C145 Lemma B.2: |a1'| <= 2 tanh(2|x|) a1 and positivity of R give |u| <= 2a tanh(2Q) r), so E_+ <= M(G + 2H). By Lemma 2.2(c), int_0^Q (G + 2H) <= (3/2)(Z(Q) - 1). With Lemma 2.3 and v(Q) <= M,

    Var_[-Q,Q] v <= 2[M + 2a I(b) + 3M(Z(Q) - 1)],

which is the first line. For the second use tanh <= 1 and Z(Q) <= B. Finally a B = pi^2 sin(eta)/(2eta) and the maximum over eta of pi sin(eta)(4 I_sup - 8) + 6 pi^2 sin(eta)/eta with I_sup = 5.9982 is 92.87, at eta = 1.16 (constants_C146.json). End.

Remarks. (1) The edge terms G, H carry the square-root edge singularity of r'/r and Z'/Z; they never need to be resolved, only integrated, and their integral is the edge value Z(Q) - 1 itself. The factor a in front makes a(Z(Q) - 1) <= a(B - 1) <= pi^2/2 bounded, which is enough here, although it is not of order eta. (2) F_1 is the log-concavity defect of the translated family t -> R[a1(. + t)](x) at t = 0. At half filling F_1 = a b^2 sech^2(bx) >= 0.

[Numerical, check_W1_*.json, 1 + Delta in {0.5, 0.1, 1e-2, 1e-4, 1e-6}, Q in {0.01, 0.1, 0.5, 1, 2, 4}, mesh N = 700] The identity of Lemma 2.1 holds to a relative error at most 1.2e-3 for 1 + Delta >= 1e-4 and 3.4e-3 at 1 + Delta = 1e-6, Q <= 1 (finite-difference v' against the right side); at 1 + Delta = 1e-6 and Q in {2, 4} the error is 0.8 and 1.3 percent, the mesh-limited cases already flagged in C145. int_0^Q (G + H) equals Z(Q) - 1 to 5e-4 relative (7e-3 in the mesh-limited case). F_1 >= 0 in every sampled case except the mesh-limited one (1 + Delta = 1e-6, Q = 4, min F_1/a = -0.008), so the Turan-type inequality (R a1')^2 >= (R a1)(R a1'') appears to hold, although it is not proved and not used. Var(v)/a <= 3.0, while the first line of Theorem 2 overestimates it by factors 1.02 to 215, the largest near Delta = -1 where Z(Q) - 1 is of order 1/eta (the true edge term is at most 0.42 a).

[Assumed, not needed] (W1') itself, the monotonicity of v_Q, remains unproved. Lemma 2.1 shows what it would take: F_1 >= E pointwise.

# Part 3. The shell inequality (G)

## 3.1 What is proved

[Proved] Proposition 3.1 (the derivative half from the comparability half). Suppose that on the shell S = {x : M_0 <= eta L r_Q(x) < M_1}, M_1 = 4C_D(eta), M_0 = 2, one has sigma_L >= c r_Q with some c in (0, 1]. Then

    eta Var_S(r_Q/sigma_L) <= (1/c)[2 eta log(M_1/M_0) + eta kappa_eta |S|],
    eta kappa_eta <= sqrt(2) max(eta/sin eta, 2 eta) <= 4.45,
    |S| <= 2 min_(d >= 0.7) [d + (1/b) log(2M_1/((1 - 4e^(-2d)) M_0))].

Proof. Derivative of sigma. sigma_L = a1 + (1/L) sum_k T(x - x_k), |a1'|/a1 = 2 sinh 2x/(cosh 2x + cos eta) <= 2 and |T'|/T = sinh(2x)/(sinh^2 x + sin^2 eta). With t = sinh|x| and s = sin eta: for t <= 1, sinh 2|x| = 2t(1 + t^2)^(1/2) <= 2^(1/2) 2t and 2t/(t^2 + s^2) <= 1/s; for t >= 1, (1 + t^2)^(1/2) <= 2^(1/2) t, so the ratio is at most 2^(3/2). Hence |sigma_L'| <= kappa_eta sigma_L with kappa_eta = 2^(1/2) max(1/sin eta, 2). Variation. g = r/sigma satisfies |g'| <= g(|r'|/r + |sigma'|/sigma) <= (1/c)(|r'|/r + kappa_eta). On each half of S (an interval, since S is the difference of two superlevel sets of the even, unimodal r), int |r'|/r = log(r_max/r_min) <= log(M_1/M_0). Length. On [0, Q], points of S at distance at least d from Q satisfy M_0/(eta L) <= r <= r_0 and r_0 <= r/kappa(d) < M_1/(kappa(d) eta L), kappa(d) = 1 - 4e^(-2d) (C145 Theorem B.1); the set of x >= 0 with r_0(x) in [alpha, beta) has length at most (1/b) log(2beta/alpha), from log y <= arccosh y <= log 2y. Add d for the points within d of Q, and double for the two sides. End.

So the derivative half is not where the difficulty sits. What remains of (G) is

[Assumed] (G_c) sigma_L(x) >= (3/4) r_Q(x) for every x in [-Q, Q] with 2 <= eta L r_Q(x) < 4C_D(eta), uniformly in eta, L, n.

The price of Proposition 3.1 is the shell length |S|, which is of order log(C_D(eta))/b, that is of order log log(1/eta). The constant it gives is C_GS(eta) = 17.7, 19.2, 20.4, 21.3, 22.7, 24.1, 25.1 at eta = 1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-9, 1e-12 (constants_C146.json, PropG1_table, with M_1 replaced by its analytic majorant (16/pi)(log(1/eta) + pi/2)).

[Numerical, from C145 check_shell*.json] (G_c) holds with min sigma/r = 0.995 against the needed 0.75 in 36 cases, L to 16384, 1 + Delta to 1e-6, and sigma >= r in all but one mesh-limited case.

## 3.2 Why (G_c) was not proved here

[Proved, the computation behind the obstruction] Comparability on S is a lower bound on sigma - r = R(T*R_L) = int calR_Q(x, .) dR_L. The first-order bound (C142 Lemma II.1) charges every cell (1/L) osc calR_Q(x, .), and the cells at distance u in (eta, 1), where calR'(u) is about -1/(pi eta u), cost (1/(pi eta L)) log(1/eta). Writing the cell error of a smooth f as e_j[f] = -f'(x_j) mu_j + (second order), with mu_j = int over X_j of (y - x_j) sigma_L dy the first moment of the cell about its root (the root is exactly the sigma-median of its cell, since z_L(x_j) = I_j/L), the logarithm is carried by the sum over sparse cells of calR'(x - x_j) mu_j. The trivial bound |mu_j| <= |X_j|/(2L) reproduces the log. A local bound |mu_j| <= C (kappa_eta |X_j|)|X_j|/L holds on any cell where sigma varies by a bounded factor, and gives a log log in place of the log, by the same mechanism as Proposition 3.1. [Heuristic] The true size is mu_j of order sigma'|X_j|^3, which would remove the log entirely, but showing that the background density near a sparse root has a gradient of order sigma rather than sigma/eta is a statement about the spacing regularity of the sparse roots, that is about their log-gas structure, and it was not proved. Nothing in the positivity structure (R >= 0, T >= 0, unimodality) gives a sign for sigma - r: the local quadrature error T*R_L is positive at sparse roots and negative between them, and only its low-frequency content matters.

# Part 4. Assembly

## 4.1 The strip

[Proved, conditional on (G) with M_0 = 2, C_GS = 1] For every -1 < Delta <= 0, L >= 2 and 0 <= n <= L/2,

    |F_L(n) - L e_bulk(Delta, n/L)| <= K/L,   K = 154.0,

with the five terms of C145 Theorem A.2: core (pi/4)(C_G1 + C_GS) = 11.15 (C_G1 = 13.19 analytic, C145 A.3), velocity variation V_W/6 = 15.48 (Theorem 2 here, replacing pi C_W/12), first-order tail 4 pi M_0 A_tail = 100.25, straddling cells 16 pi/3 = 16.76, mismatch pi^2 C_e/8 = 10.39 (Theorem 1 here, C_e = 8.43). By matched doubling (C138 (4.4), C136 Section 9, refereed),

    L(E_(2L) - 2E_L) <= C_strip M,   C_strip = 5K + pi^2/4 = 772.5 <= 773,

for every bond of the U(1) class with -1 < J_z/J <= 0, every field and DM angle, every L >= 2 (constants_C146.json, C_strip). With C_GS = 5 the same formula gives 788.3.

Proof. C145 Theorem A.2 with (W1) in the form Var(v_Q) <= V_W (Part 2) and (W3) in the form e_bulk'' <= C_e eta^2 (Part 1). The only change in the proof of Theorem A.2 is the core term, where (1/(8L))(4/3) Var(v_Q) <= V_W/(6L). End.

[Proved, conditional on (G_c) alone] The same holds with C_GS replaced by C_GS(eta) of Proposition 3.1 (c = 3/4): C_strip(eta) = 838, 844, 849, 852, 858, 863, 867 at eta = 1e-1, ..., 1e-12, growing like (5pi/4)(4/3) 2^(3/2) log log(1/eta). This is length-uniform at every eta and diverges only doubly logarithmically as Delta -> -1, but it is not a uniform constant.

## 4.2 The whole U(1) class J_z/J in [-1, 1]

| region | C in L(E_2L - 2E_L) <= C M | source | status |
| --- | --- | --- | --- |
| 0 <= J_z/J <= 1 | 4078 | C142 Corollary I.10 | proved, not refereed |
| -1 + delta <= J_z/J <= 0 | 6 * 10^15 delta^(-25/3) | C136 | proved, refereed |
| J_z/J = -1 | 0 | C129 Section 8 | proved |
| -1 < J_z/J <= 0 | 773 | Section 4.1 | conditional on (G) only |
| -1 < J_z/J <= 0 | 868 at eta >= 1e-12, order log log(1/eta) | Section 4.1 | conditional on (G_c) only |

Seams as in C145 E.2. [Proved, conditional on (G)] For the whole class, every field, DM angle and L >= 2, L(E_(2L) - 2E_L) <= 4078 M. [Proved, unconditional] Nothing on the strip itself beyond C145 (no new unconditional constant for fixed delta).

## 4.3 Consistency checks

* Free point Delta = 0 [Proved]: T = 0, Z = 1, G = H = 0, r = a1 = 1/(pi cosh 2x), v = 2a tanh 2x, and the identity of Lemma 2.1 reduces to v' = a[(a1')^2 - a1 a1'']/a1^2 = -a(log a1)'' = 4a sech^2(2x) >= 0, so v is monotone and Var = 4a tanh 2Q, inside Theorem 2.
* Half filling [Proved]: F_1 = a b^2 sech^2(bx), the edge terms vanish, Var = 2ab.
* Dilute limit [Proved]: Theorem 2's first line is of order aQ as Q -> 0 (tanh(2Q) and the region length both vanish, Z(Q) - 1 -> 0), matching Var(v) = 2v(Q) <= 4a tanh 2Q.
* Delta -> -1 [Numerical]: Var(v)/a -> 2 at Q >= 2 (C138's Var(v)/eta about 2pi); eta Z^2/tanh 2Q -> limit near 0.8 above pi/4, the dilute Lieb-Liniger value.
