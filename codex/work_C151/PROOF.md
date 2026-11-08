# C151 proof notes: the ferromagnetic strip through the Destri-de Vega representation

This run was done by a Claude Opus agent (2026-11-07), not by Codex. J = 1 throughout. Tags: [Proved] is a mathematical deduction, [Proved, cited] uses an earlier run (named, with its referee status), [Numerical] is floating point without enclosures, [Heuristic] has an unproved step, [Assumed] is a hypothesis.

Notation is that of C136, C138, C142 Part II and C145. eta = pi - gamma, Delta = -cos eta, a = pi sin eta, b = pi/gamma in [1, 2), q = 1 - 2 eta/pi, r_0 = 1/(2 gamma cosh bx), r = r_Q, Z_Q, eps_Q = (I - T_Q)^(-1)(e - h), Q = Q_L the counting boundary, lambda(x) = eta L r_Q(x), C_D(eta) <= (4/pi)(log(1/eta) + pi/2) the constant of C142 Theorem II.6, M_1(eta) = 4 C_D(eta), kappa(d) = 1 - 4e^(-2d).

Summary.

- [Proved] Part 1. An exact Destri-de Vega (DdV) representation of every root sum on [-Q_L, Q_L], at every filling, with all hypotheses verified: Re sigma_L > 0 in the strip |Im z| < min(eta, pi/4), so the counting function has no complex level points and 1 + (-1)^n e^(2 pi i theta) vanishes in the closed rectangle only at the Bethe roots. The vertical sides of the contour sit at the counting boundary, where 1 + (-1)^n e^(2 pi i theta) = 2. Harnack's inequality for the positive harmonic function Re sigma_L turns the real-axis density at a point into exponential decay of the nonlinear term above that point.
- [Proved] Part 2. Sector theorem. If sigma_L >= c r_Q on the shell S = {M_0 <= lambda < 4C_D(eta)} for some c in (0, 3/4], then |F_L(n) - L e_bulk| <= K(c, M_0)/L for all eta in (0, pi/4], L, n. With c = 3/4 and M_0 = 2, which is exactly C146's (G_c), K = 76.7 and C = 5K + pi^2/4 = 386. C146 got only a log log(1/eta) constant from (G_c), and needed the full (G) for a uniform one. No derivative bound, no (D1), no velocity variation (W1) and no variation of r/sigma is used. Any c > 0 closes the strip (K about 50/c for small c).
- [Proved] Part 2, unconditional by-product. Without any hypothesis, |F_L(n) - L e_bulk| <= K_0(eta)/L with K_0 = 202, 278, 426, 570 at eta = 1e-1, 1e-2, 1e-4, 1e-6, growing like 31 log(1/eta). C142 had about 247 log(1/eta), conditional on (W1)-(W3).
- [Numerical] Parts 3, 4. The identity holds to 1e-9 relative at every height and filling tested; the Fermi-edge pieces reproduce the CFT term -pi v_F/(6L); at the free point L(F - L e) -> -pi/6; the DdV nonlinear integral equation, solved by continuation in eta without Bethe input, reproduces the Bethe counting function to 1e-13 and the Bethe energy to 1e-9 down to eta = 0.01; the Harnack bound holds with ratio >= 1.53; the shell comparability holds with min sigma/r = 0.99 against any c > 0.
- [Assumed] Part 5. The shell lower bound itself is not proved. The DdV representation does not supply it, for three identified reasons, the last of which is the same spacing-regularity statement that C145 and C146 isolated.

# Part 1. The exact representation

## 1.1 Kernels in the strip

[Proved] Lemma 1.1. For complex z,

    T(z)  = sin(2eta) / (pi (cosh 2z - cos 2eta)) = sin(2eta) / (2 pi sinh(z - i eta) sinh(z + i eta)),
    a1(z) = sin(eta)  / (pi (cosh 2z + cos eta)),   e(z) = -a a1(z).

T has its poles at z = +-i eta (mod i pi), a1 at z = +-i gamma/2, so sigma_L = a1 + (1/L) sum_k T(. - x_k) is analytic in |Im z| < eta. With z = x + is,

    Re T(z)  = (sin 2eta/pi) (cosh 2x cos 2s - cos 2eta) / |cosh 2z - cos 2eta|^2,
    Re a1(z) = (sin eta/pi)  (cosh 2x cos 2s + cos eta) / |cosh 2z + cos eta|^2,

because Re cosh(2x + 2is) = cosh 2x cos 2s. If |s| < eta_* := min(eta, pi/4), then cosh 2x cos 2s >= cos 2s > cos 2eta and cos 2s > 0, so both real parts are strictly positive. Hence

    Re sigma_L(z) > 0   on the strip  Sigma = {|Im z| < eta_*}.                                  (1.1)

Proof. The first form of T is C138 (1.1) rewritten with sinh^2 x + sin^2 eta = (cosh 2x - cos 2eta)/2; the second uses sinh(z - i eta) sinh(z + i eta) = (cosh 2z - cos 2eta)/2. The form of a1 is C138's tan(eta/2) sech^2 x/[pi(1 + tan^2(eta/2) tanh^2 x)] = sin eta/(2 pi (sinh^2 x + cos^2(eta/2))). e = -a a1 follows from e(x) = -sin^2 eta/(2(sinh^2 x + cos^2(eta/2))). End.

[Proved] Lemma 1.2 (strip bounds). Let 0 <= s < eta <= pi/4 and x real.

(i) |T(x + is)| <= kappa_T(s) T(x) with kappa_T(s) = sin^2 eta/(sin^2 eta - sin^2 s).
(ii) |a1'(x + is)| <= (2/cos^2 2s) a1(x).

Proof. (i) From the sinh form, |T(x + is)| = sin 2eta / (2 pi [(A + alpha)(A + beta)]^(1/2)) with A = sinh^2 x, alpha = sin^2(eta - s), beta = sin^2(eta + s), and T(x) = sin 2eta/(2 pi (A + c)), c = sin^2 eta. So the squared ratio is g(A) = (A + c)^2/((A + alpha)(A + beta)). A direct computation gives g'/g = [(alpha + beta - 2c)A + 2 alpha beta - c(alpha + beta)]/[(A + c)(A + alpha)(A + beta)]. Here alpha + beta - 2c = cos 2eta (1 - cos 2s) >= 0, so the numerator is affine and nondecreasing in A. Hence g is either increasing on [0, infinity) (sup = g(infinity) = 1) or first decreasing and then increasing toward 1 (sup = max(g(0), 1)). Since alpha beta = (sin^2 eta - sin^2 s)^2 <= c^2, g(0) >= 1, and in both cases sup g = g(0) = kappa_T(s)^2. (ii) a1'(z) = -2 sin eta sinh 2z/(pi (cosh 2z + cos eta)^2), |sinh 2z| <= cosh 2x and |cosh 2z + cos eta| >= cosh 2x cos 2s + cos eta >= cos 2s (cosh 2x + cos eta), since cos eta >= 0. Divide by a1(x) = sin eta/(pi(cosh 2x + cos eta)) and use cosh 2x <= cosh 2x + cos eta. End.

[Numerical, constants_C151.py] On a grid of 24001 points in x and five values of eta in [1e-3, pi/4], the ratio of the left side of (i) to its bound is at most 1.0000, of (ii) at most 0.99999, and both real parts in Lemma 1.1 are positive for s up to 0.999 eta.

## 1.2 The counting function in the strip

Let x_1 < ... < x_n be the real ground-state roots (C136 Section 1, refereed: unique, symmetric under x -> -x), I_j = j - (n+1)/2. Define, for z in Sigma,

    theta(z) = (L/(2 pi)) p(z) + sum_k Theta(z - x_k),   p(z) = 2 atan(tan(eta/2) tanh z),   Theta(u) = (1/pi) atan(tanh u / tan eta).

[Proved] Lemma 1.3. theta is analytic on Sigma, real and odd on R, theta' = L sigma_L, theta(x_j) = I_j, and for 0 < s < eta_*,

    Im theta(x + is) = L int_0^s Re sigma_L(x + it) dt > 0.                                       (1.2)

The counting boundary Q = Q_L is the point with theta(Q) = n/2 (C138: z_L(Q_L) = nu/2; Q = infinity at even half filling). In the closed rectangle R_eps = [-Q, Q] x [-eps, eps], 0 < eps < eta_*, the function phi(z) = 1 + (-1)^n e^(2 pi i theta(z)) vanishes exactly at the n roots, each zero is simple, and phi(+-Q) = 2.

Proof. Theta is the antiderivative of T vanishing at 0 (on R, Theta = -theta_C136/(2 pi) with theta_C136(y) = -2 atan(tanh y/tan eta)), and (L/2pi)p' = L a1, so theta' = L sigma_L and theta = L z_L with the C136 counting function; theta(x_j) = I_j is the Bethe equation in C136's form. The arctangents have no branch point in Sigma: tanh(u)/tan eta is purely imaginary only for Re u = 0, where it equals i tan(Im u)/tan eta of modulus < 1, and similarly for p. (1.2) is d/ds Im F(x + is) = Re F'(x + is) for analytic F with Im F = 0 on R, together with (1.1). Zeros: |(-1)^n e^(2 pi i theta)| = e^(-2 pi Im theta) is not 1 off the real axis by (1.2) and its mirror, so zeros are real, at points where theta is in (n+1)/2 + Z. Since theta increases from -n/2 to n/2 on [-Q, Q], these are exactly I_1, ..., I_n (the endpoints +-n/2 have the other parity), simple because theta' = L sigma_L > 0. At +-Q, (-1)^n e^(+-i pi n) = 1. End.

## 1.3 The contour identity

Let ell(z) = Log(1 + w(z)), w = (-1)^n e^(2 pi i theta), principal branch, for 0 < Im z < eta_*; there |w| < 1, so ell is analytic. On the real axis ell(x + i0) is the continuous boundary value, with logarithmic singularities at the roots and Im ell(x + i0) = arg(1 + e^(i 2 pi theta(x) + i pi n)) in [-pi/2, pi/2].

[Proved] Lemma 1.4 (DdV identity). Let 0 < eps < eta_* and let f be analytic on a neighbourhood of R_eps and real on R. Let Gamma be a piecewise smooth path from Q to -Q inside {|Re z| <= Q, 0 <= Im z <= eps} whose real parts are finitely many closed intervals. Then

    sum_j f(x_j) - L int_(-Q)^Q f sigma_L dx = -(1/pi) Im int_Gamma f'(z) ell(z) dz.                (1.3)

At even half filling (Q = infinity) the same holds for f with f, f' exponentially decaying in the strip and Gamma running from +infinity to -infinity, for instance f = eps_inf = -a r_0.

Proof. By Lemma 1.3 and the argument principle, sum_j f(x_j) = (1/(2 pi i)) contour integral of f phi'/phi around R_eps, counterclockwise. On the upper part (Im z > 0) phi'/phi = ell'. On the lower part write phi = w (1 + 1/w), so phi'/phi = 2 pi i L sigma_L + ellt', with ellt = Log(1 + 1/w) analytic for -eta_* < Im z < 0. Split the boundary at the real points +-Q, where phi = 2 and both logarithms equal log 2. The upper path P_+ (Q -> Q + i eps -> -Q + i eps -> -Q) gives int f ell' = [f ell] - int f' ell, and the lower path P_- (-Q -> -Q - i eps -> Q - i eps -> Q) gives 2 pi i L int_(-Q)^Q f sigma_L (the term f theta' is analytic, deform to R) plus int f ellt' = [f ellt] - int f' ellt. The boundary terms are f(-Q) log 2 - f(Q) log 2 and f(Q) log 2 - f(-Q) log 2 and cancel. Since f and theta are real on R, ellt(zbar) = conj(ell(z)) and f'(zbar) = conj f'(z), and P_- is the mirror of P_+ run backwards, so int_(P_-) f' ellt = -conj(int_(P_+) f' ell). Hence sum_j f(x_j) - L int f sigma_L = -(1/(2 pi i)) [A - conj A] = -(1/pi) Im A, A = int_(P_+) f' ell. Finally P_+ may be replaced by Gamma: f' ell is analytic between them and continuous up to the real axis except at the roots, where |f' ell| <= C(1 + |log|z - x_j||) is integrable, so Cauchy's theorem applies after indenting by small half-discs whose contribution vanishes. At Q = infinity, ell -> log 2 as |Re z| -> infinity (theta(+-infinity) = +-n/2 at even half filling, since the vacancy count L(eta/pi + q/2) equals n) and the vertical sides at +-X vanish as X -> infinity. End.

[Proved, interpretation] On a real-axis stretch of Gamma, -(1/pi) Im(f' ell) = -(1/pi) f'(x) arg(1 + w(x)) is the Euler-Maclaurin sawtooth term, bounded by |f'(x)|/2. This is half the oscillation charge of C138 Lemma 1 and C142 Lemma II.1. On a lifted stretch the integrand is bounded by |f'| G(Psi), with Psi = -log|w| (Lemma 1.6), which is exponentially small where the roots are dense.

## 1.4 Harnack's inequality for the counting function

Write Psi_s(x) = 2 pi Im theta(x + is) = -log|w(x + is)|.

[Proved] Lemma 1.5. Let eta <= pi/4 (so eta_* = eta), x real and 0 <= s <= eta/2. Then

    Psi_s(x) >= 2 pi eta L sigma_L(x) h(s/eta),   h(t) = 2 log(1 + t) - t,
    Psi_(eta/2)(x) >= alpha_H eta L sigma_L(x),   alpha_H = 2 pi h(1/2) = 2 pi (2 log 1.5 - 1/2) = 1.95363,
    Psi_s(x) >= 4 pi h(1/2) L sigma_L(x) s = 3.90726 L sigma_L(x) s.

Proof. u = Re sigma_L is harmonic and, by (1.1), positive on the open disc of radius eta about x. Harnack's inequality gives u(x + it) >= ((eta - t)/(eta + t)) u(x) for 0 <= t < eta. Insert into (1.2): Psi_s >= 2 pi L sigma_L(x) int_0^s (eta - t)/(eta + t) dt = 2 pi L sigma_L(x) eta h(s/eta). The last line uses concavity of h, h(t) >= 2h(1/2)t on [0, 1/2]. End.

So the decay of the nonlinear term above x is controlled by the density at the same real point. No smoothness of sigma_L, no comparison with r_Q, and no information on neighbouring roots is used. (The true value is close to pi eta L sigma, Part 3.)

[Proved] Lemma 1.6. For |w| < 1, |Log(1 + w)| <= G(-log|w|), G(t) = -log(1 - e^(-t)) + arcsin(e^(-t)), G is decreasing, and int_0^infinity G = pi^2/6 + (pi/2) log 2 = 2.73373.

Proof. 1 - |w| <= |1 + w| <= 1 + |w| gives |log|1 + w|| <= -log(1 - |w|). The disc of radius |w| about 1 subtends the half-angle arcsin|w| at 0, so |arg(1 + w)| <= arcsin|w|. int_0^inf -log(1 - e^(-t)) dt = sum_m 1/m^2 = pi^2/6 and int_0^inf arcsin(e^(-t)) dt = int_0^1 arcsin(v)/v dv = (pi/2) log 2. End.

## 1.5 The dressed energy in the strip

[Proved] Lemma 1.7. Let eta <= pi/4 and 0 <= s <= eta/2. The function eps_ext(z) = e(z) - h + int_(-Q)^Q T(z - y) eps_Q(y) dy is analytic on |Im z| < eta, equals eps_Q on [-Q, Q], is real on R, and for x in [-Q, Q]

    |eps_ext'(x + is)| <= 2a kappa_eps r_Q(x),   kappa_eps = max(kappa_T(eta/2), 1/cos^2 eta) <= 2,

with kappa_eps -> 4/3 as eta -> 0. At half filling the same holds for eps_inf = -a r_0 on R.

Proof. Analyticity: the integrand is analytic in z for |Im z| < eta, uniformly on the compact y-interval. eps_Q(+-Q) = 0 (C136 Section 5, refereed) and eps_Q is Lipschitz, so integrating by parts, eps_ext'(z) = e'(z) + int_(-Q)^Q T(z - y) eps_Q'(y) dy. By Lemma 1.2(i) and |eps_Q'| <= 2a r_Q (C136 (5.2), refereed), the integral is at most 2a kappa_T int T(x - y) r_Q(y) dy = 2a kappa_T (r_Q - a1)(x), using r_Q = a1 + T_Q r_Q on [-Q, Q]. By Lemma 1.2(ii), |e'(z)| = a|a1'(z)| <= (2a/cos^2 2s) a1(x). Add, use a1 <= r_Q, and note kappa_T(s), 1/cos^2 2s are increasing in s. For eta <= pi/4, kappa_T(eta/2) = 1/(1 - 1/(4cos^2(eta/2))) <= 1.4142 and 1/cos^2 eta <= 2. At Q = infinity the proof is the same with eps_inf - T*eps_inf = e on R and |eps_inf'| = a|r_0'| <= ab r_0 <= 2a r_0. End.

# Part 2. The sector theorem

## 2.1 Statement

[Assumed] (S_c), for given c in (0, 3/4] and M_0 > 0: sigma_L(x) >= c r_Q(x) for every x in [-Q_L, Q_L] with M_0 <= lambda(x) < M_1(eta) = 4C_D(eta).

On K_1 = {lambda >= M_1(eta)} the same holds unconditionally with c = 3/4 by C142 Theorem II.6 (refereed correct in C147), since |sigma - r| <= C_D(eta)/(eta L) <= r/4 there. So under (S_c), sigma_L >= c r_Q on all of K_0 = {lambda >= M_0}. (S_{3/4}) with M_0 = 2 is C146's (G_c) word for word.

[Proved, conditional on (S_c)] Theorem 2.1. For every eta in (0, pi/4], L >= 2 and 0 <= n <= L/2,

    |F_L(n) - L e_bulk(Delta, n/L)| <= K(c, M_0)/L,
    K = 2 pi M_0 A_tail + 2 a_0/c + 4 a_0 kappa_eps I_G/(3.90726 pi c) + 4 kappa_eps ell_B sum_(k>=0) 2^(k+1) M_0 G(alpha_H c 2^k M_0) + pi^2 C_e/8,

with a_0 = pi sin(pi/4) = 2.2214 (a bound for a on (0, pi/4]), kappa_eps = 2, I_G = 2.73373, alpha_H = 1.95363, A_tail = min_d [d + pi/(2 kappa(d))] = 3.4613 (d = 1.492), ell_B = min_d [d + log(4/kappa(d))] = 3.0342 (d = 1.242), and C_e = 8.43 (C146 Theorem 1, proved, not refereed).

[Numerical evaluation, constants_C151.json] c = 3/4, M_0 = 2: K = 76.72 (tail 43.50, straddling cells 5.92, connectors 5.28, lifted core 11.62, mismatch 10.40), C = 5K + pi^2/4 = 386.1. Optimizing M_0 at c = 3/4 (M_0 = 1.55, a slightly larger hypothesis region) gives K = 74.57, C = 375.3. For smaller c: K = 106.7 (c = 1/2), 202.9 (c = 1/4), 491.7 (c = 1/10), at the optimal M_0.

## 2.2 Proof

Fix eta, L, n; write r = r_Q, f = eps_ext (eps_inf at even half filling), eps = eta/2 (< eta_* = eta). By C136 (4.3) (refereed) and Lemma 1.4,

    F_L(n) - L[E(Q) + h(Q)(nu - m(Q))] = sum_j f(x_j) - L int_(-Q)^Q f sigma_L = -(1/pi) Im int_Gamma f' ell,   (2.1)

for any admissible Gamma. K_0 = {lambda >= M_0} is a centred interval [-x_0, x_0] (r is even and nonincreasing on [0, Q], C136 Section 2). Cell boundaries are the points where theta is in n/2 + Z; they include +-Q, they are symmetric, and consecutive ones enclose one root and sigma-mass 1/L.

Choice of Gamma. Case A: K_0 is empty, or contains no cell boundary other than possibly 0 and lies inside one cell. Take Gamma = [-Q, Q] run from Q to -Q. Case B: x_0 < Q and there is a cell boundary in (0, x_0]; let x_0' be the largest one, and take Gamma: Q -> x_0' along R, up to x_0' + i eps, left to -x_0' + i eps, down to -x_0', along R to -Q. (If n is even, 0 is a cell boundary, and x_0' = 0 is allowed; the lifted part is then empty.) Case C: x_0 = Q, that is lambda(Q) >= M_0. Take Gamma: Q -> Q + i eps -> -Q + i eps -> -Q.

Real-axis pieces (Cases A, B). By the remark after Lemma 1.4 and |eps_Q'| <= 2a r,

    (1/pi)|Im int_real f' ell| <= (1/2) int_real |f'| <= a int_real r.

The real pieces are [x_0', Q] and its mirror. On (x_0, Q], r < rho_0 = M_0/(eta L). Per side, the part within d of Q has length at most d; on the rest r >= kappa(d) r_0 (C145 Theorem B.1, refereed correct in C147), so r_0 < rho_0/kappa(d) there and its r-integral is at most int_(r_0 < rho_0/kappa) r_0 <= (pi/(2b)) rho_0/kappa(d), by int_y^inf r_0 <= (gamma/2) r_0(y) (C142 Lemma I.1, refereed). With b >= 1, two sides contribute at most 2 rho_0 A_tail, and a * 2 rho_0 A_tail <= 2 pi M_0 A_tail/L since a/eta <= pi. The pieces [x_0', x_0] and its mirror lie in single cells inside K_0, where r <= sigma/c, so their r-mass is at most 1/(cL) each: total 2a/(cL). In Case A the whole interval is charged in the same way (if K_0 is inside one cell its r-mass is at most 1/(cL)). At Q = infinity the tail bound is smaller (no d-part).

Vertical pieces (connectors at +-x_0' in Case B, Fermi-edge pieces at +-Q in Case C). The base point y = x_0' or Q is a cell boundary in K_0, so w(y) = 1 and sigma(y) >= c r(y). By Lemmas 1.7, 1.6 and 1.5, each piece is at most

    (1/pi) int_0^eps |f'(y + is)| G(Psi_s(y)) ds <= (1/pi) 2a kappa_eps r(y) int_0^inf G(3.90726 L sigma(y) s) ds
                                                 = (2a kappa_eps/pi) (r(y)/sigma(y)) I_G/(3.90726 L) <= 2a kappa_eps I_G/(3.90726 pi c L).

Lifted piece (Cases B, C). On x + i eps with |x| <= x_0' (or Q), x is in K_0, so by Lemmas 1.7, 1.6, 1.5 and sigma >= c r,

    (1/pi) int |f'(x + i eps)| |ell(x + i eps)| dx <= (2a kappa_eps/(pi eta L)) int_(K_0) lambda(x) G(alpha_H c lambda(x)) dx.

Split [0, x_0] into the bands B_k = {2^k M_0 <= lambda < 2^(k+1) M_0}. Each B_k lies in {rho_k <= r < 2 rho_k}, whose length is at most ell_B: within d of Q at most d, and beyond, kappa(d) r_0 <= r <= r_0 (C136 (2.2), C145 B.1) forces r_0 in [rho_k, 2rho_k/kappa(d)), a set of length at most (1/b) log(4/kappa(d)) (C146 Proposition 3.1, from log y <= arccosh y <= log 2y). On B_k, lambda G(alpha_H c lambda) <= 2^(k+1) M_0 G(alpha_H c 2^k M_0) since G decreases. With a/eta <= pi the lifted piece is at most (4 kappa_eps ell_B/L) sum_k 2^(k+1) M_0 G(alpha_H c 2^k M_0); the series converges geometrically.

Mismatch. By convexity of e_bulk and (W3) in the proved form e_bulk'' <= C_e eta^2 (C146 Theorem 1), 0 <= L[e_bulk(nu) - E(Q) - h(nu - m)] <= C_e eta^2 L(nu - m)^2/2 <= pi^2 C_e/(8L), using L|nu - m| <= B - 1 (C136 (4.2), refereed) and (B - 1)eta <= pi/2 (C138 (4.2)). At even half filling it is zero.

Adding the pieces and bounding a <= a_0 in the straddling and connector terms gives Theorem 2.1. End.

## 2.3 Remarks

[Proved] What is not used. Theorem 2.1 uses no derivative bound on sigma_L - r_Q ((D1) is not needed anywhere), no bound on Var(v_Q) ((W1) is not needed; C146 Theorem 2 is superseded for this purpose), and no variation of r/sigma on the shell (the derivative half of C145's (G) is not needed). The inputs are C136 (2.2), (4.2), (4.3), (5.2), C142 Theorem II.6 and Lemma I.1, C145 Theorem B.1, C146 Theorem 1 and Proposition 3.1 (band length), and the new Lemmas 1.1 to 1.7.

[Proved] Unconditional version. Apply the same proof with K_0 = {lambda >= theta_0 C_D(eta)}, theta_0 > 1, where C142 Theorem II.6 gives eta L sigma >= lambda - C_D >= (1 - 1/theta_0) lambda with no hypothesis. Then |F_L - L e_bulk| <= K_0(eta)/L with K_0(eta) = min over theta_0 of K(1 - 1/theta_0, theta_0 C_D(eta)) (C_D replaced by its majorant (4/pi)(log(1/eta) + pi/2)). [Numerical evaluation] K_0 = 202.1, 278.2, 425.7, 570.0, 994.1 at eta = 1e-1, 1e-2, 1e-4, 1e-6, 1e-12, that is about 31.4 log(1/eta) + 140. This is still a logarithm, but it assumes nothing beyond proved statements (the mismatch uses C146 Theorem 1, proved but not refereed; (W1) and (D1) are not used) and its slope is about one eighth of C142's 247.

[Proved] Weaker shell hypotheses. The proof only uses sigma >= c r through G(alpha_H eta L sigma) on the lifted segment and through r/sigma at one straddling cell and one connector point per side. On the lifted segment it suffices that eta L sigma >= m(lambda) with m nondecreasing and sum_k 2^k G(alpha_H m(2^k M_0)) < infinity, which holds for m(lambda) = mu log lambda with mu > 1/alpha_H = 0.512. So, up to the two boundary cells, a logarithmic lower bound eta L sigma_L >= 0.52(1 + delta) log lambda on the shell already gives a uniform constant. I record this because it shows how little of the comparability is actually used by the energy, not because it is easier to prove (see Part 5).

[Proved] Relation to the other hypotheses. (S_c) for any c > 0 implies C145/C146's (G_c) form with 3/4 replaced by c, and implies C150's (H_I): with the band bound of 2.2, int over {2 <= lambda < 2C_D} of dy/(eta L sigma) <= (1/c) int dy/lambda <= (2 ell_B/c) sum_(k>=0) 1/(2^k * 2) = 2 ell_B/c = 6.07/c. Neither (H_I) nor the logarithmic version above implies the other. All of them follow from shell comparability with any positive constant, and none is proved here.

[Heuristic, Task 3 of the spec] The matched defect in DdV language is L D = Phi(2L)/2 - 2 Phi(L), Phi(L) = L(F_L - L e_bulk), and each Phi is the sum of the four pieces above; bulk terms are already absent from (2.1). On the 2L ring the shell sits where r is half as large, so the two DdV integrals are over different sets and I found no cancellation between them beyond what the sector bound gives (L D <= (5/2) max|Phi|).

# Part 3. Consistency checks [Numerical]

All floating point; scripts and JSON in this folder.

1. Identity at half filling (t1.py, run inline; f = eps_inf). eta in {pi/2 - 1e-9, 0.5, 0.1, 0.01}, L in {16, 64, 256}: theta(x_j) = I_j to 2e-14 (this checks the analytic continuation formula for theta), and the DdV integral at heights eta/8, eta/4, eta/2 equals L(F - L e_bulk(1/2)) from the Bethe energy to 1e-7 absolute, in most cases to 1e-8 (the residual at eta = 0.5 grows slightly with the height and comes from the quadrature panels). Height independence holds to the same accuracy.
2. Free point Delta = 0 (eta = pi/2 - 1e-9). L(F - L e) = -0.52596, -0.52375, -0.52361 at L = 16, 64, 256, tending to -pi/6 = -0.52360, the CFT value with c = 1, v = 1.
3. CFT term at half filling. At eta = 0.5, L(F - L e) = -0.14926 at L = 256 against -pi v/6 = -0.14928 with v = (pi/2) sin(gamma)/gamma.
4. Finite Fermi boundary (check_finite.py, check_finite.json; eta in {0.3, 0.1}, (L, n) in {(64,16), (64,24), (128,48), (256,64)}). f = eps_ext built from the restricted solver (an exactly analytic function). The DdV integral over the rectangle top at heights eta/4 and eta/2 equals the direct sum_j f(x_j) - L int f sigma to all nine printed digits in all eight cases. The C136 identity F - L[E + h(nu - m)] = L int eps_Q dR_L holds to 3e-4 (restricted-solver accuracy). The two Fermi-edge pieces approach -pi v_F/(6L), v_F = eps'(Q)/(2 pi sigma(Q)), when lambda(Q) is large: at eta = 0.3, L = 256, n = 64 (lambda(Q) = 4.29) the edge pieces give L * (-0.0719834), the full integral L * (-0.0719833) and the CFT value L * (-0.0719943). For lambda(Q) of order one the horizontal part is not small, as expected.
5. Harnack bound (check_shell.py, check_shell_1e-*.json). min over the grid of Psi_(eta/2)/(alpha_H eta L sigma) is at least 1.532 in all 46 cases (1 + Delta = 1e-2, 1e-4, 1e-5, 1e-6; L up to 16384; n in {L/8, L/4, 3L/8, L/2 - 2, L/2}). On the shell, Psi/(pi lambda) >= 0.991, so the true decay rate is pi eta L sigma, against the proved alpha_H = 1.954.
6. Shell margins (same files; 34 cases with S nonempty). min sigma/r over S = {2 <= lambda < 4C_D(eta)}: 1.0000 (1e-2, L <= 2048), 0.9964 (1e-4, L <= 4096), 0.9916 (1e-5, L = 8192), 1.0000 (1e-6, L = 16384). The two values below one are at n = L/2 - 2 with Q_L of order 7, where the restricted solver is mesh-limited (C145 saw the same). min of eta L sigma/log(lambda) on S is 2.71 or more. max over [0, Q] of eta L|sigma - r| is at most 0.317, as in C138 and C142.
7. DdV NLIE against Bethe roots (Part 4, nlie.py, nlie_main.json).

# Part 4. The nonlinear integral equation at half filling

[Proved] Derivation. Apply Lemma 1.4 (half-filled form, conjugate pair not combined) to f(y) = Theta(x - y) for x real, then dress with (I - T)^(-1) on R. With calR = sum_(k>=1) T^(*k) (transform That/(1 - That), analytic in |Im u| < eta since That/(1 - That) decays like e^(-eta|w|)), and the odd driving term z_0(x) = (1/pi) atan(tanh(bx/2)) (L z_0' = L r_0),

    theta(x) = L z_0(x) - (1/pi) Im int calR(x - y - i eps) ell(y + i eps) dy,     x real, 0 < eps < eta/2.

Continuing x to x + i eps (allowed since |Im| of every kernel argument stays below 2 eps < eta) gives the closed equation on the line Im z = eps,

    th(x) = L z_0(x + i eps) - (1/(2 pi i)) int calR(x - y)[ell(y) - log 2] dy + (1/(2 pi i)) int calR(x - y + 2 i eps)[conj ell(y) - log 2] dy,   (4.1)

ell = Log(1 + (-1)^n e^(2 pi i th)). Subtracting log 2 is exact because int calR(u) du = int calR(u + 2 i eps) du (contour shift), and it makes the integrands decay, since th -> +-n/2 + i0 at +-infinity. The energy is F - L e_bulk = (1/pi) Im int eps_inf'(x + i eps) ell(x) dx.

[Numerical, nlie.py, nlie_main.json] (4.1) was discretized on a uniform grid (spacing min(eta/8, gamma/(6L))) with FFT convolutions and solved without Bethe input: Picard iteration at eta >= 0.5, and for smaller eta Newton-Krylov with continuation in eta from 0.5 in steps of ratio 1.25. Plain Picard iteration diverges for eta <= 0.1, which reflects the norm pi/(2eta) of the dressing; Newton converges in every case reported below. Results for eta in {pi/2, 0.5, 0.1, 0.03} and L in {16, 64, 256}, and eta = 0.01 at L = 16, 64 (table in RESULT): at L = 256 the solution agrees with theta(x + i eps) computed from the Bethe roots to between 5e-14 and 5e-10, and the DdV energy agrees with the Bethe energy to 1e-9 absolute; at small L and eta <= 0.03 the agreement is 3e-5 to 2e-4 in theta and 3e-6 relative in the energy, limited by grid truncation (the Bethe theta has residual 2e-7 to 5e-7 in the discretized equation). Roots read off the real-axis counting function agree with the Bethe roots to the interpolation accuracy 1e-5 to 3e-7. At eta = 0.01 the continuation in steps of ratio 1.25 failed at L = 64 and succeeded with ratio 1.1 (nlie_eta001.json); L = 256 at eta = 0.01 was not run; the residual of the Bethe theta in the equation was checked only at the final continuation step of each converged run.

# Part 5. What is left, and why the DdV representation does not supply it

[Assumed] The remaining statement is the shell lower bound (S_c) for some c > 0 (or the logarithmic version of 2.3, or C150's (H_I)). [Numerical] Its margin: min sigma/r on S is 0.99 against any positive c, in every case down to 1 + Delta = 1e-6.

[Proved, the computations behind the obstruction] Three attempts to get (S_c) from the representation, and where each stops.

1. Bootstrap through the dressed NLIE. On [-Q, Q], sigma - r = (1/(pi L)) Im int d_x calR_Q(x, y + i eps) ell(y + i eps) dy plus edge terms, from Lemma 1.4 applied to f = calR_Q(x, .). For eta < |x - y| < 1, |d calR(u)| is about 1/(pi eta |u|) (C142 II.4). Bounding |ell| by G(Psi) and Psi by Lemma 1.5 under the hypothesis sigma >= c r gives, at x, a contribution (1/(pi^2 eta L)) int G(alpha_H c lambda(y)) du/|u|, whose part near x is G(alpha_H c lambda(x)) log(1/eta)/(pi^2 eta L). Returning sigma(x) >= c r(x) needs this to be below (1 - c) lambda(x)/(eta L), so the bootstrap closes only where lambda(x) >= (1/(alpha_H c)) log log(1/eta) + O(1). A doubly logarithmic threshold, but not a constant.
2. The interface with the sparse tail. On {lambda < M_0} the lifted contour is useless (Psi = O(1)), and the tail must be charged first order. Its influence on sigma at a point at distance u_t from it is (1/(pi eta L)) log(1/max(u_t, eta)) + O(1/(eta L)), which is the geometric obstruction of C145 D.1 at the inner edge of the hypothesis region: the improved bound does not return the hypothesis within e^(-pi(1-c)M_0) of that edge, a fixed layer in which the cells adjacent to x are again uncovered and cost log(1/eta).
3. Using the oscillation of ell. The log(1/eta) in item 1 is the Hilbert-transform sum over scales eta < u < 1 of a term bounded without its phase. One integration by parts in the phase 2 pi theta(y + i eps), whose frequency is 2 pi L sigma = 2 pi lambda_sigma/eta, gains a factor eta/(lambda u) per scale, but produces the amplitude derivative sigma'(y + i eps)/sigma, and Harnack gives only |sigma'| <= (2/eta) sigma. That returns log(1/eta)/lambda (and e^(-2 pi lambda) log(1/eta) after splitting off the mean). A log-free bound needs |sigma'(y + i eps)| of order sigma on the unit scale in the shell, that is regular spacing of the shell roots, which is C146's diagnosis in another form (the cell first moments mu_j of C146 Section 3.2 are exactly the non-oscillating part of ell over a cell).

[Heuristic] The representation shows why the true answer has no logarithm: for a locally periodic root lattice the mean of ell over a period is zero exactly (the cell integral of Log(1 + e^(2 pi i t)) at constant Im t vanishes), so only gradient corrections survive, and those carry an extra factor eta. What is missing is an a priori statement that excludes locally irregular configurations in the shell. Such a statement has to use the Bethe equations beyond positivity and the counting function, for instance the convexity of the Yang-Yang action that the comparison-principle route works with.

# Part 6. Assembly

[Proved, conditional on (S_{3/4}) with M_0 = 2, that is C146's (G_c) alone] On -1 < J_z/J <= 0 with eta <= pi/4 (Delta <= -cos(pi/4)), every field, DM angle and L >= 2,

    |F_L(n) - L e_bulk| <= 76.72/L,      L(E_(2L) - 2E_L) <= 386.1 M,

by matched doubling and the C129 reductions (C138 (4.4), C136 Section 9, refereed). The rest of the attractive side, -cos(pi/4) <= Delta <= 0, is covered by C146's conditional theorem (constant 773 under (G)) or unconditionally by C136 with delta_0 = 1 - cos(pi/4) (constant of order 1e20). With (S_c) for a smaller c the constant is 5K(c) + pi^2/4, for example 536 at c = 1/2 and 2461 at c = 1/10.

[Proved, unconditional] On eta <= pi/4: L(E_(2L) - 2E_L) <= (5 K_0(eta) + pi^2/4) M with K_0(eta) of order 31.4 log(1/eta) + 140, and the C138 Theorem A bound (15 pi/4) eta L. Their minimum is still unbounded.

| region | C in L(E_2L - 2E_L) <= C M | status |
| --- | --- | --- |
| 0 <= J_z/J <= 1 | 4078 | proved (C142), refereed (C147) |
| -1 + delta <= J_z/J <= 0 | 6e15 delta^(-25/3) | proved (C136), refereed |
| J_z/J = -1 | 0 | proved (C129) |
| -1 < J_z/J <= -cos(pi/4) | 386 | proved here conditional on (G_c) alone; not refereed |
| -1 < J_z/J <= -cos(pi/4) | about 5(31.4 log(1/eta) + 140) | proved here, unconditional; not refereed |
