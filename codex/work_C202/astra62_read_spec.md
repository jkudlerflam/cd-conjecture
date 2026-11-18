# Task: adversarial referee check of a proof (find errors, do not confirm)

You are an expert in Bethe ansatz for the XXZ chain and in linear integral equations. Below is DOCUMENT A, a claimed proof about the easy-axis XXZ ring (Delta = cosh(lambda) > 1), followed by DOCUMENT B, an excerpt of an earlier refereed proof for 0 <= Delta < 1 that DOCUMENT A cites "verbatim" in places (Lemmas I.5 to I.8 and Theorem I.4). You have no other files. Treat DOCUMENT B and the results it cites (C139, C140, C129) as correct in their own setting (the real line, Delta = cos(gamma)). Your job is to decide whether each step of DOCUMENT A is correct in its setting (the circle T = R/(pi Z)).

Work through these in order and give, for each, one of: CORRECT, CORRECT AFTER MINOR REPAIR (state it), GAP (state what is missing), ERROR (give the failing line and if possible a counterexample).

1. Section 1: facts (F1) to (F4). Check signs and the 2 pi shifts of p and theta.
2. Lemma 2.1 (circle convolution preserves circle-unimodality) and Lemma 2.2 (i)-(vii), Lemma 2.3 (i)-(iv). Check each formula, especially r0 = (K_e/pi^2) dn, J = a1 (*) r0, sum_m J^{(*)m} = K, ||K'||_1, M_1 bounds, c_p, kappa_A.
3. Lemma 4.1 (pairing on the circle), both cases of (a) and both bounds of (b). Check carefully that every gap is paired with a distinct root, that the monotonicity of J(x - .) on each gap is as claimed in the antipode case, and what happens when x and its antipode x* lie in the same gap, or when y* lies in out.
4. Corollary 4.2, Lemma 4.3 (including the elliptic-function inequality Loc >= 4 sin(pi/8) k K_e / pi^2), and Theorem 4.4. In Theorem 4.4 check the final bound Phi(lambda) <= lambda [2 coth(lambda) - pi J(0)] / [2 pi sin(pi/8) (1 - 2 e^{-pi^2/lambda})^2], the theta-function inequality theta_2(e^{-lambda}) >= sqrt(pi/lambda)(1 - 2 e^{-pi^2/lambda}), the bound pi J(0) >= ln 2/lambda, the monotonicity in lambda of each factor used on (0,1] and in the cell bracketing on [1, 1.36] (lambda <= l1, coth(lambda) <= coth(l0), lower bound 1/2 + 2 sum_{j<=3} t^j/(1+t^j) evaluated at l1, denominator evaluated at l1), and the direct bound for lambda >= 1.36. Recompute the numbers 0.804, 0.544, 0.939 (at least the endpoint values) and 0.0210.
5. Theorem 4.5: the continuation through Delta = 1 in XXX rapidities, the open-closed argument on lambda in (0, infinity), and the identification E = F_L(n) by the identity theorem. Is completeness of the Bethe ansatz, or nonvanishing of the Bethe vector, used anywhere on Delta > 1? Is anything missing for odd L, or for n = L/2, or n = (L-1)/2?
6. Section 3 lemmas (3.1 to 3.4): are the bounds uniform in Q in (0, pi/2], including Q = pi/2?
7. Sections 5 and 6: Lemma 5.1 (sawtooth with s = 3/2), Lemma 5.2 (g = 13.5 s), Lemma 5.3, the core, the ends, the empty core, the curvature term. Check that the transfer from DOCUMENT B is legitimate on the circle and that sigma_L possibly being negative on (x_n, Q] does no harm.
8. Section 7 arithmetic: recompute every number, in particular whether the pieces sum to at most 218.4 ab, and 218.4 * 2 pi^2 <= 4311, (5/2) * 4311 <= 10778, plus pi^2/4 <= 10781.

Show the computations that support each verdict. Be concrete. If you believe a step is wrong, try to build an explicit counterexample (small n, L). End with a short list of all problems found, most serious first. Plain prose, no filler.

=====================================================================
DOCUMENT A
=====================================================================
# C148 proof: the easy-axis U(1) class, J_z/J > 1

2026-11-07. Run by a Claude Opus agent, not by Codex. J = 1 in all analytic sections (restore J by homogeneity). Tags: [Proved] is a mathematical deduction; [Proved, cited] uses a refereed earlier result (named); [Proved, computer-checked] is a deduction in which one explicit inequality between monotone elementary functions is verified by finite interval bracketing in floating point (script named, margins at least 6 percent, far above rounding); [Numerical] is floating point without enclosures; [Heuristic]; [Assumed].

Notation. Delta = cosh(lambda), lambda > 0, b = pi/lambda, a = pi sinh(lambda), M_0 = (1+Delta)/4, t = e^(-2 lambda), q = e^(-lambda). The circle is T = R/(pi Z), represented by (-pi/2, pi/2]; d(x,y) is circular distance. Convolution on T is (f (*) g)(x) = int_T f(x-y) g(y) dy, Fourier coefficients fhat_k = int_T f e^(-2ikx) dx, so f = (1/pi) sum_k fhat_k e^(2ikx) and (f (*) g)^ = fhat ghat. A function on T is "circle-unimodal" if it is even and nonincreasing in d(x,0).

## 0. Results

* Theorem A (Section 6) [Proved, with the computer-checked step of Theorem 4.4]. For every lambda > 0, L >= 2, 0 <= n <= L/2,

      |F_L(n) - L e_bulk(Delta, n/L)| <= K(lambda)/L,   K(lambda) <= 218.4 ab <= 4311 M_0,

  with K <= 218.4 ab and ab/M_0 = 4 pi^2 tanh(lambda/2)/lambda, so K/M_0 <= 218.4 * 4 pi^2/lambda -> 0 as Delta -> infinity.
* Corollary B (Section 8) [Proved]. For every bond of the U(1) class with J_z/J > 1, every field, DM angle and L >= 2,

      L(E_(2L) - 2E_L) <= 10781 M.

  The exact lambda-dependent constant evaluated on a grid gives 10386 (constants.json). Together with C142 (4078 on 0 <= J_z/J <= 1) the single constant 10781 holds on all of J_z/J >= 0.
* No region is left open on J_z/J > 1. The numerical margin is large: the largest computed L(F_(2L)(2n) - 2F_L(n))/M_0 is 3.9993 (ED, L = 2, Delta = 1.001), and L|F_L - L e_bulk|/M_0 <= 1.82 in 261 Bethe cases, against the proved 4311.

## 1. Bethe equations on the circle

[Proved, cited: C129 Sections 2-6] The C129 reductions put the target on the sector energies F_L(n) of K_L = sum_x [-(S^x S^x + S^y S^y) + Delta S^z S^z] (flux L pi), n down spins, 0 <= n <= L/2, with the field cancelling in matched sectors (2L, 2n) versus (L, n), the twist cost at most pi^2 M/(4L) paid on the doubled ring, and M >= M_0 = (1 + Delta)/4 for Delta >= 1. The sector Hamiltonian has off-diagonal entries -1/2 in the S^z basis, so it is stoquastic and irreducible for 0 < n < L, F_L(n) is a simple eigenvalue for every real Delta, and Delta -> F_L(n) is real-analytic on R.

Set (all functions of x in T, with p, theta continued so that p(x + pi) = p(x) + 2pi, theta(x + pi) = theta(x) + 2pi)

    p(x) = 2 atan[coth(lambda/2) tan x],   theta(x) = 2 atan[coth(lambda) tan x],
    a1 = p'/(2pi) = sinh(lambda)/[pi(cosh(lambda) - cos 2x)],   K := a2 = theta'/(2pi) = sinh(2lambda)/[pi(cosh(2lambda) - cos 2x)],
    e(x) = -a a1(x) = -sinh(lambda)^2/(cosh(lambda) - cos 2x).

The Bethe equations with centered quantum numbers are L p(x_j) - sum_k theta(x_j - x_k) = 2pi I_j, I_j = j - (n+1)/2, and the Bethe energy is E = Delta L/4 + sum_j e(x_j). [Check: one magnon has e(0) = -(1 + Delta), the k = 0 magnon of K_L.] With mu = (1/L) sum delta_(x_j), nu = n/L,

    sigma_L = a1 - K (*) mu,   z_L(x) = p(x)/(2pi) - (1/(2pi L)) sum_k theta(x - x_k),   z_L' = sigma_L.

[Proved] Three exact facts, for an ordered reflection-symmetric solution x_1 < ... < x_n in (-pi/2, pi/2):
(F1) z_L(x_j) = I_j/L, so each inner gap g_j = [x_j, x_(j+1)] has int_(g_j) sigma_L = 1/L.
(F2) z_L(pi/2) = (1 - nu)/2. Indeed p(pi/2) = pi, and theta(pi/2 - y) + theta(pi/2 + y) = 2pi (theta odd, theta(x + pi) = theta(x) + 2pi), so by symmetry sum_k theta(pi/2 - x_k) = n pi. Hence int_T sigma_L = 1 - nu.
(F3) The out-arc mass is int_(x_n)^(pi - x_n) sigma_L = 1 - 2nu + 1/L >= 1/L, and int_(x_n)^(pi/2) sigma_L = (1 - 2nu)/2 + 1/(2L).
(F4) No collision across +-pi/2: x_1 + pi = x_n would give equal left sides minus 2pi(L - n) against 2pi(I_1 - I_n) = -2pi(n - 1), i.e. n = (L + 1)/2 > L/2. Distinct coordinates follow as in C139 Section 2.1.

[Numerical, ed_checks.json] Bethe energies equal the ED sector minima at flux pi per bond in all 324 sectors 1 <= n <= L/2, L = 2..12, Delta in {1.001, 1.01, 1.1, 1.5, 2, 3, 5, 10, 20}, maximum error 9.2e-14 (once Theorem 4.5 is proved this is a check, not an input).

Continuum. For Q in (0, pi/2], in = [-Q, Q], out = T \ in (an arc of length ell = pi - 2Q), K_Q the operator with kernel 1_in(x) K(x-y) 1_in(y). Let r_Q = (I + K_Q)^(-1) a1, Z_Q = (I + K_Q)^(-1) 1 on in, m(Q) = int_in r_Q, E(Q) = Delta/4 + int_in e r_Q, h(Q) = -a r_Q(Q)/Z_Q(Q), eps_Q = (I + K_Q)^(-1)(e - h), so eps_Q(+-Q) = 0. I + K_Q is invertible on L^2(in) because Khat_k = t^|k| > 0 makes K_Q a positive operator. The C136/C139 identities (sign-free, they use only self-adjointness and eps_Q(+-Q) = 0) give m'(Q) = 2 r_Q(Q) Z_Q(Q) > 0 (Lemma 3.2), e_bulk(m(Q)) := E(Q), e_bulk' = h, e_bulk'' = eps_Q'(Q)/[2 r_Q(Q) Z_Q(Q)^2]. m maps [0, pi/2] onto [0, 1/2]; Q = pi/2 is half filling. At Q = pi/2 everything is a Fourier multiplier and e_bulk(1/2) = Delta/4 - sinh(lambda)[1/2 + 2 sum_(k>=1) 1/(e^(2k lambda) + 1)], the known value [checked numerically to 1e-14]. The lower critical field is h_c1 = -e_bulk'(1/2^-) = 2a r0(pi/2) = (2/pi) sinh(lambda) K(k) k' (r0 from Section 2) [Numerical: equals F_L(L/2 - 1) - F_L(L/2) at L = 1024 to 4e-6 at Delta = 2, 5; ed_checks.json].

## 2. Kernels on the circle

Lemma 2.1 (circle convolution) [Proved]. If f, g >= 0 are circle-unimodal, so is f (*) g.
Proof. By the layer-cake representation it suffices to take f = 1_A, g = 1_B, A = {d(y,0) <= alpha}, B = {d(y,0) <= beta}. Then F(c) = (f (*) g)(c) = |A cap (B + c)| and, for 0 < c < pi/2, F'(c) = 1_A(c + beta) - 1_A(c - beta). For 0 <= c, beta <= pi/2, d(c - beta, 0) = |c - beta| <= min(c + beta, pi - c - beta) = d(c + beta, 0) (the second comparison is c <= pi/2 if c >= beta and beta <= pi/2 otherwise). So 1_A(c + beta) <= 1_A(c - beta) and F' <= 0. Evenness is clear. End.

Lemma 2.2 (kernels) [Proved].
(i) ahat1_k = q^|k|, Khat_k = t^|k|, so int_T a1 = int_T K = 1 (q := ||K||_1 = 1, versus 1 - 2gamma/pi < 1 on the line). a1 and K are circle-unimodal (cos 2x decreases in d(x,0)).
(ii) r0 := (I + K(*))^(-1) a1 has r0hat_k = 1/(2cosh(k lambda)), int_T r0 = 1/2, and r0 = (K_e/pi^2) dn(2K_e x/pi, k), where k, K_e are the elliptic modulus and complete integral of nome q (DLMF 22.11.3: (2K_e/pi) dn = 1 + 4 sum q^n cos(2nx)/(1 + q^(2n))). dn decreases on [0, K_e], so r0 is circle-unimodal, r0(0) = K_e/pi^2, r0(pi/2) = k' K_e/pi^2. Also r0 = sum_m R(x + m pi) with R(x) = 1/(2 lambda cosh(bx)) (Poisson summation, since Rhat(w) = 1/(2cosh(w lambda/2)) at w = 2k). Since |R'| <= bR and |R''| <= b^2 R pointwise, |r0'| <= b r0 and |r0''| <= b^2 r0.
(iii) The hole kernel J = (I + K(*))^(-1) K(*) has Jhat_k = t^|k|/(1 + t^|k|) = ahat1_k r0hat_k, so J = a1 (*) r0 >= 0 is circle-unimodal (Lemma 2.1), int_T J = 1/2 =: q_J, |J'| <= bJ, |J''| <= b^2 J (differentiate r0). (I + K(*))^(-1) = I - J(*) on bounded functions.
(iv) sum_(m>=1) J^((*)m) = K, since Jhat/(1 - Jhat) = t^|k|. So sum_(m>=2) J^((*)m)(0) = K(0) - J(0). All J^((*)m) have nonnegative coefficients, so their sup is at 0.
(v) K(0) = coth(lambda)/pi, ||K'||_(L1(T)) = 2(K(0) - K(pi/2)) = 4/(pi sinh 2lambda), J(0) = (1/pi) sum_k t^|k|/(1 + t^|k|). Poisson summation of the line kernel with transform 1/(1 + e^(|w| lambda)) (positive images) gives pi J(0) >= ln2/lambda; dropping terms gives pi J(0) >= 1/2 + 2 sum_(j=1..3) t^j/(1 + t^j).
(vi) |K'| is unimodal on (0, pi/2): its derivative has the sign of c^2 + Cc - 2 (c = cos 2x, C = cosh 2lambda), which is negative at c = -1, positive at c = 1, so changes sign once. Hence ||K''||_1 = 4 M_1, M_1 = max|K'|. From the Fourier series K' = -(4/pi) sum k t^k sin 2kx, M_1 <= (4/pi) t/(1-t)^2, and M_1/||K'||_1 <= coth(lambda)/2.
(vii) ell(x) := r0'/r0 = -(2K_e/pi) k^2 sn cn/dn has |ell| <= b, and Var_[-Q,Q](ell) <= 4b: with s = sn^2, d/du(sn cn/dn) has the sign of 1 - 2s + k^2 s^2, which has the single root s = 1/(1 + k') in (0,1), so sn cn/dn is unimodal on [0, K_e]. Also |ell'| <= |r0''|/r0 + ell^2 <= 2b^2.

Lemma 2.3 (tails and reciprocals) [Proved]. For 0 <= A <= pi/2:
(i) int_A^(pi/2) r0 <= min(lambda/2, pi/2 - A) r0(A), so int_out r0 <= min(lambda, ell) r0(Q).
(ii) r0 <= c_p R on [-pi/2, pi/2], c_p = 2 + 4/(e^(pi^2/lambda) - 1), hence r0(Q)/r0(x) <= 2c_p e^(-b(Q - |x|)) for |x| <= Q.
(iii) int_(-A)^A dx/r0 <= 2 kappa_A/(b r0(A)), kappa_A = min(pi^2/(2lambda), c_p) <= 2.63.
(iv) T(Q) := int_Q^(pi/2) r0 = (pi/2 - am(u_Q))/(2pi), u_Q = 2K_e Q/pi.
Proof. (i) int over the arc of the periodization is at most int_(|y|>A) R = 2 int_A^inf R <= (pi/b) R(A) <= lambda r0(A), using C142 Lemma I.1 (int_A^inf sech(by) dy <= (pi/(2b)) sech(bA)) and R <= r0; the other bound is monotonicity. (ii) For A in [0, pi/2], cosh(bA)/cosh(b(A + m pi)) is <= 1 for m = 0, -1, <= 2e^(-b m pi) for m >= 1 and <= 2e^(-b(|m| - 1) pi) for m <= -2; sum. Then R(Q)/R(x) <= 2e^(-b(Q - |x|)). (iii) Either 2A/r0(A) <= pi/r0(A), or int dx/R = (4lambda/b) sinh(bA) <= 2/(b R(A)) <= 2c_p/(b r0(A)). The bound 2.63 holds since c_p <= 2 + 4/(e^2 - 1) for lambda <= pi^2/2 and pi^2/(2lambda) < 1 otherwise. (iv) int dn = am. End.

## 3. The continuum on the circle

Fix Q in (0, pi/2]. Let J_out f = J (*) (1_out f) restricted to out, and omega_Q = sup_x int_out J(x - y) dy = int_(d(u,0) < ell/2) J <= 1/2 (rearrangement, circle-unimodality), which bounds J_out on L^inf and on L^1.

Lemma 3.1 (duality) [Proved; C139 Lemma 3.2 verbatim on T]. s_Q := a1 - K (*) (1_in r_Q) on T satisfies s_Q = r_Q on in and s_Q = r0 + J (*) (1_out s_Q). Likewise Ztilde = 1/2 + J (*) (1_out Ztilde) (since (I + K(*))^(-1) 1 = 1/2) and epstilde = -a r0 - h/2 + J (*) (1_out epstilde), with epstilde(+-Q) = 0.

Lemma 3.2 (density) [Proved]. All Neumann series on out are positive, so on T: r0 <= s_Q; on out s_Q <= r0(Q)/(1 - omega_Q) <= 2 r0(Q); h_Q := int_out s_Q <= int_out r0/(1 - omega_Q) <= 2 min(lambda, ell) r0(Q); m(Q) = 1/2 - h_Q/2 (integrate the dual equation over T). On in, U1 := J (*) (1_out s_Q) <= 2r0(Q) int_out J(x - y) dy <= r0(Q) <= r0(x), so

    r0 <= r_Q <= c_r r0,  c_r = 2,     U1/r0 <= w(x) := min(1, 2c_p e^(-b(Q - |x|))),   int_in w <= (2/b) kappa_w,

kappa_w = min(pi^2/(2lambda), 1 + ln(2c_p)) <= 2.66. Differentiating, |r_Q'| <= b r_Q, |r_Q''| <= b^2 r_Q. Also 1/2 <= Ztilde <= 1.

Lemma 3.3 (Z and the restricted inverse) [Proved]. On in, Z' = J' (*) (1_out Ztilde), so Var_in(Z) <= int_in int_out |J'(x - y)| dy dx <= int_T |u||J'(u)| du = 1/2 - pi J(pi/2) <= 1/2, using |in cap (out + u)| <= d(u, 0). If (I + K_Q) f = g on in, F = g 1_in - K (*) (1_in f) equals f on in and F = (I - J(*))(g 1_in) + J (*) (1_out F), so ||(I + K_Q)^(-1)||_(inf -> inf) <= (1 + 1/2)/(1 - 1/2) = 3.

Lemma 3.4 (dressed energy) [Proved; C139 Lemma 3.5 on T]. epstilde' = -a r0' + J (*) (1_out epstilde') (boundary terms vanish), so |epstilde'| <= ab s_Q, and |eps_Q'| <= ab r_Q on in. For 0 < x, y < pi/2, d(x, y) <= d(x, -y), so J(x - y) - J(x + y) >= 0 and the odd reduction gives eps_Q' >= 0 on (0, Q]. Hence 0 <= e_bulk'' <= ab r/(2r (1/2)^2) = 2ab =: a2. For v = eps_Q'/r_Q on in: v = -a ell + w_v, w_v = (U2 + a ell U1)/r_Q, U2 = J (*) (1_out epstilde'), |U2| <= ab U1, |U2'| <= ab^2 U1, |U1'| <= b U1. So |w_v| <= 2ab w and |w_v'| <= (ab^2 + 2ab^2 + ab^2 + 2ab^2) U1/r <= 6ab^2 w, and

    |v| <= ab,   Var_in(v) <= 4ab + 6ab^2 (2/b) kappa_w = ab(4 + 12 kappa_w) =: V_v.

[Numerical] The restricted solver (circle.py, Nystrom) reproduces r0 at Q = pi/2 to 1e-12 and the closed-form e_bulk(1/2) to 1e-12 at Delta = 1.001, 1.1, 2, 20.

## 4. Positivity of the counting density and identification

Let x be the roots of an ordered, reflection-symmetric real solution with n >= 2, Q = x_n, in = [x_1, x_n], v = sigma_L on out. Hypothesis (P): sigma_L >= 0 on in.

Lemma 4.1 (pairing on the circle) [Proved, under (P)]. Write sigma_L = r0 + J (*) (sigma_L dx) - J (*) mu (from sigma + K(*)sigma = a1 + K(*)(sigma dx - mu), apply I - J(*)).
(a) For x in in, sigma_L(x) >= r0(x) + int_out J(x - y) v(y) dy - (1/L) J(0) - (1/L) 1_case2 J(x - x_1) [case 2 means |x| + Q > pi/2, stated for x > 0, mirror otherwise], and in case 2 J(x - y) >= J(x - x_1) for all y in out.
(b) For y in out: r0 + J_out v - psi/L <= v <= r0 + J_out v, psi(y) = J(y - x_1) + J(y - x_n).
Proof. J(x - .) is nonincreasing in d(x, .). Case 1 (|x| + Q <= pi/2): in lies within distance pi/2 of x on the real line, so circular distance is |x - y| there and C140 Lemma 3.1(a) applies verbatim: gaps left of x pair with their left endpoint, gaps right of x with their right endpoint, g_k containing x with both endpoints; each pair is >= 0 by (P) and (F1), and g_k loses at most max(J(x - x_k), J(x - x_(k+1)))/L <= J(0)/L. Case 2 (x > 0): the antipode x* = x - pi/2 lies in in. On [-Q, x*], J(x - .) decreases in y (distance pi - (x - y) grows), on [x*, x] it increases. Pair the gaps in [-Q, x*] with their right endpoints and those in [x*, x] with their left endpoints, the gap containing x* with nothing (contributes >= 0), g_k with both endpoints and the gaps right of x with right endpoints. The root x_1 is left over. For y in out = (Q, pi - Q), d(x, y) = y - x <= pi - Q - x = d(x, x_1), since x + Q >= pi/2. (b) For y in out the antipode y* lies in out or in a gap g_m. Along each of the two arcs from y to y*, J(y - .) is monotone. For the lower bound pair every gap with its farther endpoint, leaving the root(s) nearest to y (one or both of x_1, x_n) and possibly the gap g_m (>= 0); for the upper bound pair every gap with its nearer endpoint (g_m with the endpoint of larger J), leaving one root with a minus sign. End.

Corollary 4.2 (outer Neumann bounds) [Proved, under (P)]. With phit = (I - J_out)^(-1) psi >= 0: s_Q - phit/L <= v <= s_Q on out. phit is even, ||phit||_1 <= ||psi||_1/(1 - omega_Q) <= 1/(1/2) = 2 (each J(. - x_i) has mass at most 1/2 on out), so the negative mass of sigma_L on each half (x_n, pi/2], [-pi/2, -x_n) is at most 1/L. Also J (*) (1_out phit) <= sum_(k>=0) [J^((*)(k+2))(. - x_n) + J^((*)(k+2))(. - x_1)] <= 2(K(0) - J(0)) by positivity and Lemma 2.2(iv).

Lemma 4.3 (location of the outermost root) [Proved, under (P)]. For n >= 2,

    L r0(x_n) >= Loc := max( 4 r0(Q), (1 - omega_Q) r0(Q)/(2 T(Q)) ) >= 4 sin(pi/8) k K_e/pi^2.

Proof. L >= 2n >= 4 gives the first entry. By (F3) and v <= s_Q, 1/L <= 1 - 2nu + 1/L = int_out v <= h_Q <= int_out r0/(1 - omega_Q) = 2T(Q)/(1 - omega_Q), the second entry. With phi = am(u_Q), r0(Q) = (K_e/pi^2) dn, and Lemma 2.3(iv), Loc >= (K_e/pi^2) dn(phi) max(4, pi/(pi - 2phi)) (using omega_Q <= 1/2). For phi <= 3pi/8 this is >= 4 dn(3pi/8) = 4 sqrt(1 - k^2 sin^2(3pi/8)) >= 4 sqrt(1 - sin^2(3pi/8)) = 4 sin(pi/8). For phi > 3pi/8 put psi_a = pi/2 - phi <= pi/8: dn^2 = k'^2 + k^2 sin^2 psi_a >= k^2 sin^2 psi_a, and dn pi/(2psi_a) >= k (pi/2) sin(pi/8)/(pi/8) = 4k sin(pi/8). End.

Theorem 4.4 (positivity with margin) [Proved, computer-checked on 1 <= lambda <= 1.36]. Along any continuous family of ordered symmetric solutions on which (P) holds, in fact

    sigma_L >= (1 - Phi(lambda)) r0(x_n) on in,   Phi(lambda) <= 0.939   (0 < lambda <= 1.36),
    sigma_L >= [tanh(lambda/2) - coth(lambda)/2]/pi >= 0.021/pi > 0 on T, for every configuration   (lambda >= 1.36).

Proof. For lambda >= 1.36 and any n <= L/2 points: sigma_L >= min a1 - nu max K >= a1(pi/2) - K(0)/2 = [tanh(lambda/2) - coth(lambda)/2]/pi, which is increasing in lambda and equals 0.0210/pi at 1.36. For lambda <= 1.36 and n >= 2: in case 1 of Lemma 4.1(a) insert v >= s_Q - phit/L and drop U1 >= 0; in case 2 use U1 >= J(x - x_1) h_Q >= J(x - x_1)/L (Lemma 4.1(a), h_Q >= 1/L from the proof of Lemma 4.3) to cancel the left-over root. Either way, with r0(x) >= r0(x_n) on in and Corollary 4.2, sigma_L(x) >= r0(x_n) - Lambda/L, Lambda = J(0) + 2(K(0) - J(0)) = 2K(0) - J(0). So sigma_L >= (1 - Phi) r0(x_n) with Phi = Lambda/(L r0(x_n)) <= (2K(0) - J(0)) pi^2/(4 sin(pi/8) k K_e) by Lemma 4.3. Now k K_e = (pi/2) theta_2(q)^2 and, by Poisson summation with alternating decreasing terms, theta_2(e^(-lambda)) >= sqrt(pi/lambda)(1 - 2e^(-pi^2/lambda)). So

    Phi(lambda) <= lambda [2coth(lambda) - pi J(0)] / [2pi sin(pi/8) (1 - 2e^(-pi^2/lambda))^2].

On (0, 1] use pi J(0) >= ln2/lambda: the bound is increasing in lambda, equal to 0.804 at lambda = 1 and to (2 - ln2)/(2pi sin(pi/8)) = 0.544 as lambda -> 0. On [1, 1.36], on each of 36 cells [l0, l1] bound lambda by l1, coth(lambda) by coth(l0), pi J(0) from below by 1/2 + 2 sum_(j<=3) t^j/(1+t^j) at l1 (decreasing in lambda) and the denominator from below at l1; the largest cell value is 0.939 (pos_rigorous.py, pos_rigorous.json). n = 1: x_1 = 0 and a1/K >= Delta/(1 + Delta) > 1/2 >= 1/L, so sigma_L > 0 on T. End.

[Numerical, pos_check.json] The sharper exact inequality (with (K - J)(min(Q, pi - 2Q)) in place of K(0) - J(0) and the exact omega_Q) has sup ratio 0.53 to 0.74 on 0 < lambda <= 1.5. The actual min sigma_L/r0 over in is >= 0.974 and over all of T >= 0.969 in all 261 computed cases (scan_small.jsonl), and the actual L r0(x_n) >= 0.64.

Theorem 4.5 (branch and identification) [Proved]. For every L >= 2 and 1 <= n <= L/2 there is a real-analytic family Delta -> x(Delta), Delta in (1 - eps, infinity), of ordered reflection-symmetric real solutions with the ground-state I_j, with positive definite Yang-Yang Hessian and (P) strict, whose Bethe energy equals F_L(n). Moreover the negative mass of sigma_L on each half of out is at most 1/L, and is zero for lambda >= 1.36.
Proof. (1) Through Delta = 1. Use XXX rapidities u: x = gamma u for Delta = cos gamma < 1, x = lambda u for Delta = cosh lambda > 1. With s = gamma^2 (s = -lambda^2 for Delta > 1), tanh(gamma u)/tan(gamma/2) and tanh(gamma u)/tan(gamma) are even in gamma and become tan(lambda u)/tanh(lambda/2), tan(lambda u)/tanh(lambda) at gamma = i lambda, so the Bethe equations and the Bethe energy are real-analytic in (s, u) near s = 0 for bounded u, and s -> Delta = cos(sqrt s) is an analytic bijection near 0. C142 Theorem I.4 gives on 0 < gamma <= pi/2 a branch with sigma >= theta_1 r0 on [x_1, x_n], theta_1 >= 0.168, and Lemma I.2 there gives cosh(pi u_n) <= L, so the u-roots stay bounded. Take gamma_k -> 0 with u(gamma_k) -> u*. u* is an ordered XXX critical point with distinct coordinates, and in u-variables sigma >= theta_1 r0 passes to the limit, so the XXX Hessian 2 pi L diag(sigma(u_j)) + Theta (Theta positive definite by Bochner) is positive definite. The implicit function theorem at (s, u) = (0, u*) and the open-closed argument of C142 Step 3 identify C142's branch with the IFT branch near s = 0 and continue it to small s < 0. (2) Along lambda in (0, infinity) in the circle variables x = lambda u. The set of lambda_1 such that the branch extends analytically to (0, lambda_1] with (P) strict and positive definite Hessian is open (IFT; strict (P) at the roots makes Y'' = 2 pi L diag(sigma(x_j)) + Theta positive definite, Theta_jk = theta'(x_j - x_k), positive definite since Khat_k > 0). It is closed: T is compact, limit points are ordered critical points with distinct coordinates and no collision across +-pi/2 (F4), (P) holds non-strictly in the limit, and Theorem 4.4 upgrades it to strict with margin, so the Hessian is positive definite at the limit and the C142 Step 3 argument continues the branch. Hence the branch exists on (0, infinity). Reflection symmetry follows from local uniqueness. (3) Identification. E(Delta) = Delta L/4 + sum e(x_j(Delta)) is real-analytic on (1 - eps, infinity) and equals F_L(n) on (1 - eps, 1) (C142 Theorem I.4); F_L(n) is real-analytic on R (Section 1). By the identity theorem E = F_L(n) on (1 - eps, infinity). The negative-mass statement is Corollary 4.2 and the first line of Theorem 4.4. End.

## 5. Sawtooth, mismatch, cumulative error, density

Fix lambda > 0, L >= 2, 1 <= n <= L/2 and the roots of Theorem 4.5. Q = Q_L is the first x >= x_n with z_L(x) = nu/2; it exists in (x_n, pi/2] since z_L(x_n) = nu/2 - 1/(2L) and z_L(pi/2) = (1 - nu)/2 >= nu/2. r = r_Q, Z, eps, m, h as in Section 1, R = mu - 1_in sigma_L dx, Cum R(x) = R([-Q, x]) = C_mu(x) - z_L(x) on in.

Lemma 5.1 (sawtooth) [Proved]. sup_in |Cum R| <= s/L with s = 3/2 for lambda < 1.36 and s = 1/2 for lambda >= 1.36. Proof. C142 Lemma I.5 verbatim: on [x_1, x_n] z_L is increasing (Theorem 4.4) and Cum R runs from 1/(2L) to -1/(2L) on each gap; on (x_n, Q], 0 <= Cum R = nu/2 - z_L <= 1/(2L) + (negative mass on (x_n, pi/2]) <= 1/(2L) + 1/L, and zero extra for lambda >= 1.36. End.

Lemma 5.2 (mismatch and G) [Proved]. (I + K_Q) nu' = R on in as measures (nu' = mu - r dx), |nu - m| <= Var(Z) sup|Cum R| <= s/(2L), and G(x) = nu'([-Q, x]) satisfies ||G||_inf <= g/L, g = 13.5 s. Proof. C142 Lemma I.6 on T with k(y) = int_0^y K: integrating over [-Q, x] gives G + K_Q G = Cum R + c_0 - (nu - m) psi_G, psi_G(x) = k(x - Q) + k(2Q) in [0, k(2Q)] subset [0, 1]. So G = g_1 + c_0 Z with ||g_1|| <= 3(s + s/2)/L (Lemma 3.3), |c_0| Z(x) <= ||g_1|| Z(x)/Z(-Q) <= 2||g_1||, ||G|| <= 3||g_1|| <= 13.5 s/L. End.

Lemma 5.3 (density) [Proved]. On in, sigma_L - r = -K (*) nu' = -K(x - Q)(nu - m) - int_in K'(x - u) G(u) du, so

    |sigma_L - r| <= D/L,  D = g||K'||_1 + K(0) s/2,     |sigma_L' - r'| <= D1/L,  D1 = M_1(4g + s/2).

[Numerical, scan_small.jsonl] L sup|Cum R| = 0.500, L||G|| <= 0.568, L|nu - m| <= 0.25, lambda L sup|sigma - r| <= 0.079, against s = 3/2, g = 20.25.

## 6. Quadrature and the sector theorem

Exact identity [Proved; C139 (5.1), sign-free, rederived for T]. F_L(n) - L[E(Q) + h(nu - m)] = L int eps dR = -L int_in eps' Cum R dx. (Sum e over roots, write e - h = (I + K_Q) eps on in, use self-adjointness, K_Q mu = (a1 - sigma) 1_in dx, and E(Q) - Delta/4 = int_in eps a1 + hm.)

Core [Proved]. If r0(0) >= 2D/L let A in [0, pi/2] be the largest point with r0(A) >= 2D/L (A = pi/2 if r0(pi/2) >= 2D/L). On [-A, A] cap in, sigma_L >= r - D/L >= r/2 > 0. The composite trapezoid rule for f = eps o z_L^(-1) with step 1/L gives an error at most Var(f')/(8L), f' = v g_r, g_r = r/sigma_L in (0, 2]. |g_r'| <= 4(bD + D1)/(L r), and by Lemma 2.3(iii) int_(-A)^A dx/r0 <= kappa_A L/(bD), so Var(g_r) <= 4 kappa_A (1 + D1/(bD)), and

    V_f = ab Var(g_r) + 2 V_v = ab[4 kappa_A(1 + D1/(bD)) + 8 + 24 kappa_w].

Ends [Proved]. As in C142 Lemma I.8, each end is at most s ab int_(x_j2)^Q r, int_(x_j2)^A r = G(x_j2^+) - G(A) <= 2g/L, and int_A^Q r <= c_r int_A^(pi/2) r0 <= 2 min(pi/(2b), pi/2) r0(A) = 4 min(pi/(2b), pi/2) D/L (Lemma 2.3(i)). Both ends: 2s ab[2g + 4 min(pi/(2b), pi/2) D]/L. Empty core (r0(0) < 2D/L): |L int eps dR| <= s ab int_in r <= s ab c_r 2 min(pi/(2b), pi/2) r0(0) <= the ends bound, so no threshold on L is needed. Curvature: |e_bulk(nu) - e_bulk(m) - h(nu - m)| <= a2 (nu - m)^2/2 <= ab s^2/(4L^2).

Theorem A [Proved]. For every lambda > 0, L >= 2, 0 <= n <= L/2,

    |F_L(n) - L e_bulk(Delta, n/L)| <= K(lambda)/L,
    K(lambda) = V_f/8 + 2s ab[2g + 4 min(pi/(2b), pi/2) D] + ab s^2/4.

n = 0 is exact. The case Q_L = pi/2 (even half filling with no earlier crossing) is included: then out is empty, h = -2a r0(pi/2), eps(+-pi/2) = 0.

## 7. Constants

[Proved] Analytic majorant, with s = 3/2, g = 20.25 at every lambda:
* D/b-type factor: min(pi/(2b), pi/2) D <= g (lambda/2)||K'||_1 + (s/2) K(0) min(lambda/2, pi/2) <= g/pi + (s/2)(pi/2)coth(pi)/pi <= 6.446 + 0.376 = 6.822 (2lambda <= sinh 2lambda; lambda coth(lambda)/2 increasing, and (pi/2)coth(lambda) <= (pi/2)coth(pi) for lambda >= pi).
* Ends <= 2(3/2) ab [40.5 + 4(6.822)] = 203.4 ab.
* kappa_A D1/(bD) <= (4 + s/(2g)) kappa_A M_1/(b||K'||_1) <= (4 + 1/27) kappa_A lambda coth(lambda)/(2pi) <= 4.037 * 2.07 = 8.35 (kappa_A lambda coth(lambda)/(2pi) <= 2.626 * 4.935 coth(4.935)/(2pi) = 2.07 for lambda <= pi^2/2, and <= (pi/4)coth(lambda) < 0.79 beyond).
* V_f <= ab[4(2.63) + 4(8.35) + 8 + 24(2.66)] = 115.8 ab, V_f/8 <= 14.5 ab.
* Curvature ab s^2/4 = 0.5625 ab.
* K <= 218.4 ab, and ab/M_0 = 4pi^2 tanh(lambda/2)/lambda <= 2pi^2, so K <= 4311 M_0.

[Numerical evaluation of the exact formula, constants.py, constants.json] With the exact lambda-dependent s, g, D, D1, kappa_A, kappa_w on 670 points of lambda in [1e-4, 200], sup K/M_0 = 4153.4, approached as lambda -> 0 (ends 3953, core 189, curvature 11). For lambda >= 1.36 (s = 1/2) K/M_0 <= 462 and decays like 1/lambda.

## 8. Assembly and the glue at Delta = 1

Corollary B [Proved]. Matched doubling: L[F_(2L)(2n) - 2F_L(n)] <= L[K/(2L) + 2K/L] = (5/2)K <= 10778 M_0 (the extensive terms cancel since (2L, 2n) and (L, n) have the same density). With the C129 reductions (field cancellation in matched sectors, twist cost pi^2 M/(4L) on the doubled ring, M >= M_0, every DM angle, every L >= 2, sign-free and refereed in C136):

    L(E_(2L) - 2E_L) <= (10778 + pi^2/4) M <= 10781 M    for all J_z/J > 1.

(Exact formula: 2.5 * 4153.4 + pi^2/4 = 10386.)

The other regions: J_z/J = 1 is covered by C142 Corollary I.10 (4078) and also by continuity of the present bound at fixed L. The saturation region |B| >= J + J_z has zero defect (C129 Section 9); it is consistent with, not needed by, Corollary B. Glue at Delta = 1: both C142 and Theorem A are uniform up to Delta = 1 from their sides, with the same C129 reductions and the same M >= M_0 (M_0 = J/2 at Delta = 1 on both sides), and the identification is one analytic branch through Delta = 1 (Theorem 4.5). So on J_z/J >= 0 the single constant max(4078, 10781) = 10781 holds; combined with C142's assembly, on J_z/J >= -1 + delta the constant is max(10781, 6e15 delta^(-25/3)).

The four regions named in the spec needed no separate theorem: (a) Delta -> 1+: all constants have finite limits as lambda -> 0, where they reduce to the gamma -> 0 constants of C142 except s (3/2 here against 1 there, from the two-sided out arc) and c_r (2 against 1 + 4/pi). (b) Near h_c1 (Q near pi/2, dilute holes): Lemmas 3.2 to 6 are uniform in Q including Q = pi/2; the dilute-hole region is the case of a short out arc, where omega_Q and h_Q are small and Lemma 4.3 is strongest. C129's dilute bound is not needed. (c) Half filling: for even L, Q_L may equal pi/2 and Theorem A holds as stated; for odd L the sector n = (L - 1)/2 has one vacancy (a spinon at x = pi/2). The defect there is positive (Section 9). (d) Large Delta: K/M_0 <= 218.4 * 4 pi^2/lambda -> 0.

## 9. Consistency checks

* Ising limit, lambda -> infinity [Proved]: p(x) -> 2x, theta -> 2x, the Bethe equations become 2(L - n) x_j + 2 sum x_k = 2pi I_j, so x_j = pi I_j/(L - n), and e(x) = -cosh(lambda) - cos 2x + O(e^(-lambda)): free hard-core fermions on L - n sites, as expected from the exclusion picture. Theorem 4.4's direct bound becomes sigma_L >= (1 - nu)/pi - o(1). K/M_0 -> 0. [Numerical] The defect/M_0 is 0.17 at Delta = 20 (L = 128) and 0.32 at Delta = 10, decreasing like J/M_0.
* Heisenberg, lambda -> 0 [Proved]: in u-variables a1, K, r0, J converge to the XXX kernels; Phi -> 0.544, and kK_e/pi^2 ~ 1/(2lambda) reproduces C142's location scale b/(2pi). [Numerical] The ED end-to-end ratio at Delta = 1.001, L = 2 is 3.9993, continuous with the Heisenberg supremum 4; the half-filled Bethe defect at Delta = 1.001, L = 128 is 2.4711, against 2.4689 at Delta = 1, L = 1024 (C139).
* Saturation and n = 0: exact, zero defect.
* Near h_c1 [Numerical]: h_c1 = 2a r0(pi/2) equals F_L(L/2 - 1) - F_L(L/2) at L = 1024 to 4e-6 (Delta = 2, 5) and converges at Delta = 1.5 (0.0866 against 0.0873). At Delta = 1.01, 1.1, h_c1 = 9e-15, 2e-4 is far below 1/L, so finite rings there behave like XXX (finite-size gaps 0.0046, 0.0056 at L = 1024).
* Half filling parity (task (c)) [Numerical, ed_checks.json]: even L, n = L/2: the defect is exponentially small in L once L exceeds the correlation length (7e-5 at Delta = 2, L = 130; 1e-11 at Delta = 5, of either sign at roundoff). Odd L, n = (L-1)/2: the defect is positive, L(F_(2L)(L - 1) - 2F_L((L-1)/2))/M_0 = 1.82, 1.38, 0.51 at Delta = 2 and L = 9, 33, 129, and 0.80, 0.29, 0.08 at Delta = 5. [Heuristic] The L-ring puts its single spinon at the band minimum x = pi/2; the doubled ring holds two spinons with distinct quantum numbers, which cannot both sit at the minimum, so it pays a positive kinetic cost of order 1/(m* L^2) plus a scattering shift. So the odd-ring spinon contributes to the defect with a positive sign, but its contribution vanishes like a power of 1/L, and in any case Theorem A bounds it.
* Bethe vs ED: 324 sectors, error 9.2e-14. End-to-end ED with field (43 values in [0, 1.05 h_c2], including fields just above h_c1) and DM in {0, 0.3, 1.0}, L = 2..6, nine Delta: 5805 samples, max L(E_2L - 2E_L)/M = 3.9993, min -17.6 (negative defects occur at level crossings in field). Bethe doubling to 2L = 1024: 621 cases (scan_small.jsonl, scan_large.jsonl), max 2.79 (L = 4), 2.47 at L >= 128, no solver failures.

## 10. What remains

Nothing on J_z/J > 1 is open at the level of Corollary B. The one computer-checked step is the cell bracketing of Phi(lambda) on [1, 1.36] (36 cells of explicit monotone elementary functions, floating point, largest value 0.939). It could be replaced by interval arithmetic without change of method. The constant 10781 is crude; the losses are s = 3/2 (two-sided out arc in Corollary 4.2) and the first-order ends. No independent referee pass has been made.

=====================================================================
DOCUMENT B (excerpt, Part I of the earlier proof, real line, Delta = cos gamma)
=====================================================================
# C142 proof: the two U(1) strips

2026-11-07. Run by a Claude Opus agent, not by Codex. J=1 in all analytic sections (restore J by homogeneity). Tags: [Proved] is a mathematical deduction, [Proved, cited] uses a result of an earlier run that has had a referee pass (named), [Numerical] is floating point without enclosures, [Heuristic] has an unproved step, [Assumed] is a hypothesis.

Part I closes the repulsive strip 63/64 < Delta <= 1 and redoes C139 Section 8.2 line by line. Part II proves a version of (D) and (D1) near Delta = -1 that is weaker by a factor log(1/eta), and isolates what remains. Part III is the assembly.

# Part I. The repulsive side, 0 <= Delta <= 1

Notation of C139 Section 0 and C140: Delta = cos gamma, 0 < gamma <= pi/2, b = pi/gamma, a = pi sin gamma, ab = pi^2 sin(gamma)/gamma <= pi^2, q = 1 - 2gamma/pi, q_J = (pi-2gamma)/(2(pi-gamma)) < 1/2, r0(x) = 1/(2 gamma cosh(bx)), K the positive kernel, J = (I+K*)^(-1)K* >= 0 the hole kernel (even, unimodal, ||J||_1 = q_J), s_Q = r0 + J*(1_out s_Q), out = {|x| > Q}, h_Q = int_out s_Q. Ground-state roots x_1 < ... < x_n, mu = (1/L) sum delta_(x_j), sigma = sigma_L = a1 - K*mu, counting function z_L with z_L' = sigma, z_L(x_j) = I_j/L, I_j = j - (n+1)/2, nu = n/L, 1 <= n <= L/2.

## I.1 Where C140 lost the factor

[Proved] Lemma I.1 (sharp tail integral). For every Q >= 0,

    int_Q^inf r0(x) dx <= (pi/(2b)) r0(Q),

with equality at Q = 0, and the ratio of the two sides tends to 2/pi as Q -> infinity.

Proof. With t = e^(-bQ) in (0,1], int_Q^inf sech(bx) dx = (2/b) atan(t), so int_Q^inf r0 = atan(t)/(gamma b), while r0(Q) = t/(gamma(1+t^2)). The ratio is f(t)/b with f(t) = atan(t)(1+t^2)/t. Then f'(t) = [t - (1-t^2) atan(t)]/t^2 >= 0 on (0,1], since atan(t) <= t and 1-t^2 <= 1. So f is increasing, f(1) = pi/2, f(0+) = 1. End.

C140 Lemma 3.3 used atan(t) <= t <= 2t/(1+t^2), i.e. the constant 2 in place of f(t) <= pi/2. This factor 4/pi is the whole 6 percent shortfall.

[Proved, under (P): sigma >= 0 on [x_1, x_n]] Lemma I.2 (location of the outermost root).

    L r0(x_n) >= b(1-q_J)/pi = b/(2(pi-gamma)) >= b/(2pi) = 0.15915 b.

Proof. C140 Lemma 3.3 verbatim up to its last line (refereed in C140 REFEREE Section 1): by (F3), Corollary 3.2 and C139 Lemma 3.3, 1/(2L) <= int_(x_n)^inf sigma <= int_Q^inf s_Q = h_Q/2 <= int_Q^inf r0/(1-q_J), with Q = x_n. Lemma I.1 gives 1/(2L) <= (pi/(2b)) r0(Q)/(1-q_J). Finally 1-q_J = pi/(2(pi-gamma)). End.

The needed value is lambda_* b = ((2-ln2)/pi^2) b = 0.13241 b. Since 2pi(2-ln2) = 8.211 < pi^2 = 9.870, the bound exceeds it at every gamma, by a factor at least pi^2/(2pi(2-ln2)) = 1.2020 (at gamma -> 0). So (L*) holds, strictly, along the branch.

[Proved, under (P)] Theorem I.3 (uniform margin). For x in [x_1, x_n],

    sigma(x) >= theta_1(gamma) r0(x),   theta_1(gamma) = 1 - 2(2-ln2)(pi-gamma)/pi^2 >= 1 - 2(2-ln2)/pi = 0.16803.

Proof. C140 Theorem 3.4 gives sigma(x) >= r0(x) - (2-ln2)/(pi gamma L). Since r0(x) >= r0(x_n) on [x_1, x_n] (symmetry and monotonicity of r0), Lemma I.2 gives (2-ln2)/(pi gamma L) <= (2-ln2) r0(x) 2(pi-gamma)/(pi gamma b) = r0(x) 2(2-ln2)(pi-gamma)/pi^2, using pi gamma b = pi^2. End.

[Numerical, check_repulsive.json, 160 cases, Delta in {0.5, 0.9, 0.99, 0.999}, L from 4 to 256, seven fillings per L] L r0(x_n)/b >= 0.442 against the proved 0.159 to 0.239, and min sigma/r0 on [x_1, x_n] >= 0.974 against theta_1 between 0.168 and 0.45.

## I.2 The continuation argument on all of 0 < gamma <= pi/2

[Proved] Theorem I.4. For every 0 < gamma <= pi/2 (0 <= Delta < 1), every L >= 2 and 1 <= n <= L/2, the branch of real, strictly ordered critical points of the Yang-Yang functional continued from the free solution at gamma = pi/2 exists on [gamma, pi/2], is real-analytic and reflection symmetric, its Hessian is positive definite, and

    (i)  sigma_L >= theta_1 r0 on [x_1, x_n],
    (ii) int_(x > x_n) (sigma_L)_- dx <= (pi-2gamma)/(2pi L), and the same on x < -x_n.

Its Bethe energy is the sector ground energy F_L(n) at flux pi per bond, at both parities of L.

Proof. n = 1 is C140 Theorem 4.1 (x_1 = 0, sigma > 0 on R). Let n >= 2 and fix gamma_0 in (0, pi/2]. At gamma = pi/2, K = 0, sigma = a1 > 0, and the Hessian Y'' = 2pi L diag(sigma(x_j)) + Theta (C139 Lemma 2.2, Theta positive semidefinite by Bochner) is positive definite. By the implicit function theorem there is a unique real-analytic branch x(gamma) of critical points near pi/2. Let gamma_* be the infimum of the gamma_1 in [gamma_0, pi/2] such that the branch extends real-analytically, through critical points with positive definite Hessian, to (gamma_1, pi/2]. Suppose gamma_* >= gamma_0.

Step 1 ((P) along the branch). On the connected set (gamma_*, pi/2] the branch is continuous, strictly ordered, and reflection symmetric (the reflected family solves the same equations, agrees with x(gamma) at pi/2, and the coincidence set is open by local uniqueness and closed by continuity). The set where (P) holds is closed in (gamma_*, pi/2], since sigma_L(x; gamma) is jointly continuous and [x_1, x_n] moves continuously. Where (P) holds, Theorem I.3 gives sigma >= theta_1 r0 > 0 on the compact [x_1, x_n], and by continuity (P) holds strictly nearby, so the set is open. It contains pi/2, hence it is all of (gamma_*, pi/2], and (i) holds there.

Step 2 (limit point). The roots stay in the a priori ball of C139 Section 2.1, whose radius is bounded for gamma in [gamma_0, pi/2]. Take gamma_k decreasing to gamma_* with x(gamma_k) -> x_*. Then x_* is a critical point at gamma_* (grad Y is jointly continuous), with distinct coordinates (C139 Section 2.1), ordered and symmetric. By (i) and Lemma I.2, sigma(x_j(gamma_k)) >= theta_1 r0(x_n(gamma_k)) >= theta_1 b/(2(pi-gamma_k)L), which passes to the limit, so the Hessian at x_* is positive definite.

Step 3 (no change of branch). The implicit function theorem at (x_*, gamma_*) gives a neighbourhood U of x_*, an eps > 0, and a real-analytic x~ on (gamma_* - eps, gamma_* + eps) such that x~(gamma) is the only critical point in U at each such gamma. The set I = {gamma in (gamma_*, gamma_* + eps) : x(gamma) = x~(gamma)} contains gamma_k for k large (x(gamma_k) is in U), is closed (continuity), and is open (if x(gamma) = x~(gamma) is in U, then x(gamma') is in U for gamma' near gamma, hence equals x~(gamma')). So I = (gamma_*, gamma_* + eps), and x~ continues the branch analytically below gamma_*, with positive definite Hessian near gamma_*. This contradicts the definition of gamma_* unless gamma_* < gamma_0. Hence the branch exists on [gamma_0, pi/2], with (i).

(ii) From C140 Corollary 3.2, (sigma)_- <= phit/L on out, phit = (I-J_out)^(-1) psi. psi and J_out commute with the reflection, so phit is even and int_(y > x_n) phit = ||phit||_1/2 <= ||psi||_1/(2(1-q_J)) = q_J/(2(1-q_J)) = (pi-2gamma)/(2pi). (C140 stated the two-sided total 1/L.)

Identification. The facts used, all sign-free (C139 Section 1 item 2, C134 Section 2, C140 REFEREE Section 2): the sector Hamiltonian at flux pi per bond (total L pi, at either parity) has off-diagonal entries -1/2 in the S^z basis, so it is stoquastic and irreducible and F_L(n)(gamma) is a simple, real-analytic eigenvalue; along the analytic branch the coordinate Bethe vector is an analytic eigenvector with eigenvalue E(gamma) = Delta L/4 + sum e(x_j) wherever it is nonzero; at gamma = pi/2 it is the free-fermion Slater determinant of the symmetric Fermi sea, nonzero and the sector ground state. So E = F_L(n) near pi/2, and E and F_L(n) are both real-analytic on [gamma_0, pi/2], hence equal there by the identity theorem. End.

[Numerical] The integrated negative mass beyond x_n is zero to grid accuracy in all 160 cases (sigma_L > 0 on R), consistent with C140's sigma_L >= 0.952 r0 on R.

## I.3 Theorem C of C139, redone line by line with the proved sawtooth

Fix 0 < gamma <= pi/2, L >= 2, 1 <= n <= L/2, and the roots of Theorem I.4. Let Q = Q_L be the first x >= x_n with z_L(x) = nu/2 (Q_L > x_n because z_L(x_n) = nu/2 - 1/(2L); Q_L is finite for nu < 1/2 because z_L(inf) = c > nu/2; at even half filling Q_L = infinity is allowed). r = r_Q, Z = Z_Q, eps = eps_Q, m = m(Q), h = h(Q) are the restricted continuum objects of C139 Section 3, R = mu - 1_[-Q,Q] sigma dx, Cum R(x) = R((-inf, x]) = C_mu(x) - z_L(x).

[Proved, cited] Continuum inputs, all uniform in Q and gamma (C139 Lemmas 3.3 to 3.5, refereed in C140 REFEREE Sections 3 and 4): r0 <= r <= c_r r0 with c_r = 1 + 4/pi; |r'| <= b r; 1/(1+q) <= Z <= 1; Var_[-Q,Q](Z) <= q_J; ||(I+K_Q)^(-1)||_(inf->inf) <= (1+q_J)/(1-q_J) <= 3; |eps'| <= ab r; v = eps'/r has |v| <= ab and Var_[-Q,Q](v) <= V_v = (2+40/pi) ab; 0 <= e_bulk'' <= a2 = ab(1+q)^2/2 <= 2pi^2; eps(+-Q) = 0. At Q = infinity the same hold with r = r0, Z = 1/(1+q), v = ab tanh(bx).

[Proved, cited] Exact identity (C139 (5.1), sign-free, rederived in C140 REFEREE): with nu' = mu - 1_Q r dx, (I + K_Q) nu' = R on [-Q,Q] as measures, and

    F_L(n) - L[E(Q) + h(nu - m)] = L int eps dR = -L int_(-Q)^Q eps'(x) Cum R(x) dx.        (I.1)

The second form is integration by parts using eps(+-Q) = 0 and Cum R(+-Q) = 0.

### Lemma I.5 (sawtooth) [Proved]

    sup_[-Q,Q] |Cum R| <= s/L,   s = 1 - gamma/pi in [1/2, 1).

Proof. On (x_j, x_(j+1)), C_mu = j/L - nu/2 and z_L increases (Theorem I.4(i)) from I_j/L to I_j/L + 1/L, so Cum R decreases from j/L - nu/2 - I_j/L = 1/(2L) to -1/(2L). On (x_n, Q], Cum R = nu/2 - z_L(x) >= 0 by the definition of Q, and z_L(x) >= z_L(x_n) - int_(x_n)^x (sigma)_- >= nu/2 - 1/(2L) - (pi-2gamma)/(2pi L) by Theorem I.4(ii), so Cum R <= (1/2 + (pi-2gamma)/(2pi))/L = s/L. The left side is the mirror image. At Q = infinity the same argument applies on (x_n, inf). End.

[Numerical] L sup|Cum R| is between 0.4997 and 0.5000 in all 160 cases (no negative mass occurs numerically), against s in [0.67, 0.9997].

### Lemma I.6 (mismatch and cumulative error) [Proved]

    |nu - m| <= q_J s/L,
    ||G||_inf <= g/L,   g = (2+q) [(1+q_J)/(1-q_J)] s (1 + q q_J/2) <= 45/4,

where G(x) = nu'((-inf, x]) = C_mu(x) + nu/2 - int_(-Q)^x r.

Proof. Mismatch: nu - m = int Z dR (apply the self-adjoint (I+K_Q) to Z = (I+K_Q)^(-1)1), and int Z dR = -int Z' Cum R dx, so |nu - m| <= Var(Z) sup|Cum R|.

Cumulative error. Let k(y) = int_0^y K (odd, |k| <= q/2). Integrate (I+K_Q)nu' = R over [-Q, x]. Integrating by parts in u, with G(-Q) = 0 and G(Q) = nu - m,

    int_(-Q)^x (K_Q nu')(y) dy = int [k(x-u) - k(-Q-u)] dnu'(u)
                               = (K_Q G)(x) - (K_Q G)(-Q) + (nu - m)[k(x-Q) + k(2Q)].

So G + K_Q G = S + c_0 - (nu - m) psi on [-Q,Q], with S = Cum R, c_0 = (K_Q G)(-Q), psi(x) = k(x-Q) + k(2Q) in [0, k(2Q)] subset [0, q/2]. Hence G = g_1 + c_0 Z with g_1 = (I+K_Q)^(-1)[S - (nu-m)psi], and G(-Q) = 0 gives c_0 = -g_1(-Q)/Z(-Q), |c_0| Z(x) <= ||g_1|| (1+q). So ||G|| <= (2+q)||g_1|| and ||g_1|| <= [(1+q_J)/(1-q_J)] (s + q_J s q/2)/L. The value 45/4 = 3 * 3 * 5/4 is the gamma -> 0 limit, and each factor is monotone in gamma (constants_rep.py evaluates g(gamma); sup = 11.24999). End.

This replaces C139 8.2's ||G|| <= 27/(4L) (sawtooth 1/2, with the cruder psi <= q) and C140's unchecked rescaling.

[Numerical] L||G|| is between 0.500 and 0.532, L|nu - m| <= 0.057 (proved q_J s up to 0.49).

### Lemma I.7 (density) [Proved]

On [-Q,Q], sigma - r = -K*nu' = -int G(y) K'(x-y) dy (G = 0 below -Q and G = nu - m above Q, and K vanishes at infinity), so

    |sigma - r| <= ||K'||_1 ||G|| <= D/L,      D := 2g/(pi gamma)  (>= 2K(0) g = 2 cot(gamma) g/pi),
    |sigma' - r'| <= ||K''||_1 ||G|| <= D1/L,   D1 := 4 M_1 g,

with ||K'||_1 = 2K(0) = 2cot(gamma)/pi and ||K''||_1 = 4 max|K'| <= 4M_1, M_1 = (2/pi)(0.325/sin^2 gamma + 0.25/sin gamma) (C139 Lemma 5.0). D is chosen as an upper bound that does not vanish at gamma = pi/2, so that the ratio D1/D stays bounded there.

[Numerical] L sup|sigma - r| <= 0.92 against D up to 152 at Delta = 0.999, and L sup|sigma' - r'| <= 35 against D1 up to 4576.

### Lemma I.8 (quadrature) [Proved]

Assume L >= 4 gamma D = 8g/pi (so L >= 29 suffices at every gamma). Let A >= 0 solve r0(A) = 2D/L and call [-A,A] the core.

Core. On [-A,A] intersected with [-Q,Q], r >= r0 >= 2D/L, so sigma >= r - D/L >= r/2 > 0. Let x_(j1), ..., x_(j2) be the roots there (consecutive, z_L(x_j) = I_j/L) and f = eps o z_L^(-1). The composite trapezoid rule with step 1/L and the Peano kernel (t-a)(b-t)/2 <= 1/(8L^2) give

    |sum_(j1)^(j2) eps(x_j) - [eps(x_j1) + eps(x_j2)]/2 - L int_(x_j1)^(x_j2) eps sigma dx| <= V_f/(8L),

where V_f bounds the variation of f' = eps'/sigma = v g, g = r/sigma in (0, 2]. By the product rule for variation, Var(vg) <= ||v|| Var(g) + ||g|| Var(v) <= ab Var(g) + 2V_v. Writing r'sigma - r sigma' = r'(sigma - r) - r(sigma' - r') and using sigma^2 >= r^2/4, |r'| <= br and Lemma I.7, |g'| <= 4(bD + D1)/(L r). Since int_(-A)^A dx/r0 = (4gamma/b) sinh(bA) <= 2/(b r0(A)) = L/(bD),

    Var(g) <= 4 + 4 D1/(bD),   V_f = ab Var(g) + 2V_v <= 8ab + 80ab/pi + 4ab D1/(bD).

Here 4ab D1/(bD) = 4ab * 4M_1 g/((pi/gamma)(2g/(pi gamma))) = 8 ab gamma^2 M_1 = 8pi^2 gamma sin(gamma) M_1, and gamma sin(gamma) M_1 <= (2/pi)(0.325 gamma/sin gamma + 0.25 gamma) <= 0.575 (gamma/sin gamma <= pi/2 and gamma <= pi/2). So V_f <= 8pi^2 + 80pi + 45.4 <= 375.7, and the core contributes at most 46.97/L.

Ends. Integrating by parts on (x_j2, Q] with Cum R(x_j2^+) = 1/(2L) and eps(Q) = 0,

    sum_(j > j2) eps(x_j) - L int_(x_j2)^Q eps sigma dx + eps(x_j2)/2 = -L int_(x_j2)^Q eps' Cum R dx,

and similarly on the left (Cum R(x_j1^-) = -1/(2L)). By Lemma I.5 and |eps'| <= ab r, each end is at most s ab int_(x_j2)^Q r. If A < Q, there is no root in (x_j2, A], so int_(x_j2)^A r = G(x_j2^+) - G(A) <= 2g/L, and int_A^Q r <= c_r int_A^inf r0 <= c_r (pi/(2b)) r0(A) = c_r pi D/(bL) by Lemma I.1. If A >= Q, then j2 = n and int_(x_n)^Q r = G(x_n^+) - G(Q) <= 2g/L. With D/b = 2g/pi^2, both ends together are at most

    2 s ab g (2c_r/pi + 2)/L.

Empty core. If no root lies in [-A,A] (then A < Q), (I.1) gives |L int eps dR| <= s ab int_(-Q)^Q r <= s ab [2g + 2c_r pi D/b]/L, using int_(-A)^A r = G(-A) - G(A), which is below the ends bound. So no lower threshold on L beyond L >= 8g/pi is needed (C139 used L >= 77 to guarantee a nonempty core).

Curvature. |e_bulk(nu) - e_bulk(m) - h(m)(nu-m)| <= a2 (nu-m)^2/2 <= a2 (q_J s)^2/(2L^2). At Q = infinity there is no mismatch.

### Theorem I.9 (repulsive sector theorem, uniform through Heisenberg) [Proved]

For every 0 <= Delta < 1, L >= 2 and 0 <= n <= L/2,

    |F_L(n) - L e_bulk(Delta, n/L)| <= K_U/L,
    K_large(gamma) = V_f/8 + 2 s ab g (2c_r/pi + 2) + a2 (q_J s)^2/2  (for L >= 29),
    K_large(gamma) <= 46.96 + 765.50 + 2.47 <= 815 (analytic), sup over a 200001-point grid = 811.30,
    K_U = max(29^2, 815) = 841.

Proof. n = 0 is exact. For L >= 29 add the core, the two ends (or the empty-core bound) and the curvature, using F_L(n) - L e_bulk(nu) = [F_L - L(E(Q) + h(nu-m))] - L[e_bulk(nu) - e_bulk(m) - h(m)(nu - m)] with E(Q) = e_bulk(m), h = e_bulk'(m). For L < 29 use |F_L - L e_bulk| <= L (bond spectral diameter 1, e_bulk in the bond interval, as in C139 Section 5.4). The analytic majorant bounds each term separately: V_f/8 <= 375.7/8 = 46.96 (Lemma I.8), ends <= 2 pi^2 (45/4)(2c_r/pi + 2) = 765.50 (s <= 1, ab <= pi^2, g <= 45/4), curvature <= 2pi^2 (1/2)^2/2 = pi^2/4. constants_rep.py evaluates the exact gamma-dependent expression on 200001 points of (0, pi/2]; its supremum, 811.30, is approached as gamma -> 0 (where V_f/8 -> 43.42). Pointwise: K_large = 53.6 at Delta = 0, 376 at gamma = 0.5, 629 at Delta = 63/64, 811.3 at Delta -> 1. End.

[Numerical] L|F_L(n) - L e_bulk| <= 0.91 in all 160 cases, a factor 890 below K_large.

### Corollary I.10 (repulsive U(1) class including Heisenberg) [Proved]

For every bond of the U(1) spin-1/2 class with 0 <= J_z/J <= 1, every field, DM angle and L >= 2,

    L(E_(2L) - 2E_L) <= C_rep M,   C_rep = 5 * 815 + pi^2/4 = 4077.5 <= 4078.

(With the exact gamma-dependent K_large evaluated on the grid, 5 * 811.30 + pi^2/4 = 4058.96; the stated 4078 uses only the analytic majorant.)

Proof. Matched doubling (2L, 2n): for L >= 29, L[F_(2L)(2n) - 2F_L(n)] <= L[K/(2L) + 2K/L] = (5/2)K, i.e. 5K in units of M_0 = J/2. For L < 29 seam gluing gives at most 4L <= 112 in units of M_0 (C136 Section 9). The C129/C136 reductions (DM angle to flux, field cancellation in matched sectors, flux cost pi^2 M/4, M >= M_0 for |Delta| <= 1) are sign-free and refereed (C136 REFEREE, Assembly). At Delta = 1, for fixed L, E_L, E_(2L) and M are continuous in the bond, so the bound extends from Delta < 1. End.

This supersedes C139 Theorem A's 345600 on [0, 63/64] and C140's sketched 19935 + pi^2/4. The constant is larger than C139's conditional 2215 + pi^2/4 because the proved sawtooth is s <= 1 instead of the assumed 1/2. In C139 8.2 the factor 27 in the end integral was 4||G||L; the direct bound above uses 2||G||L.

### I.11 Which variables are uniform at Delta -> 1

[Proved, bookkeeping] All statements are in the rapidity x of C139, where the kernel width is sin gamma and the root density scale is b = pi/gamma. Every constant above is expressed through gamma-free combinations: L r0(x_n)/b, D/b, gamma D, ab, gamma sin(gamma) M_1, q, q_J, s. The rescaled rapidity lambda = x/gamma (XXX rapidities at gamma -> 0) makes r0 and K converge, but no step needs that limit: Delta = 1 itself is reached by finite-matrix continuity at fixed L.

### I.12 Consistency checks

* Delta = 0 [Proved]: K = J = 0, q = q_J = 0, s = 1/2, g = 2 * 1 * (1/2) = 1, sigma = r0 = a1 and Cum R is the exact sawtooth 1/(2L). K_large(pi/2) = 53.6, consistent with (coarser than) the exact free constant 8 - 4 sqrt 2.
* gamma -> 0 [Proved]: theta_1 -> 1 - 2(2-ln2)/pi = 0.168 > 0, L r0(x_n)/b >= 1/(2pi); no constant diverges.
* L = 2 [Numerical, ed_check.json]: the end-to-end ED defect L(F_(2L)(2n) - 2F_L(n))/M_0 is 3.987, 3.9987, 4.000 at Delta = 0.99, 0.999, 1 (L <= 6), below 4078 and equal to the known supremum 4 at Heisenberg.
* ED [Numerical]: Bethe energies equal ED sector energies at flux pi per bond to 2.3e-14 for all 1 <= n <= L/2, L = 2 to 12, Delta = 0.99, 0.999.
* Casimir and Delta = 1/2: unchanged from C139 (no new energy identity here).

