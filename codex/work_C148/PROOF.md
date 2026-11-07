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
