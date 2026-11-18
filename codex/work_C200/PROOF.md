# C200 proof notes: the ferromagnetic strip, by a counting bound for the Bethe roots at every filling

2026-11-18. Manager: a Claude Fable 5.2 agent. Derivations and referee passes delegated to Astra (gpt-6.2-astra via the shim, seven calls, transcripts `astra_A.out` to `astra_G.out`, each ending with the trailer model=gpt-6.2-astra). No step in this file comes from gpt-6-astra. J = 1 throughout.

Tags. [rigorous] is a complete mathematical deduction written out here. [rigorous given cited X] uses an earlier run X with the status stated. [numerical] is floating point without enclosures. [physics-level] has an unproved step (none is used below). Each step also carries its source: (Astra) produced by gpt-6.2-astra, (Fable) produced by the managing agent, (checked) means the managing agent verified the step line by line, and where stated also numerically with its own script. Referee passes by separate Astra 6.2 instances are named where they apply (astra_C.out for Parts 1 to 3, astra_E.out for 4.1 to 4.4, astra_F.out for 4.5 and Corollary 4.12, astra_G.out for Corollary 4.13).

Notation (that of C149 and C151). eta = pi - gamma, Delta = -cos eta, b = pi/gamma, a = 2/b = 2 - 2 eta/pi in Part 2 only (elsewhere a = pi sin eta, as in C136), q = 1 - 2 eta/pi.

    p(x)     = 2 atan(tan(eta/2) tanh x),      theta(x) = -2 atan(tanh x/tan eta),
    a1       = p'/(2 pi) = sin(eta)/(pi (cosh 2x + cos eta)),
    T        = -theta'/(2 pi) = sin(2 eta)/(pi (cosh 2x - cos 2 eta)) = sin(2 eta)/(2 pi (sinh^2 x + sin^2 eta)),
    Theta(x) = int_0^x T = (1/pi) atan(tanh x/tan eta),   Theta(inf) = q/2,   int_R T = q,   T(0) = cot(eta)/pi.

Bethe map on R^n, I_j = j - (n+1)/2: F_j(x) = L p(x_j) - sum_k theta(x_j - x_k) - 2 pi I_j. Its unique real zero x* is the ground state of the sector (C134, C136 refereed; uniqueness also C149 Corollary 1.2). sigma_L = a1 + (1/L) sum_j T(. - x*_j) is the counting density, N_x(v) = (1/L) #{j : x*_j <= v}. Sawtooth s#(u) = u - floor(u) - 1/2.

Half filling means L even and n = L/2. There the continuum density is r_0(v) = b/(2 pi cosh(b v)) and N(v) = int_{-inf}^v r_0 = (1/pi) atan(e^{b v}), with values in (0, 1/2).

## Summary

- Route. C149's discrete comparison principle, with two new devices: trial nodes at +infinity (Lemma 1.1), which removes C149's constraint kappa0 < 1/2, and exterior trial nodes at explicit positions on the scale of the last continuum cell (Theorem 4.11).
- Referee status. An independent Claude Opus 5.5 agent refereed this run as C204 (codex/work_C204/REFEREE.md): no error, Theorem 4.11 and Corollaries 4.12 to 4.14 correct, Remark 4.15 correct as a sketch, constants universal and not explicit. Its four requested repairs (citation of C136 Section 2 for the monotonicity of r_Q, retagging of C149 Theorem 1.1 as background, renaming of the exterior node index, the justification of the even reflection in Theorem 4.2) are applied in this file.
- Theorem H [rigorous; refereed by a second Astra instance]. At half filling, for every eta in (0, pi/2) and every even L >= 2: |N_x(v) - N(v)| <= 1/L for all real v, and eta L |sigma_L - r_0| <= (2/pi) eta cot(eta) <= 2/pi at every point.
- Theorem 4.2 [rigorous; refereed]. |r_Q'| <= C r_Q/rho, |r_Q''| <= C r_Q/rho^2 with rho = min(1, eta + distance to the Fermi point), from a harmonic representation of the continuum equation.
- Theorem 4.3 [rigorous; refereed]. The sawtooth quadrature error of the quantile lattice of r_Q is at most C eta/rho, with no logarithm.
- Theorem 4.11 [rigorous; refereed]. At every filling, |N_x - N_0| <= K/L with a universal K (not explicit), for eta <= 1/64. Hence (Corollary 4.12) the logarithm-free (D) and (D1): eta L |sigma_L - s| <= 2K/pi and eta^2 L |sigma_L' - s'| <= 16K/pi^2 at every point, relative to the continuum Fermi point Q_nu.
- Corollary 4.13 [rigorous given the cited, refereed (I1) and (I3); refereed]. |F_L(n) - L e_bulk| <= K'/L with K' = 35 + 181 K + 9 K^2, every L and n. Its input Var(v_Q) <= 92.9 (C146 Theorem 2, which was not refereed) is proved again here, independently, by Astra.
- Corollary 4.14 [rigorous given 4.13 and the refereed doubling and reductions]. The strip is closed: L(E_2L - 2E_L) <= (5K' + pi^2/4) M for -1 < J_z/J <= -cos(1/64).
- Half filling by the C151 route has the explicit constant 70.2 (Corollary 3.2), and sectors of bounded magnetization h have K(h) <= 81 + 10h (Corollary 4.9); both rest on C151 Theorem 2.1, not refereed.
- Part 5 assembles the U(1) qubit class.

# Part 1. The comparison principle with nodes at infinity

[background, not used] (C149 Theorem 1.1, restated; C149 is not refereed, and nothing below depends on it, because Lemma 1.1 is proved here in full and is the only comparison statement used) If x, y in R^n and F_j(y) >= F_j(x) for all j, then y_j >= x_j for all j.

[rigorous] (Fable; Astra confirmed; checked) Lemma 1.1. Let a trial configuration have finite nodes y_j, j in a set Fin, and the remaining n_inf = n - |Fin| nodes at +infinity. For j in Fin put

    Ft_j = L p(y_j) - sum_{k in Fin} theta(y_j - y_k) - n_inf (pi - 2 eta) - 2 pi I_j.

If Ft_j >= 0 for every j in Fin, then x*_j <= y_j for every j in Fin.

Proof. Suppose not, and let j in Fin maximize x*_j - y_j > 0 over Fin. For k in Fin, x*_j - x*_k >= y_j - y_k, and theta is nonincreasing, so theta(x*_j - x*_k) <= theta(y_j - y_k). For k not in Fin, theta(x*_j - x*_k) <= sup theta = pi - 2 eta. p is strictly increasing, so p(x*_j) > p(y_j). Hence 0 = F_j(x*) > Ft_j >= 0, a contradiction. No condition is needed at the infinite nodes. End.

# Part 2. The quadrature lemma at half filling

For y in R and m >= 0 let w_+(m), w_-(m) in [0, inf] be the matched distances (w_+-(0) = 0),

    N(y + w_+) = N(y) + m  if m < 1/2 - N(y), else w_+ = inf;      N(y - w_-) = N(y) - m  if m < N(y), else w_- = inf,

and D_y(m) = Theta(w_+(m)) - Theta(w_-(m)), with Theta(inf) = q/2. D_y(0) = 0 and D_y(m) = 0 for m >= 1/2. Put alpha = pi N(y) = atan(e^{b y}) in (0, pi/2).

[rigorous] (Astra; checked line by line, and numerically in check_lemma.py) Lemma 2.1. For every eta in (0, pi/4] and y in R, D_y has one sign, |D_y| is nondecreasing and then nonincreasing in m (it vanishes identically at y = 0, and is zero after both distances saturate), and

    TV(D_y) = 2 sup_m |D_y(m)| <= |cos 2 alpha| sin(2 eta)/pi <= 2 eta/pi.

Proof. By the reflection y -> -y (which exchanges w_+ and w_-) take y <= 0, so 0 < alpha <= pi/4. Let phi = pi m in (0, pi/2), S = sin 2 alpha, d = cos 2 alpha in [0, 1), z = S cot phi, which decreases from +inf to 0.

Step 1 (coordinates). From N(v) = (1/pi) atan(e^{b v}), e^{b w_+} = tan(alpha + phi)/tan(alpha) while w_+ is finite, so

    coth(b w_+/2) = (tan(alpha + phi) + tan alpha)/(tan(alpha + phi) - tan alpha) = sin(2 alpha + phi)/sin(phi) = z + d,

using tan A +- tan B = sin(A +- B)/(cos A cos B). In the same way coth(b w_-/2) = sin(2 alpha - phi)/sin(phi) = z - d. At phi = pi/2 - alpha, z + d = cos(alpha)/cos(alpha) = 1, and at phi = alpha, z - d = 1, so w_+- = inf exactly when z +- d <= 1. Hence w_+- = a arccoth(z +- d) when z +- d > 1.

Step 2 (integral form). Define H(x) = Theta(a arccoth x) for x > 1 and H(x) = q/2 for x <= 1. H is continuous. For x > 1 put u = arccoth x, so du/dx = -sinh^2 u and

    h(x) := -H'(x) = a T(a u) sinh^2 u = (a sin(2 eta)/(2 pi)) sinh^2 u/(sinh^2(a u) + sin^2 eta),            (2.1)

and h := 0 for x <= 1. Since a > 1, h -> 0 as x -> 1+ (u -> inf) and as x -> inf (u -> 0), so h is continuous on R, and

    -D_y(m) = H(z - d) - H(z + d) = int_{z-d}^{z+d} h(x) dx >= 0.                                             (2.2)

Step 3 (h is unimodal). With t = sinh^2 u and g(t) = sinh^2(a u), g'(t) = a sinh(2 a u)/sinh(2 u). Its logarithmic derivative in u is 2a coth(2 a u) - 2 coth(2u) > 0, because r coth r is strictly increasing on r > 0 (its derivative has numerator sinh r cosh r - r > 0) and a > 1. So g is strictly convex in t. The ratio t/(g(t) + sin^2 eta) has derivative of the sign of g + sin^2 eta - t g', which equals sin^2 eta > 0 at t = 0 and has derivative -t g'' < 0. The ratio is positive and tends to zero at both ends, so it has exactly one critical point, a maximum. Hence h increases on (1, x_0] and decreases on [x_0, inf) for some x_0 > 1.

Step 4 (D_y is unimodal). B(z) := int_{z-d}^{z+d} h has B'(z) = h(z + d) - h(z - d), which is >= 0 for z + d <= x_0, <= 0 for z - d >= x_0, and decreasing in z in between (there h(z + d) is nonincreasing and h(z - d) nondecreasing). So B' changes sign at most once, from plus to minus. B(0) = 0 since h = 0 on [-d, d] (d < 1), and B -> 0 as z -> inf. So B, and with it |D_y| as a function of m, is nondecreasing up to a maximum and nonincreasing afterwards (for d = 0, B = 0). Therefore TV(D_y) = 2 sup|D_y| and sup|D_y| <= 2 d sup h.

Step 5 (the bound on h). Claim: sinh^2(a u) + sin^2 eta >= 2 a sinh^2 u for u >= 0. Expand sinh^2 x = sum_{k>=1} (2x)^{2k}/(2 (2k)!). The coefficient of u^{2k} in sinh^2(a u) - 2 a sinh^2 u is 2^{2k}(a^{2k} - 2a)/(2 (2k)!), which is -a(2 - a) for k = 1, (a^4 - 2a)/3 for k = 2, and nonnegative for k >= 2 since a^{2k} >= a^4 > 2a for a >= 3/2. Hence

    sinh^2(a u) - 2 a sinh^2 u >= -a(2 - a) u^2 + ((a^4 - 2a)/3) u^4 >= -a(2 - a) u^2 + (a^2/4) u^4,

the last step because a^2 - 2/a >= 9/4 - 4/3 = 11/12 > 3/4 for a >= 3/2. Concavity of sine gives sin eta >= 2 eta/pi = 2 - a, so

    sinh^2(a u) + sin^2 eta - 2 a sinh^2 u >= (2 - a)^2 - a(2 - a) u^2 + (a^2/4) u^4 = (2 - a - a u^2/2)^2 >= 0.

With (2.1), sup h <= sin(2 eta)/(4 pi), and sup|D_y| <= d sin(2 eta)/(2 pi). End.

Remarks. (i) The hypothesis eta <= pi/4 enters only through a >= 3/2 in Step 5. (ii) The constant is sharp: along alpha = eta -> 0 with z = eta^(-1/2), h/eta is asymptotic to 2/(pi(4 + eta^2 z^2)) on the window and |D_y|/eta -> 1/pi (referee computation, astra_C.out; tvD.py gives 0.31831). The limit depends on the rate: with z = eta^(-2) it is 0. (iv) [rigorous] (Astra referee, astra_C.out; checked) Steps 1 to 4 use only a > 1, so for every eta in (0, pi/2), TV(D_y) = 2 sup|D_y| <= q, because |D_y| <= q/2 trivially. For eta >= pi/4, q <= 2 eta/pi. Hence TV(D_y) <= 2 eta/pi and |J(y)| <= eta/pi (Lemma 2.2) hold for all eta in (0, pi/2), and Theorem H below holds on the whole attractive side 0 < eta < pi/2 with the same proof. The sharper factor |cos 2 alpha| sin(2 eta)/pi is proved only for eta <= pi/4. (iii) [numerical] (Fable, check_lemma.py) Step 5: the minimum of the ratio of the two sides over 400 values of eta and 4000 values of u is 1.0000000068. sup h/(sin 2 eta/(4 pi)) <= 0.9999993 and h is unimodal in 200 cases. (2.2) agrees with the direct formula to 9e-13 eta in 300 random cases.

[rigorous] (Astra and Fable; checked) Lemma 2.2. For every y in R and every real L > 0,

    J(y) := int_R s#(L (N(v) - N(y))) T(y - v) dv   satisfies   |J(y)| <= (1/2) TV(D_y) <= |cos 2 alpha| sin(2 eta)/(2 pi) <= eta/pi.

Proof. Split at v = y and use that T is even. With m_+(w) = N(y + w) - N(y) and m_-(w) = N(y) - N(y - w), both strictly increasing in w >= 0, and s#(-r) = -s#(r) except at integers (a finite set of w, of measure zero),

    J(y) = int_0^inf T(w) [s#(L m_+(w)) - s#(L m_-(w))] dw.

The image of the measure T(w) dw under m_+- has cumulative mass Theta(w_+-(m)) up to the end of the range of m_+-, and q/2 afterwards, with no atom. Subtracting, J(y) = int_{[0,1/2]} s#(L m) dD_y(m), a Riemann-Stieltjes integral with a continuous integrator of bounded variation and an integrand with finitely many jumps. Since |s#| <= 1/2, |J| <= TV(D_y)/2, and Lemma 2.1 applies. End.

[numerical] (Fable, check_lemma.py) J evaluated directly (cell by cell) in 2000 random cases (eta from 1e-3 to pi/4, L from 2 to 3000, alpha from 1e-4 to pi/2 - 1e-4): max |J|/(|cos 2 alpha| sin(2 eta)/(2 pi)) = 0.618.

# Part 3. Theorem H

[rigorous] (Fable for the construction; Astra checked and wrote out Claims 3 to 8 of spec_A; checked) Theorem H. Let eta in (0, pi/4] (and, by Remark (iv) of Part 2, any eta in (0, pi/2)), L >= 2 even, n = L/2. Then for all real u, v,

    |N_x(v) - N(v)| <= 1/L,        |sigma_L(u) - r_0(u)| <= 2 cot(eta)/(pi L).

Proof. Step 1 (the trial). For j = 1, ..., n - 1 let y_j be defined by L N(y_j) = j, and let y_n = +infinity. This is the quantile configuration of C149 with shift kappa0 = 1/2, which C149 excluded because its top level L N = n is not attained. Lemma 1.1 removes that restriction. (For L = 2 there is no finite node and Step 3 is vacuous.)

Step 2 (residual identity). For u real put z(u) = p(u)/(2 pi) + (1/L) sum_{k<n} Theta(u - y_k) - q/(2L), so that Ft_j = 2 pi L (z(y_j) - I_j/L) for j < n (theta = -2 pi Theta, pi - 2 eta = pi q, n_inf = 1). Let N_y(v) = (1/L) #{k < n : y_k <= v}. Since int_{y_k}^inf T(u - v) dv = q/2 + Theta(u - y_k),

    z(u) = p(u)/(2 pi) + (T*N_y)(u) - q/4.

Integrating r_0 = a1 + T*r_0 from -inf to u (Tonelli, p(-inf) = -eta) gives N(u) = (p(u) + eta)/(2 pi) + (T*N)(u), that is p/(2 pi) + T*N - q/4 = N - 1/4, because eta/(2 pi) + q/4 = 1/4. At a finite node, L N(y_j) - L/4 - I_j = j - n/2 - j + (n+1)/2 = 1/2. Therefore

    Ft_j/(2 pi L) = 1/(2L) - [T*(N - N_y)](y_j),     j < n.                                                   (3.1)

Step 3 (supersolution). With t(v) = L N(v) in (0, n), L N_y(v) = #{k <= n - 1 : k <= t(v)} = floor(t(v)), so L(N - N_y)(v) = t - floor(t) = 1/2 + s#(t(v)) for every v. At a finite node t(y_j) = j is an integer, so s#(t(v)) = s#(t(v) - t(y_j)) and, with int T = q and Lemma 2.2,

    L [T*(N - N_y)](y_j) = q/2 + J(y_j) <= q/2 + eta/pi = 1/2.

By (3.1), Ft_j >= 0 for all j < n, and Lemma 1.1 gives x*_j <= y_j for j < n.

Step 4 (cumulative bound). A finite node counted by N_y(v) has its root at or below v, so N_x(v) >= N_y(v) = floor(L N(v))/L >= N(v) - 1/L, and the same holds for the left limit N_x(v-) since N is continuous. The roots are symmetric, so N_x(v) = 1/2 - N_x((-v)-), and N(-v) = 1/2 - N(v). Applying the lower bound at -v gives N_x(v) <= N(v) + 1/L.

Step 5 (density). sigma_L - r_0 = T*(mu_x - r_0 dv) = int T'(u - v)(N_x - N)(v) dv after one integration by parts (the boundary terms vanish). int |T'| = 2 T(0) = 2 cot(eta)/pi, so |sigma_L - r_0| <= 2 cot(eta)/(pi L). End.

[numerical] (Fable, check_roots.py, check_roots.json) For 1 + Delta in {0.29, 1e-1, 1e-2, 1e-4, 1e-6, 1e-8} and L in {2, ..., 2048}, with Bethe roots from the C138 solver: the trial of Step 1 has min_j Ft_j/(2 pi eta) between 0.159 and 0.318 (proved >= 0); sup_v L|N_x - N| is between 0.50 and 0.765 (proved <= 1); eta L sup|sigma_L - r_0| is at most 0.318 (proved <= 0.637); min sigma_L/r_0 on {lambda >= 2} is 1.0000 to 1.0036.

[rigorous] (Fable) Corollary 3.1 (shell comparability at half filling). With lambda(u) = eta L r_0(u), sigma_L(u) >= (1 - 2/(pi lambda(u))) r_0(u). In particular sigma_L >= (1 - 1/pi) r_0 >= 0.68 r_0 on {lambda >= 2}, and sigma_L >= (3/4) r_0 on {lambda >= 8/pi}. So C145/C146's (G_c), C150's (H_I) with C_I <= 3.18/0.68, and C151's (S_c) hold at half filling, uniformly in eta in (0, pi/4] and L, with no logarithm.

[rigorous given cited C151 Theorem 2.1; C151 is not refereed, and the managing agent read its proof but did not re-derive its constants] Corollary 3.2 (sector theorem at half filling). For eta in (0, pi/4] and every even L >= 2,

    |F_L(L/2) - L e_bulk(Delta, 1/2)| <= K_H/L,     K_H = 70.2.

Derivation. C151's proof uses the shell hypothesis only through eta L sigma_L >= m(lambda) on the lifted segment and through r/sigma at one straddling cell and one connector per side (C151 Section 2.3). Theorem H gives eta L sigma_L >= lambda - 2/pi on all of {lambda >= M_0} and r_0/sigma_L <= 1/(1 - 2/(pi M_0)) there, so C142 Theorem II.6 is not needed. At half filling Q = infinity, the mismatch term of C151 vanishes and r_Q = r_0. Evaluating C151's formula with these inputs (constants_C200.py, a copy of C151's K_of with m(lambda) = lambda - 2/pi) gives K = 70.15 at M_0 = 2 (tail 43.50, straddling cells 6.52, connectors 5.81, lifted core 14.33) and 69.98 at the optimal M_0 = 2.1. [numerical evaluation of explicit expressions]

[rigorous given cited C151 Part 6, C138 (4.4), C136 Section 9 and the C129 reductions; these were not rechecked here] Corollary 3.3. For the matched pair of half-filled sectors, L |F_{2L}(L) - 2 F_L(L/2)| <= (5/2) K_H <= 176, uniformly in eta in (0, pi/4] and even L. In the normalization of C151 Part 6 this is the constant 5 K_H + pi^2/4 <= 354 for the zero-magnetization sectors on even rings.

What Theorem H does not give. The c-d bound needs every filling, because the field is arbitrary. Theorem H covers the sector n = L/2 only. Part 4 treats the other fillings.

# Part 4. General filling

Here 1 <= n < L/2, nu = n/L, and eta in (0, 1/64] (the restriction is for convenience in the disc arguments of 4.2, see Remark 4.15). Q = Q_nu is the continuum Fermi point, m(Q) = nu, r = r_Q = (I - T_Q)^(-1) a1 on [-Q, Q], and s(u) = a1(u) + int_{-Q}^{Q} T(u - v) r(v) dv is its extension to R (s = r on [-Q, Q]; this is C149's s_nu). N_0 is the cumulative function of r 1_[-Q,Q], Z_0 = p/(2 pi) + T*N_0 - q nu/2 the continuum counting function (Z_0' = s, Z_0 = N_0 - nu/2 on [-Q, Q]), t(v) = L(Z_0(v) + nu/2), so t(-Q) = 0, t(Q) = n. B(d) = int_d^inf T = q/2 - Theta(d). delta(y) = Q - |y|, rho(y) = min(1, eta + delta(y)). C denotes universal constants, not explicit.

## 4.1 The trial and the exact criterion

[rigorous] (Fable; Astra rechecked, astra_B.out Section 1; checked) Fix an integer K0 >= 0, kappa = K0 + 1/2, and suppose n > K0 (if n <= K0 the bound |N_x - N_0| <= n/L <= K0/L is trivial). Trial: y_j is the point with t(y_j) = j + K0 when that level is attained by t, and +infinity otherwise. Nodes j <= n - K0 lie in [-Q, Q] (node n - K0 is at Q), the K0 top nodes lie beyond Q or at infinity. With N_y counting finite nodes and A = L(N_0 - N_y), the argument of Part 3, Step 2 gives for every finite node

    Ft_j/(2 pi) = L Z_0(y_j) - I_j - (T*A)(y_j) = kappa - (T*A)(y_j),

and direct counting gives A = 0 for v < -Q, A = K0 + {t} - (K0 - floor(t))_+ on [-Q, Q], and A = K0 - min(floor(t - n), K0) for v > Q. Put, for a node y in [-Q, Q],

    S(y) = int_{-Q}^{Q} T(y - v) s#(t(v) - t(y)) dv.

Since int T = q, Ft_j >= 0 at a node y = y_j in [-Q, Q] follows from

    (**)   S(y) <= (2 eta/pi) kappa + kappa B(Q + y) + int_{v > Q} T(y - v) [1/2 + min(floor(t(v) - n), K0)] dv.

At a finite node y >= Q (including y = Q) the pairing v = y +- w gives A(y - w) - kappa <= 1/2 and A(y + w) - kappa <= -1/2, so Ft_j/(2 pi) >= (1 - q) kappa > 0 with no further condition. If (**) holds at every node in [-Q, Q), Lemma 1.1 gives x*_j <= y_j at the finite nodes, hence N_x >= N_y >= N_0 - (K0 + 1)/L, and by reflection |N_x - N_0| <= (K0 + 1)/L.

[rigorous] (Fable and Astra; checked) As in Lemma 2.2, with w_+-(m) the matched distances for r (equal r-mass m on the two sides of y, capped at the distances Q - y and Q + y to the edges) and D(m) = Theta(w_+(m)) - Theta(w_-(m)),

    S(y) = int_0^inf s#(L m) dD(m),   so |S(y)| <= TV(D)/2,

and TV(D) depends on eta, Q, y only, not on L.

## 4.2 A harmonic representation and the regularity of r_Q

[rigorous] (Astra, astra_B.out Section 2; checked line by line; numerically in check_reg.py) Lemma 4.1. Let P_h(x) = sin(2h)/(pi (cosh 2x - cos 2h)) for 0 < h < pi/2. Then P_h(x) = -(1/pi) Im coth(x + ih) is positive and harmonic in (x, h), int P_h dx = 1 - 2h/pi, P_h -> delta as h -> 0, P_eta = T, and P_{pi/2 - eta} = sin(2 eta)/(pi (cosh 2x + cos 2 eta)). With f = r 1_[-Q,Q] and g(x) = 1/(2 eta cosh(pi x/eta)), the function

    U(x, h) = (P_h * f)(x) + (P_{pi/2 - h} * g)(x)

is positive and harmonic in the strip 0 < h < pi/2, U(x, h) -> f(x) as h -> 0 for |x| < Q, and U(x, eta) = s(x).

Proof. coth(x + ih) = (sinh 2x - i sin 2h)/(cosh 2x - cos 2h) gives the first claims. For the last, P_eta * f = T * f, and in Fourier variables the transform of P_{pi/2 - eta} is sinh(eta k)/sinh(pi k/2) (the a1 transform with 2 eta in place of eta) and that of g is 1/(2 cosh(eta k/2)), whose product is sinh(eta k/2)/sinh(pi k/2), the transform of a1. So P_{pi/2 - eta} * g = a1 and U(., eta) = a1 + T*f = s. End. [numerical] The identity P_{pi/2 - eta} * g = a1 holds to 4e-16 relative at five points for three values of eta.

So the continuum equation r = s on [-Q, Q] reads U(x, 0) = U(x, eta) for |x| < Q.

[rigorous] (Astra; checked) Theorem 4.2 (regularity on the distance-to-edge scale). For eta <= 1/64, every Q in (0, inf) and every y in [-Q, Q],

    |r'(y)| <= C r(y)/rho(y),      |r''(y)| <= C r(y)/rho(y)^2,      and |s'| <= C s/eta on R.

Proof. V(x, h) = (1/eta) int_h^{h + eta} U(x, z) dz is positive and harmonic for 0 < h < pi/2 - eta, and eta V_h(x, h) = U(x, h + eta) - U(x, h). As h -> 0, V and V_h extend continuously to the open segment {h = 0, |x| < Q}, with V_h(x, 0) = (s(x) - f(x))/eta = 0. So the even reflection Vt(x, h) = V(x, |h|) is weakly harmonic: for a smooth test function phi supported in a disc centred on the segment, Green's identity on {h > e} gives int_{h>e} V Lap(phi) = int_{h=e} (phi V_h - V phi_h) dx, and on the reflected half the same with the opposite sign of the phi_h term; as e -> 0 the V_h terms vanish because V_h -> 0 locally uniformly on the segment, and the two V phi_h terms cancel because V is continuous up to the segment, so int Vt Lap(phi) = 0. By Weyl's lemma Vt is harmonic, and it is positive, in a neighbourhood of that segment that contains every disc of radius R <= min(delta(y), 1/4) centred at (y, 0). Two standard facts for a positive harmonic W on a disc of radius R are used: W is comparable to its central value on the concentric disc of radius R/2 (Harnack), and |grad^j W| <= C_j R^(-j) W(centre) on the disc of radius R/4.

Fix y and R = min(delta(y)/2, 1/4). If R < 16 eta then rho(y) <= 33 eta, and the disc estimates for U on the disc of radius eta/2 about (y, eta), which lies in the strip, give |d_x^j s(y)| <= C eta^(-j) s(y); this also proves |s'| <= C s/eta on R. If R >= 16 eta: V(y, eta) is comparable to s(y), because U(y, z) for eta <= z <= 2 eta is comparable to U(y, eta). Harnack and the derivative estimates for the reflected V on the disc of radius R about (y, 0) give, for 0 <= h <= R/4,

    V(y, h) <= C s(y),      |d_x^j d_h V(y, h)| <= C R^(-j-1) s(y),   j = 1, 2.

Let N = floor(R/(4 eta)), so R/8 <= N eta <= R/4. The same vertical comparison gives U(y, N eta) <= C V(y, N eta) <= C s(y), and the disc estimates for U about (y, N eta) with radius N eta give |d_x^j U(y, N eta)| <= C R^(-j) s(y). Finally U(x, eta) = U(x, N eta) - eta sum_{l=1}^{N-1} V_h(x, l eta) identically in x; differentiating j times at x = y, the sum costs at most N eta C R^(-j-1) s(y) <= C R^(-j) s(y). Since R is comparable to rho(y) in this case, the theorem follows. End.

[numerical] (Fable, check_reg.py) For 1 + Delta in {1e-2, 1e-3, 1e-4} and Q in {0.3, 1, 3}: sup |r'| rho/r <= 1.19 and sup |r''| rho^2/r <= 1.71.

Consequence. Integrating the first bound, |log(r(v)/r(y))| <= C |int_y^v dx/rho(x)| for v, y in [-Q, Q]; in particular r(v) is comparable to r(y) when |v - y| <= rho(y)/4 (rho is 1-Lipschitz, so rho >= 3 rho(y)/4 on that segment). Integrating |s'| <= C s/eta, s(Q + u) >= e^(-C) r(Q) for 0 <= u <= eta. (An earlier draft stated a global ratio bound in terms of rho(v)/rho(y) alone; the referee, astra_E.out, showed it is false when rho saturates at 1. It is not used.)

## 4.3 The logarithm-free quadrature bound

[rigorous] (Astra, astra_B.out Section 3; checked line by line; numerically in tv_gen.py) Theorem 4.3. For eta <= 1/64, Q in (0, inf) and y in [-Q, Q], TV(D) <= C eta/rho(y), hence |S(y)| <= C_S eta/rho(y) for every L.

Proof. If delta(y) < eta, each branch Theta(w_+-(m)) is nondecreasing up to its cap, so TV(D) <= Theta(Q - y) + Theta(Q + y) <= q, and eta/rho >= 1/2. Let delta(y) >= eta, r = r(y), rho = rho(y) <= 2 delta(y), and use the rescaled mass u = m/r. Both segments of length rho/4 from y stay in [-Q, Q], and r is comparable to r(y) on them by the consequence of Theorem 4.2, so each contains r-mass at least c' r(y) rho. Hence there is a universal c > 0 such that for 0 <= u <= u_0 = c rho neither branch is capped and w(u) = w_+-(r u) satisfies w(0) = 0, w' = r/r(y +- w) comparable to 1, w'' = -r^2 a'/a^3, w''' = -r^3 a''/a^4 + 3 r^3 a'^2/a^5 with a(w) = r(y +- w), so |w''| <= C/rho, |w'''| <= C/rho^2 and

    w = u + A u^2 + O(u^3/rho^2),   w' = 1 + 2 A u + O(u^2/rho^2),   A = A_+- = -+ r'(y)/(2 r(y)),   |A| <= C/rho.

For 0 <= u <= 1, |T^(j)(u)| <= C T(u) (u + eta)^(-j), j = 1, 2 (the denominator sinh^2 u + sin^2 eta is comparable to (u + eta)^2, with first derivative O(u) and second derivative O(1)). Taylor expansion gives

    T(w) w' = T(u) + A (u^2 T(u))' + O(u^2 T(u)/rho^2).                                                       (4.1)

The two first-order corrections, from the kernel and from the Jacobian, are kept together. Now u^2 T(u) = (sin(2 eta)/(2 pi)) (u/sinh u)^2 [sinh^2 u/(sinh^2 u + sin^2 eta)] is a product of a decreasing and an increasing function with values in [0, 1], so its total variation is at most sin(2 eta)/pi, and int_0^{u_0} u^2 T <= (sin(2 eta)/(2 pi)) u_0. The T(u) terms of the two branches cancel in dD, so the variation of D over 0 <= u <= u_0 is at most |A_+ - A_-| sin(2 eta)/pi + C eta u_0/rho^2 <= C eta/rho. Beyond u_0 each branch has variation at most its remaining kernel mass, at most int_{c' rho}^inf T <= C eta/rho, by T(w) <= eta/(pi w^2). End.

The logarithm that every earlier route met comes from estimating 2 u T and u^2 T' separately in (4.1); together they are the derivative of a function of bounded variation.

[numerical] (Fable, tv_gen.py) For 1 + Delta in {1e-2, 1e-3, 1e-4}, Q in {0.3, 1, 2.5, 5} and 160 points y per case: TV(D) min(1, delta)/eta <= 0.88.

## 4.4 Regimes in which the counting bound is proved

[rigorous] (Astra, astra_B.out Section 4; checked) Lemma 4.4. For d >= 0, eta/min(1, eta + d) <= C [eta + B(d + eta)]. (For d + eta < 1, B(d + eta) >= int_{d+eta}^{2(d+eta)} T >= c eta/(d + eta).)

[rigorous] (Astra; checked) Theorem 4.5 (dense edge). There are universal K0 and Lambda such that for eta <= 1/64 and all L, n with lambda_Q := eta L r_Q(Q) >= Lambda,

    |N_x(v) - N_0(v)| <= (K0 + 1)/L   for all real v.

Proof. By Theorem 4.3 and Lemma 4.4, K0 can be chosen so that |S(y)| <= kappa [2 eta/pi + B(Q + y) + B(Q - y + eta)] for all y in [-Q, Q] (the nearer edge is at distance delta(y), and B is decreasing). By the consequence of Theorem 4.2, t(Q + eta) - n = L int_Q^{Q+eta} s >= e^(-C) lambda_Q. With Lambda = e^C (K0 + 1) the bracket in (**) equals kappa for v >= Q + eta, so the right exterior contributes at least kappa B(Q - y + eta), and (**) holds at every node in [-Q, Q]. Apply 4.1. End.

[rigorous] (Astra; checked) Proposition 4.6 (tiny Fermi sea). If Q <= eta, then |N_x - N_0| <= K/L for a universal K and every L. (B(Q + y) >= B(2 eta) >= c and |S| <= q/2.)

[rigorous] (Astra and Fable independently; checked) Proposition 4.7 (interlacing with half filling). For every eta in (0, pi/2), L >= 2 and 1 <= n <= L/2, with M = ceil(L/2), h = M - n, and X_1 < ... < X_M the half-filled roots at length 2M (at n = L/2 put Q = infinity and N_0 = N; then h = 0, x* = X, and the statement is Theorem H):

    X_j <= x*_j <= X_{j+h}   (j = 1..n),      sup_v |N_x(v) - N_0(v)| <= (1 + h)/L.

Proof. Use y_j = X_{j+h} as a trial for (L, n). Let h >= 1. Subtracting the half-filled equation for the index j + h (the quantum numbers differ by h/2, I^(M)_{j+h} = I^(n)_j + h/2),

    F_j(y) = (L - 2M) p(y_j) + sum_{k=1}^h theta(y_j - X_k) + pi h >= -(2M - L) eta + 2 h eta >= 0   for h >= 1,

because p < eta, 2M - L is 0 or 1, and theta >= -(pi - 2 eta). So x*_j <= X_{j+h} by the comparison principle, and x*_j >= X_j by reflection (x*_j = -x*_{n+1-j} >= -X_{n+1-j+h} = X_j). Counting: L N_x(v) >= #{X <= v} - h >= 2M N(v) - 1 - h >= L N_0(v) - 1 - h, by Theorem H at length 2M (valid for all eta < pi/2 by Remark (iv)) and r_Q <= r_0 (C136 (2.2), refereed). Reflection gives the upper bound. End.

[rigorous] (Fable) Corollary 4.8 (density near half filling). In the setting of Proposition 4.7, for every real u, sigma_L(u) >= r_0(u) - (h + 2) cot(eta)/(pi L).

Proof. #{x* <= v} - 2M N(v) lies in [-h - 1, 1] for every v (from X_j <= x*_j <= X_{j+h} and Theorem H), an interval of length h + 2. One integration by parts, with the constant -h/2 subtracted (int T' = 0), gives |(1/L) sum_j T(u - x*_j) - (2M/L)(T*r_0)(u)| <= (h + 2) T(0)/L. Finally a1 + (2M/L) T*r_0 >= a1 + T*r_0 = r_0. End.

[rigorous given cited C151 Theorem 2.1, not refereed] Corollary 4.9 (sectors of bounded magnetization). Since r_0 >= r_Q for every Q, Corollary 4.8 gives eta L sigma_L >= lambda - (h + 2)/pi with lambda = eta L r_{Q_L}, which is the shell hypothesis of C151 with any M_0 > (h + 2)/pi. C151's formula (constants_C200.py) then gives, for eta in (0, pi/4], every L and h = ceil(L/2) - n,

    |F_L(n) - L e_bulk(Delta, n/L)| <= K(h)/L,   K(0) = 80.4, K(1) = 90.6, K(2) = 100.4, K(3) = 109.9, K(5) = 128.3, K(10) = 172.3, K(20) = 256.0,

so K(h) <= 81 + 10 h on this range. The bound is uniform in eta and L but grows with the magnetization h, so by itself it does not close the strip for fields that produce a finite magnetization density. Theorem 4.11 below does, with a constant that is not explicit; this corollary keeps an explicit one for small h.

## 4.5 Every filling: exterior nodes on the scale of the last cell

Theorem 4.5 used the trial nodes beyond Q at the quantile positions of t, which are far from Q when the edge is sparse, so that the right exterior contributes only B(Q - y)/2. The managing agent took this for an obstruction (spec_D.md, Section 4) on the strength of a miscount in the top node's own inequality. It is not one: the exterior nodes may be placed at explicit positions on the scale of the last continuum cell, and they are supersolution components there because the continuum counting function grows logarithmically beyond Q on the range d <= x <= min(1, Q)/4 that the construction uses (it saturates at t(inf) - n = (L - 2n) eta/(2 pi) far away).

[rigorous] (Astra, astra_D.out Section 1; checked) Lemma 4.10 (arbitrary exterior positions). Let k = K0 >= 0, kappa = k + 1/2, n > k. Keep the interior nodes t(y_j) = j + k, 1 <= j <= n - k, place m <= k nodes at points Q < z_1 < ... < z_m, and the other k - m nodes at +infinity. Then at every interior node y = y_j,

    Ft_j/(2 pi) >= (1 - q) kappa + kappa B(Q + y) + (1/2) B(Q - y) + sum_{l=1}^m B(z_l - y) - S(y).              (4.2)

Proof. The residual identity of 4.1 holds for any finite positions: Ft_j/(2 pi) = t(y_j) - n/2 - I_j - int T(y_j - v) D(v) dv with D = L N_0 - #{finite nodes <= v} (it follows from L p(u)/(2 pi) + int Theta(u - v) L r(v) 1_[-Q,Q] dv = t(u) - n/2, both sides having derivative L s and vanishing at 0, and Theta(-inf) = -q/2 for the nodes at infinity). At an interior node the first three terms give kappa. D = 0 below -Q, D <= kappa + s#(t(v)) on [-Q, Q], and D = k - sum_l 1{v >= z_l} above Q, with int_{v > Q} T(y - v) 1{v >= z_l} dv = B(z_l - y). Substituting, and using that t(y_j) is an integer, gives (4.2). End.

So each exterior node at a finite point contributes its full tail B(z_l - y).

[rigorous] (Astra and checked) Kernel bounds, with universal b, c > 0:

    b eta/x <= B(x) <= eta/(pi x)  (eta <= x <= 2; the upper bound for all x > 0),      (1 - q) + B(x) >= c eta/min(1, eta + x)  (x >= 0).      (4.3)

(T(w) <= eta/(pi w^2) for all w, T(w) >= c_1 eta/w^2 for eta <= w <= 4, integrate over [x, 2x]; for x < eta use B(x) >= B(eta) >= b, for x >= 1 use 1 - q = 2 eta/pi.)

[rigorous] (Astra, astra_D.out Sections 2 to 5; checked line by line; numerically in ext_trial.py) Theorem 4.11 (counting bound at every filling). There is a universal K such that for every eta in (0, 1/64], L >= 2 and 1 <= n < L/2,

    |N_x(v) - N_0(v)| <= K/L   for all real v.

Proof. Universal constants A, M, k are chosen in this order, with k >= M and k >= 2; C_S is the constant of Theorem 4.3. If n <= k the bound is trivial with K = k. Let n > k, so n >= 3 and the last continuum cell lies in the right half: define ell by L int_{Q - ell}^{Q} r = 1, and d = max(eta, ell), H = A d, R = min(1, Q).

Case 1: M H <= R/4. Use m = M exterior nodes at the real points z_l = Q + l H, l = 1, ..., M (l is an integer index, not the imaginary unit), and k - M nodes at infinity.

(a) Density at depth. r is nonincreasing on [0, Q] (C136 Section 2, refereed), so r(Q - ell) is at least the average of r over the last cell, 1/(L ell), and L r(Q - u) >= 1/ell >= 1/d for ell <= u <= Q.

(b) Logarithmic growth of the counting function beyond Q. Fix l, x = l H <= R/4 <= Q/4, E_l = t(Q + x) - n = L int_0^x s(Q + w) dw. Since a1 > 0 and s(Q + w) >= int_d^x T(u + w) r(Q - u) du, and for d <= u <= x, 0 <= w <= u one has eta <= u + w <= 2u <= 1/2 and T(u + w) >= c_1 eta/(4 u^2),

    E_l >= int_d^x L r(Q - u) int_0^u T(u + w) dw du >= (c_1 eta/(4 d)) int_d^x du/u = c_2 (eta/d) log(A l).

(c) The exterior nodes are supersolution components. The node z = z_l has index j = n - k + l. Writing Theta(w) = q/2 - B(w) at positive separations, the constant parts of the interaction sum to q I_j and

    Ft_j/(2 pi) = L p(z)/(2 pi) - (1 - q) I_j - sum_{a < j} B(z - y_a) + sum_{a > j, finite} B(y_a - z).

By the identity in the proof of Lemma 4.10 at u = z > Q, L p(z)/(2 pi) = E_l + (1 - q) n/2 + int B(z - v) L r(v) 1_[-Q,Q] dv. The function v -> B(z - v) is increasing, so on each of the n unit-mass cells of [-Q, Q] its value at the right endpoint exceeds its cell integral by at most its increment over the cell; summing, the sum over all n right endpoints exceeds the integral by at most B(z - Q) = B(l H). The interior nodes are the right endpoints of the cells k + 1, ..., n, so their sum is no larger. Dropping the nonnegative terms from nodes to the right, and using n/2 - I_j = k - l + 1/2 > 0 and B(l' H) <= eta/(pi l' A d),

    Ft_j/(2 pi) >= E_l - B(l H) - sum_{l'=1}^{l-1} B(l' H) >= (eta/d) [c_2 log(A l) - (1 + log l)/(pi A)].

Choose A >= e with 1/(pi A) <= c_2. The bracket is at least c_2 (1 + log l) - c_2 (1 + log l) = 0 for every l >= 1, independently of M.

(d) Interior nodes. Let 0 < y < Q be an interior node, delta = Q - y. Since t(Q) - t(y) is a positive integer, delta >= ell. If delta < 1, then rho(y) = min(1, eta + delta) >= delta and rho(y) >= d (d <= 1 by the case assumption), and eta <= z_l - y = delta + l A d <= 5/4, so by (4.3)

    sum_{l=1}^M B(z_l - y) >= b eta sum_{l=1}^M 1/(delta + l A d) >= (b eta/rho(y)) sum_{l=1}^M 1/(1 + l A).

The harmonic sum diverges, so M can be chosen (A being fixed) with b sum_{l <= M} 1/(1 + l A) >= C_S. Then S(y) <= C_S eta/rho(y) is covered in (4.2). If delta >= 1, rho(y) = 1 and (1 - q) kappa >= C_S eta for kappa >= pi C_S/2. If y <= 0, rho(y) = min(1, eta + Q + y) and kappa [(1 - q) + B(Q + y)] >= kappa c eta/rho(y) >= C_S eta/rho(y) for kappa >= C_S/c, by (4.3). At y = Q, S(Q) <= (1/2) int_{-Q}^{Q} T(Q - v) dv <= B(0)/2, which is the third term of (4.2).

Case 2: R < 4 A M d. Put all k exterior nodes at infinity (m = 0). The arguments for y <= 0, for delta >= 1 and for y = Q are unchanged. For 0 < y < Q with delta < 1, again delta >= ell, so rho(y) >= min(1, d), while min(1, eta + Q + y) <= min(1, eta + 2R) <= (1 + 8 A M) min(1, d) <= (1 + 8 A M) rho(y). By (4.3), kappa [(1 - q) + B(Q + y)] >= (c kappa/(1 + 8 A M)) eta/rho(y), which covers C_S eta/rho(y) once c kappa >= (1 + 8 A M) C_S. (When the Fermi sea is short on the scale of its last cell, the left exterior alone suffices.)

Choose k >= M with all the stated lower bounds on kappa. Then Ft_j >= 0 at every finite trial node, Lemma 1.1 gives x*_j <= y_j there, and N_x >= N_y >= N_0 - (k + 1)/L (on [-Q, Q] the trial count is max(0, floor(t) - k), above Q it is at least n - k). Reflection, with left limits at the jumps, gives the upper bound. K = k + 1. End.

[numerical] (Fable, ext_trial.py) With k = M = 4 and A = 2, exterior nodes at Q + 2 l d, l = 1..4, for 1 + Delta in {1e-3, 1e-4, 1e-5}, L in {64, 256, 1024} and five fillings (lambda_Q from 1.2 down to 0.0001, last cell from 0.7 eta to 240 eta): min over exterior nodes of Ft/(2 pi eta) is between 0.945 and 46.6, and min over interior nodes between 2.80 and 4.43. So the construction is a supersolution in every sampled case, with small explicit constants.

Together with Theorem H (K = 1, even L, n = L/2), the counting bound holds for every L and every 1 <= n <= L/2.

## 4.6 Consequences: logarithm-free (D) at every filling, and the strip

[rigorous] (Fable; this is C149 Theorem 3.3 with kappa0 + 1/2 replaced by K, re-derived) Corollary 4.12. For eta in (0, 1/64], every L and 1 <= n <= L/2, and every real u,

    eta L |sigma_L(u) - s(u)| <= 2 K eta cot(eta)/pi <= 2K/pi,      eta^2 L |sigma_L'(u) - s'(u)| <= 4 K eta^2 T_A(0) T_B(0) <= 16 K/pi^2,

where s = r_{Q_nu} on [-Q_nu, Q_nu]. Proof: sigma_L - s = T*(mu_x - mu_0) = int T'(u - v)(N_x - N_0)(v) dv and int |T'| = 2 T(0). For the derivative, T = T_A * T_B with T_A, T_B positive, even and unimodal (C149 Lemma 3.2: the Fourier transform sinh((pi - 2 eta) w/2)/sinh(pi w/2) factors as [sinh((pi - 2 eta) w/2)/sinh(gamma w/2)] [sinh(gamma w/2)/sinh(pi w/2)], each factor the transform of a kernel of the form sin(c)/(2 beta (cosh(pi x/beta) - cos c)); checked), so sigma_L' - s' = T_A' * [T_B * (mu_x - mu_0)], and each convolution costs twice the central value. This is the logarithm-free (D) and (D1) of C138, relative to the continuum Fermi point Q_nu, with sigma_L >= (1 - 2K/(pi lambda)) r_{Q_nu} wherever lambda = eta L r_{Q_nu} > 2K/pi.

Inputs for the energy, all continuum or cited statements. e(x) = -a a1(x) with a = pi sin(eta), so e = -sin^2(eta)/(cosh 2x + cos eta) and |e'| <= 2a a1.
 (I1) [cited, C136 Sections 1 and 5, refereed] The sector ground energy is F_L(n) = Delta L/4 + sum_j e(x*_j), and e_bulk(Delta, nu) = Delta/4 + int_{-Q}^{Q} e r_Q at m(Q) = nu.
 (I2) [rigorous, checked] With Z_Q = (I - T_Q)^(-1) 1 >= 1 and h = h(Q) = -a r_Q(Q)/Z_Q(Q), the function eps = eps_Q := (I - T_Q)^(-1)(e - h) = -a r_Q - h Z_Q is even, C^1, vanishes at +-Q, and satisfies |eps'| <= 2a r_Q on [-Q, Q]: differentiating eps = e - h + T_Q eps and integrating by parts (the boundary terms vanish because eps(+-Q) = 0) gives eps' = (I - T_Q)^(-1) e', and |e'| <= 2a a1 with the positivity of the resolvent gives the bound. So v := eps'/r_Q has |v| <= 2a.
 (I3) [cited, C145 Theorem B.1, refereed in C147; re-derived independently by the Astra referee, astra_G.out] 0 <= r_0 - r_Q <= r_0(2Q - x) + r_0(2Q + x) on [-Q, Q], hence r_Q(x) >= (1 - E) r_0(x) with E = 4 e^(-2b d), d = Q - |x|, and in particular r_Q >= kappa(d) r_0, kappa(d) = 1 - 4 e^(-2d); kappa := kappa(1) = 1 - 4 e^(-2) = 0.45865....
 (I4) [rigorous] (Astra, astra_G.out, an independent proof along the route of C146 Theorem 2, which was not refereed; checked line by line) Var_{[-Q,Q]}(v) <= 24 a + 4 a B <= 88 pi/3 < 92.2 for every eta in (0, pi/2] and Q in (0, inf], where B = 1/(1 - q) = pi/(2 eta); for eta <= pi/4 the bound is at most 73.1. Put V_W := 92.9, the value of C146, which is therefore valid.

Proof of (I4). Write R = (I - T_Q)^(-1), u = eps' = -a R a1', G = R[T(. - Q)], H = R[T(. + Q)] (nonnegative), Z = R 1. (i) For f = src + T_Q f, (I - T_Q) f' = src' - f(Q) T(. - Q) + f(-Q) T(. + Q). Applied to the odd u and the even r: u' = -a R a1'' - u(Q)(G + H), r' = R a1' - r(Q)(G - H), and therefore

    v' = F_1 - (r(Q)/r) [v(Q)(G + H) - v (G - H)],      F_1 = a [(R a1')^2 - r R a1'']/r^2.

(ii) H(x) = G(-x). On odd functions T_Q acts through the kernel T(x - y) - T(x + y) >= 0 on [0, Q]^2, so G - H >= 0 and u >= 0 on [0, Q]; with |a1'| <= 2 tanh(2|x|) a1, 0 <= v <= M_v := 2a tanh(2Q) <= 2a there. By self-adjointness int_0^Q (G + H) = int_{-Q}^{Q} G = <Z, T(. - Q)> = Z(Q) - 1, and Z <= B gives Z(Q) - 1 = (T_Q Z)(Q) <= B int_0^{2Q} T <= B q/2 = (B - 1)/2. (iii) a1''/a1 = 4(C^2 - cC - 2)/(C + c)^2 <= 4 with C = cosh 2x, c = cos eta, so R a1'' <= 4 r and F_1 >= -4a. Also, with w = r_0 and delta_k = w^(k) - R a1^(k) = R[int_{|y|>Q} T(. - y) w^(k)(y) dy], |w'| <= b w and |w''| <= b^2 w give |delta_k| <= b^k delta_0, and 0 <= delta_0 = w - r <= E w by (I3). Expanding N = (R a1')^2 - r R a1'' = N_0 - 2 w' delta_1 + delta_1^2 + delta_0 w'' + r delta_2 with N_0 = w'^2 - w w'' >= 0 gives N >= -4 b^2 w delta_0, so (F_1)_- <= 4 a b^2 E/(1 - E)^2 when E < 1. Using 4a for d <= d_0 = log(12)/(2b) and the second bound beyond (E from 1/3 down to 0), (1/a) int_0^Q (F_1)_- <= 2 log(12)/b + b < 6 for 1 <= b <= 2. (iv) v is odd with v(0) = 0, so Var_{[-Q,Q]} v = 2 v(Q) + 4 int_0^Q (v')_-, and (v')_- <= (F_1)_- + M_v (G + H) on [0, Q] (r(Q)/r <= 1 by the monotonicity of r on [0, Q], C136 Section 2, refereed; v(G - H) >= 0). Hence Var v <= 2 M_v + 24 a + 4 M_v (B - 1)/2 <= 24 a + 4 a B = 24 pi sin(eta) + 2 pi^2 sin(eta)/eta. At Q = inf, v = a b tanh(b x) and Var v = 2ab. End.

[rigorous given the cited (I1) and (I3), both refereed, and Theorems H and 4.11] (Fable, following C138 Theorem B, refereed in C145 and C147, and C149 Theorem 4.3, not refereed; written out in full, checked, and refereed by an Astra instance, astra_G.out, verdict: accept after the arithmetic and convention repairs, which are applied below) Corollary 4.13 (sector theorem on the strip). For eta in (0, 1/64], every L >= 2 and 0 <= n <= L/2,

    |F_L(n) - L e_bulk(Delta, n/L)| <= K'/L,      K' = 35 + 181 K + 9 K^2,

with K the constant of Theorem 4.11 (K >= 1).

Proof. n = 0 is exact. Let Q = Q_nu, r = r_Q, s its extension, and extend eps to R by epsh(u) = e(u) - h + int_{-Q}^{Q} T(u - v) eps(v) dv, which equals eps on [-Q, Q] (at half filling Q = inf, h = 0, eps = -a r_0, and the steps about Q_L and the boundary below are empty).

Step 1 (energy identity at Q_nu; C149 Lemma 4.1). On R, e - h = epsh - T*(eps 1_[-Q,Q]). Summing over the roots and using T*mu_x = sigma_L - a1, sum_j (e - h)(x*_j) = sum_j epsh(x*_j) - L int_{-Q}^{Q} eps (sigma_L - a1). Also L int_{-Q}^{Q} (e - h) r = L int [(I - T_Q) eps] r = L int eps a1, by the symmetry of T_Q. Subtracting, the h terms cancel because n = L m(Q), and with (I1)

    F_L(n) - L e_bulk(nu) = sum_j epsh(x*_j) - L int_{-Q}^{Q} eps sigma_L dx.

Step 2 (the extension). epsh' = e' + T*(eps' 1_[-Q,Q]) (integration by parts, eps(+-Q) = 0), so |epsh'| <= 2a (a1 + T*(r 1)) = 2a s on R. With Z_0 the continuum counting function (Z_0' = s, Z_0(+-Q) = +-nu/2), |epsh(u)| <= 2a |Z_0(u) - nu/2| for u >= Q, and |eps(u)| <= 2a (nu/2 - Z_0(u)) for 0 <= u <= Q; similarly on the left.

Step 3 (cells). z_L(u) = p(u)/(2 pi) + (T*N_x)(u) - q nu/2 is the counting function of the roots: z_L' = sigma_L > 0, z_L(x*_j) = I_j/L. Since z_L(+inf) - nu/2 = (eta/(2 pi))(1 - 2 nu), for nu < 1/2 there is a unique finite Q_L > 0 with z_L(Q_L) = nu/2; at nu = 1/2 put Q_L = inf, and the two outer cells extend to infinity (they are treated by the first bound of Step 5). The cells X_j = z_L^(-1)([I_j/L - 1/(2L), I_j/L + 1/(2L)]) tile [-Q_L, Q_L], each has sigma_L-mass 1/L, and x*_j is its midpoint in the variable t = z_L. Since z_L - Z_0 = T*(N_x - N_0) and |N_x - N_0| <= K/L,

    |z_L - Z_0| <= q K/L on R,   hence   |Z_0(Q_L) - nu/2| <= q K/L   and   |z_L(Q) - nu/2| <= q K/L.                 (4.4)

Write F_L - L e_bulk = sum_j [epsh(x*_j) - L int_{X_j} epsh sigma_L] + B_L with B_L = L int_{-Q_L}^{Q_L} epsh sigma_L - L int_{-Q}^{Q} eps sigma_L.

Step 4 (boundary term). B_L is an integral over the two intervals between Q and Q_L (in either order) and their mirrors. There |epsh| <= 2a |Z_0(Q_L) - nu/2| <= 2a q K/L by Step 2 and the monotonicity of Z_0, and the sigma_L-mass of each interval is |z_L(Q_L) - z_L(Q)| <= q K/L by (4.4). So |B_L| <= 2 L (2a q K/L)(q K/L) = 4 a q^2 K^2/L.

Step 5 (two bounds per cell). With f(t) = epsh(z_L^(-1)(t)), the cell term is f(t_j) - L int_cell f dt over a cell of length 1/L centred at t_j. It is at most the oscillation of f on the cell, hence at most int_{X_j} |epsh'| dx. It is also at most Var_cell(f')/(8L) (Peano kernel of the midpoint rule), and f' = epsh'/sigma_L in the variable x.

Step 6 (core; skipped if K_c is empty). Let C_D = 2K eta cot(eta)/pi and C_D1 = 4K eta^2 T_A(0) T_B(0) (Corollary 4.12), M = 4 C_D, rho_0 = M/(eta L), and K_c = {x in [-Q, Q] : r(x) >= rho_0} = [-x_c, x_c] (r is even and nonincreasing on [0, Q], C136 Section 2, refereed). On K_c, |sigma_L - r| <= C_D/(eta L) <= r/4, so g := r/sigma_L lies in [4/5, 4/3]. Core cells are those with X_j inside K_c. On them epsh'/sigma_L = v g, so the sum of their second bounds is at most (1/(8L)) Var_{K_c}(v g) <= (1/(8L)) [2a Var_{K_c}(g) + (4/3) V_W], by (I2) and (I4). Next g' = [r'(sigma_L - r) - r(sigma_L' - r')]/sigma_L^2 and sigma_L^2 >= (9/16) r^2, so

    Var_{K_c}(g) <= (16/9) [ (C_D/(eta L)) int_{K_c} |r'|/r^2 + (C_D1/(eta^2 L)) int_{K_c} dx/r ].

int_{K_c} |r'|/r^2 = 2 [1/r(x_c) - 1/r(0)] <= 2/rho_0 (monotonicity of r on [0, Q], C136 Section 2). For int_{K_c} dx/r, split [0, x_c] at Q - 1: the part with Q - x < 1 has length at most 1 and r >= rho_0 there; on the rest r >= kappa r_0 by (I3), int_0^y dx/r_0 = (2 gamma/b) sinh(b y) <= 1/(b r_0(y)), and r_0(x_c) >= r(x_c) >= rho_0, so that part is at most 1/(kappa b rho_0). Hence int_{K_c} dx/r <= 2 (1 + 1/(kappa b))/rho_0 and

    (2a/8) Var_{K_c}(g) <= (2a/9) [1 + C_D1 (1 + 1/(kappa b))/(eta C_D)],      a C_D1/(eta C_D) = 2 pi^2 sin(eta) tan(eta) T_A(0) T_B(0) <= 8 tan(eta)/eta.

Step 7 (the other cells). Their first bounds sum to at most int_U |epsh'|, U the union of the non-core cells. U is contained in the tail [-Q, Q] minus K_c, the K_c-parts of at most two cells that cross an endpoint of K_c (if one cell contains all of K_c it is charged once), and the intervals between Q and Q_L when Q_L > Q. On a straddling cell's K_c-part, r <= (4/3) sigma_L and the cell has sigma_L-mass 1/L, so int |eps'| <= 8a/(3L); two cells give 16a/(3L). On the tail r < rho_0: per side, the part within distance 1 of Q contributes at most rho_0, and on the rest r_0 < rho_0/kappa by (I3) and int_y^inf r_0 <= e^(-b y)/(gamma b) <= 2 r_0(y)/b, so int_tail r <= 2 rho_0 (1 + 2/(kappa b)) and int_tail |eps'| <= 4a rho_0 (1 + 2/(kappa b)). Beyond Q, int |epsh'| <= 2a int s = 4a |Z_0(Q_L) - nu/2| <= 4 a q K/L. If K_c is empty every cell is treated in this way and the same tail bound applies.

Step 8 (sum). Multiplying by L,

term by term, for eta <= pi/4 (tan(eta)/eta <= 4/pi, a <= pi sin(pi/4) < 2.2215, b >= 1, 1/kappa < 2.1803, q <= 1):
 core, g-part: (2a/9) + (2/9)(8 tan(eta)/eta)(1 + 1/kappa) <= 0.4937 + 7.1987;
 core, v-part: V_W/6 <= 15.484;
 tail: 4a M (1 + 2/(kappa b))/eta = 32 K cos(eta)(1 + 2/(kappa b)) <= 32 K (1 + 2/kappa) <= 171.54 K;
 straddling cells: 16a/3 <= 11.848;
 beyond Q and boundary: 4 a q K + 4 a q^2 K^2 <= 8.886 (K + K^2).
The sum is at most 35.1 + 180.44 K + 8.9 K^2 <= 35 + 181 K + 9 K^2, using K >= 1. End.

Remarks. (i) The core terms do not depend on K; the linear term comes from the first-order tail below the threshold lambda = M = 4 C_D and from the cells beyond Q. (ii) K is universal but not explicit, because the Harnack constants of 4.2 and the constants A, M, k of Theorem 4.11 were not tracked, so K' is finite and not explicit. In the cases computed, L sup|N_x - N_0| is at most 0.77 (check_roots.py, check_gen_roots.py); this is an observation, not a bound. (iii) At half filling K = 1 and the formula gives K' <= 225; Corollary 3.2, by the C151 route, gives 70.2.

[rigorous given Corollary 4.13 and the cited matched doubling and reductions (C138 (4.4), C136 Section 9, C129; refereed)] Corollary 4.14 (the strip). For every bond of the U(1) qubit class with -1 < J_z/J <= -cos(1/64), every field, DM angle and L >= 2,

    L (E_{2L} - 2 E_L) <= (5 K' + pi^2/4) M.

So the strip -1 < Delta < -1 + delta is closed for delta = 1 - cos(1/64) = 1.22e-4, with a universal constant that is finite and not explicit.

[rigorous as a sketch] (Fable; not seen by the Astra referees; refereed as a sketch by the independent Claude Opus 5.5 agent, C204, verdict correct) Remark 4.15 (the threshold eta <= 1/64 is a convenience). Lemma 4.1 holds for every eta in (0, pi/2). For eta in [1/64, pi/4] the disc of radius eta/2 about (y, eta) lies in the strip 0 < h < pi/2, so |d_x^j s(y)| <= C_j (eta/2)^(-j) s(y) <= 128^j C_j rho(y)^(-j) s(y) for every real y, because rho <= 1. So Theorem 4.2 holds for eta in (0, pi/4] with universal constants. The proofs of Theorem 4.3, of (4.3) and of Theorem 4.11 use the smallness of eta only through sin(2 eta) >= (4/pi) eta, sin(eta) >= (2 sqrt 2/pi) eta and eta <= 1, which hold on (0, pi/4]; in Theorem 4.11 large eta falls in Case 2. Corollaries 4.12 to 4.14 then hold for eta in (0, pi/4], that is for -1 < J_z/J <= -cos(pi/4) = -0.707. This extension changes only the size of the complementary constant in Part 5.

# Part 5. The U(1) qubit class, assembled

Bond h = J_perp(S^x S^x + S^y S^y) + D(S^x S^y - S^y S^x) + J_z S^z S^z - (B/2)(S^z x 1 + 1 x S^z), J = (J_perp^2 + D^2)^(1/2), Delta = J_z/J, M the gauge-optimized bond half-spread. For every field B, DM angle and length L >= 2:

| Delta = J_z/J | bound on L(E_2L - 2E_L)/M | source | status |
|---|---|---|---|
| Delta <= -1 | 0 | C129 | proved |
| -1 < Delta <= -cos(1/64) = -1 + 1.22e-4 | 5 K' + pi^2/4, K' = 35 + 181 K + 9 K^2, K the universal constant of Theorem 4.11 (finite, not explicit) | C200 Theorems H, 4.2, 4.3, 4.11 and Corollaries 4.12 to 4.14 | proved here: lemmas by Astra 6.2, constructions and assembly by Fable, every step checked by the managing agent and refereed by separate Astra 6.2 instances (verdicts correct, minor repairs applied); rests only on the refereed C129 (reductions), C136 (energy formulas, doubling, monotonicity of r_Q) and C145 Theorem B.1; refereed by an independent Claude Opus 5.5 agent (C204), correct after minor repairs |
| -cos(1/64) <= Delta <= 0 | 6e15 delta^(-25/3) at delta = 1.22e-4, that is 2.5e48 | C134, C136 | refereed |
| 0 <= Delta <= 1 | 4078 | C139, C140, C142 | refereed (C140 referee, C147) |
| Delta > 1 | 10781 (the referee's recomputed total is 10779.2, from K <= 4310.68 M_0) | C148 | refereed (C202), correct after minor repairs |

With Remark 4.15 the second row extends to -1 < Delta <= -cos(pi/4) and the third row becomes -cos(pi/4) <= Delta <= 0 with 6e15 (0.293)^(-25/3) = 1.7e20.

Consequence. C(2) restricted to the U(1) qubit class is finite, given the rows as stated: L(E_2L - 2E_L) <= C M with C the largest of the constants above. No row is conditional on an unproved hypothesis any more. The value of C is not explicit, because K is not, and it is dominated either by 5K' or by the C136 constant on the complementary range. Making K explicit requires tracking the Harnack and derivative constants in Theorem 4.2 and the constants A, M, k in Theorem 4.11; the numerical value of L sup|N_x - N_0| is at most 0.77 in every case computed, and at half filling K = 1 is proved.

Referee status of the chain: every row is now refereed. C148 on Delta > 1 was refereed as C202 (correct after minor repairs; largest computed sector error 2.45 M_0 at L = 2, not 1.82), and this run was refereed by an independent Claude Opus 5.5 agent as C204 (no error), in addition to the Astra referee passes recorded above. Remark 4.15 is refereed as a sketch only. The earlier unrefereed runs C146, C149, C150 and C151 are no longer needed for the strip: C146 Theorem 2 is re-proved as (I4), C149 Theorems 3.3 and 4.3 are re-derived as Corollaries 4.12 and 4.13, and C150 and C151 enter only the explicit side results (Corollaries 3.2, 3.3 and 4.9).
