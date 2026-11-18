# Task C: adversarial referee of a proof (half-filled XXZ Bethe roots near the ferromagnet)

You are the referee. Below is a self-contained statement and proof. Your job is to find errors: wrong signs, wrong ranges, missing cases, illegitimate interchanges, off-by-one mistakes in counting, edge cases (L = 2, L = 4, eta = pi/4, alpha = pi/4, integer points of the sawtooth, coincident points), and any step that does not follow. Re-derive every identity and inequality yourself rather than trusting the text. For each numbered step say "correct" or give the precise defect and, if you can, the repair. Finish with a verdict: correct as stated, correct after listed repairs, or wrong. Do not be polite about it. Tags such as [rigorous], (Astra), (Fable), (checked) in the text are bookkeeping and can be ignored.

Also answer two side questions with proofs or counterexamples:
 (Q1) Is the constant in Lemma 2.1 optimal, and is Theorem H's cumulative bound 1/L optimal in order (is sup_v L|N_x - N| bounded below by a positive constant uniformly)?
 (Q2) Does Lemma 2.1 survive for eta in (pi/4, pi/2) (a in (1, 3/2)) with the same or another explicit constant? Give the best statement you can prove. (This is optional; spend at most a fifth of your effort on it.)

## Definitions and cited facts (assume these)

eta in (0, pi/4], gamma = pi - eta, b = pi/gamma in (1, 4/3], a = 2/b = 2 - 2 eta/pi in [3/2, 2), q = 1 - 2 eta/pi.

    p(x)     = 2 atan(tan(eta/2) tanh x),   strictly increasing, p(+-inf) = +-eta
    theta(x) = -2 atan(tanh x/tan eta),     odd, strictly decreasing, theta(+-inf) = -+(pi - 2 eta)
    a1       = p'/(2 pi) = sin(eta)/(pi (cosh 2x + cos eta))
    T        = -theta'/(2 pi) = sin(2 eta)/(pi (cosh 2x - cos 2 eta)) = sin(2 eta)/(2 pi (sinh^2 x + sin^2 eta)),  even, decreasing in |x|
    Theta(x) = int_0^x T = (1/pi) atan(tanh x/tan eta),  Theta(inf) = q/2,  int_R T = q,  T(0) = cot(eta)/pi.

L >= 2 even, n = L/2, I_j = j - (n+1)/2 for j = 1..n. Bethe map on R^n: F_j(x) = L p(x_j) - sum_{k=1}^n theta(x_j - x_k) - 2 pi I_j.
Cited fact 1: F has a unique real zero x* = (x*_1 < ... < x*_n), and the set of roots is symmetric under x -> -x.
Cited fact 2: r_0(v) = b/(2 pi cosh(b v)) satisfies r_0 = a1 + T * r_0 on R (convolution), and N(v) = int_{-inf}^v r_0 = (1/pi) atan(e^{b v}).
sigma_L(u) = a1(u) + (1/L) sum_j T(u - x*_j), N_x(v) = (1/L) #{j : x*_j <= v}, s#(u) = u - floor(u) - 1/2.
You may check cited fact 2 by Fourier transform if you wish (T^(k) = sinh((pi - 2 eta)k/2)/sinh(pi k/2), a1^(k) = sinh(eta k/2)/sinh(pi k/2)).

## The text to referee

# Part 1. The comparison principle with nodes at infinity

[rigorous] (C149 Theorem 1.1, restated; checked) If x, y in R^n and F_j(y) >= F_j(x) for all j, then y_j >= x_j for all j.

[rigorous] (Fable; Astra confirmed; checked) Lemma 1.1. Let a trial configuration have finite nodes y_j, j in a set Fin, and the remaining n_inf = n - |Fin| nodes at +infinity. For j in Fin put

    Ft_j = L p(y_j) - sum_{k in Fin} theta(y_j - y_k) - n_inf (pi - 2 eta) - 2 pi I_j.

If Ft_j >= 0 for every j in Fin, then x*_j <= y_j for every j in Fin.

Proof. Suppose not, and let j in Fin maximize x*_j - y_j > 0 over Fin. For k in Fin, x*_j - x*_k >= y_j - y_k, and theta is nonincreasing, so theta(x*_j - x*_k) <= theta(y_j - y_k). For k not in Fin, theta(x*_j - x*_k) <= sup theta = pi - 2 eta. p is strictly increasing, so p(x*_j) > p(y_j). Hence 0 = F_j(x*) > Ft_j >= 0, a contradiction. No condition is needed at the infinite nodes. End.

# Part 2. The quadrature lemma at half filling

For y in R and m >= 0 let w_+(m), w_-(m) in (0, inf] be the matched distances,

    N(y + w_+) = N(y) + m  if m < 1/2 - N(y), else w_+ = inf;      N(y - w_-) = N(y) - m  if m < N(y), else w_- = inf,

and D_y(m) = Theta(w_+(m)) - Theta(w_-(m)), with Theta(inf) = q/2. D_y(0) = 0 and D_y(m) = 0 for m >= 1/2. Put alpha = pi N(y) = atan(e^{b y}) in (0, pi/2).

[rigorous] (Astra; checked line by line, and numerically in check_lemma.py) Lemma 2.1. For every eta in (0, pi/4] and y in R, D_y has one sign, |D_y| increases and then decreases in m, and

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

Step 4 (D_y is unimodal). B(z) := int_{z-d}^{z+d} h has B'(z) = h(z + d) - h(z - d), which is >= 0 for z + d <= x_0, <= 0 for z - d >= x_0, and decreasing in z in between (there h(z + d) is nonincreasing and h(z - d) nondecreasing). So B' changes sign at most once, from plus to minus. B(0) = 0 since h = 0 on [-d, d] (d < 1), and B -> 0 as z -> inf. So B, and with it |D_y| as a function of m, rises to a single maximum and falls back to zero. Therefore TV(D_y) = 2 sup|D_y| and sup|D_y| <= 2 d sup h.

Step 5 (the bound on h). Claim: sinh^2(a u) + sin^2 eta >= 2 a sinh^2 u for u >= 0. Expand sinh^2 x = sum_{k>=1} (2x)^{2k}/(2 (2k)!). The coefficient of u^{2k} in sinh^2(a u) - 2 a sinh^2 u is 2^{2k}(a^{2k} - 2a)/(2 (2k)!), which is -a(2 - a) for k = 1, (a^4 - 2a)/3 for k = 2, and nonnegative for k >= 2 since a^{2k} >= a^4 > 2a for a >= 3/2. Hence

    sinh^2(a u) - 2 a sinh^2 u >= -a(2 - a) u^2 + ((a^4 - 2a)/3) u^4 >= -a(2 - a) u^2 + (a^2/4) u^4,

the last step because a^2 - 2/a >= 9/4 - 4/3 = 11/12 > 3/4 for a >= 3/2. Concavity of sine gives sin eta >= 2 eta/pi = 2 - a, so

    sinh^2(a u) + sin^2 eta - 2 a sinh^2 u >= (2 - a)^2 - a(2 - a) u^2 + (a^2/4) u^4 = (2 - a - a u^2/2)^2 >= 0.

With (2.1), sup h <= sin(2 eta)/(4 pi), and sup|D_y| <= d sin(2 eta)/(2 pi). End.

Remarks. (i) The hypothesis eta <= pi/4 enters only through a >= 3/2 in Step 5. (ii) The constant is sharp: as eta -> 0 with alpha -> 0 and z -> inf, |D_y|/eta -> 1/pi (Astra's computation, and tvD.py gives 0.31831). (iii) [numerical] (Fable, check_lemma.py) Step 5: the minimum of the ratio of the two sides over 400 values of eta and 4000 values of u is 1.0000000068. sup h/(sin 2 eta/(4 pi)) <= 0.9999993 and h is unimodal in 200 cases. (2.2) agrees with the direct formula to 9e-13 eta in 300 random cases.

[rigorous] (Astra and Fable; checked) Lemma 2.2. For every y in R and every real L > 0,

    J(y) := int_R s#(L (N(v) - N(y))) T(y - v) dv   satisfies   |J(y)| <= (1/2) TV(D_y) <= |cos 2 alpha| sin(2 eta)/(2 pi) <= eta/pi.

Proof. Split at v = y and use that T is even. With m_+(w) = N(y + w) - N(y) and m_-(w) = N(y) - N(y - w), both strictly increasing in w >= 0, and s#(-r) = -s#(r) except at integers (a finite set of w, of measure zero),

    J(y) = int_0^inf T(w) [s#(L m_+(w)) - s#(L m_-(w))] dw.

The image of the measure T(w) dw under m_+- has cumulative mass Theta(w_+-(m)) up to the end of the range of m_+-, and q/2 afterwards, with no atom. Subtracting, J(y) = int_{[0,1/2]} s#(L m) dD_y(m), a Riemann-Stieltjes integral with a continuous integrator of bounded variation and an integrand with finitely many jumps. Since |s#| <= 1/2, |J| <= TV(D_y)/2, and Lemma 2.1 applies. End.

[numerical] (Fable, check_lemma.py) J evaluated directly (cell by cell) in 2000 random cases (eta from 1e-3 to pi/4, L from 2 to 3000, alpha from 1e-4 to pi/2 - 1e-4): max |J|/(|cos 2 alpha| sin(2 eta)/(2 pi)) = 0.618.

# Part 3. Theorem H

[rigorous] (Fable for the construction; Astra checked and wrote out Claims 3 to 8 of spec_A; checked) Theorem H. Let eta in (0, pi/4], L >= 2 even, n = L/2. Then for all real u, v,

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

