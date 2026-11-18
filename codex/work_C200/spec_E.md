# Task E: adversarial referee of a proof (XXZ Bethe roots near the ferromagnet, general filling)

You are the referee. Below are definitions, cited facts, and a text (Part 4) with statements and proofs. Find errors: wrong signs, wrong ranges, missing cases, constants that secretly depend on eta, Q, L or n, illegitimate uses of Harnack or of the reflection principle, boundary regularity gaps, off-by-one mistakes in counting, edge cases (nodes at Q, nodes beyond Q, nodes at infinity, n <= K0, delta(y) < eta, Q small or large, y = 0), and any step that does not follow. Re-derive each identity and estimate yourself. The harmonic-function argument of 4.2 and the Taylor/variation argument of 4.3 are the load-bearing new steps; spend most of your effort there, and say explicitly whether every constant called C is independent of eta in (0, 1/64], Q in (0, inf), y, L and n. For each numbered item say "correct" or give the precise defect and, if you can, the repair. Finish with a verdict per item: correct as stated, correct after listed repairs, or wrong. Tags such as [rigorous], (Astra), (Fable), (checked), [numerical] are bookkeeping and can be ignored; numerical remarks are not part of the proof.

## Definitions and cited facts (assume these)

eta in (0, 1/64] unless a statement says otherwise, gamma = pi - eta, q = 1 - 2 eta/pi.

    p(x)     = 2 atan(tan(eta/2) tanh x),   strictly increasing, p(+-inf) = +-eta
    theta(x) = -2 atan(tanh x/tan eta),     odd, strictly decreasing, theta(+-inf) = -+(pi - 2 eta)
    a1       = p'/(2 pi) = sin(eta)/(pi (cosh 2x + cos eta)),   Fourier transform sinh(eta k/2)/sinh(pi k/2)
    T        = -theta'/(2 pi) = sin(2 eta)/(pi (cosh 2x - cos 2 eta)) = sin(2 eta)/(2 pi (sinh^2 x + sin^2 eta)),  even, decreasing in |x|,
               Fourier transform sinh((pi - 2 eta) k/2)/sinh(pi k/2), int T = q, T(0) = cot(eta)/pi
    Theta(x) = int_0^x T,  Theta(inf) = q/2.

Integers L >= 2, 1 <= n <= L/2, I_j = j - (n+1)/2. Bethe map F_j(x) = L p(x_j) - sum_{k=1}^n theta(x_j - x_k) - 2 pi I_j on R^n.
Cited fact 1: F has a unique real zero x* (increasing), symmetric under x -> -x. N_x(v) = (1/L) #{j : x*_j <= v}, sigma_L = a1 + (1/L) sum_j T(. - x*_j).
Cited fact 2: for 0 < Q < inf, the equation r = a1 + int_{-Q}^{Q} T(. - v) r(v) dv on [-Q, Q] has a unique continuous solution r = r_Q, which is positive, even, nonincreasing on [0, Q], and satisfies r_Q <= r_0 := 1/(2 gamma cosh(pi x/gamma)); m(Q) = int_{-Q}^{Q} r_Q increases continuously from 0 to 1/2. r_0 solves the same equation on R, and N(v) := int_{-inf}^v r_0.
Cited fact 3 (Theorem H, half filling, proved separately for all eta in (0, pi/2)): for L even and n = L/2, |#{j : x*_j <= v} - L N(v)| <= 1 for all real v.
Sawtooth s#(u) = u - floor(u) - 1/2; {u} = u - floor(u).

## Part 1 of the text (for reference; already refereed, you need not recheck it)

# Part 1. The comparison principle with nodes at infinity

[rigorous] (C149 Theorem 1.1, restated; checked) If x, y in R^n and F_j(y) >= F_j(x) for all j, then y_j >= x_j for all j.

[rigorous] (Fable; Astra confirmed; checked) Lemma 1.1. Let a trial configuration have finite nodes y_j, j in a set Fin, and the remaining n_inf = n - |Fin| nodes at +infinity. For j in Fin put

    Ft_j = L p(y_j) - sum_{k in Fin} theta(y_j - y_k) - n_inf (pi - 2 eta) - 2 pi I_j.

If Ft_j >= 0 for every j in Fin, then x*_j <= y_j for every j in Fin.

Proof. Suppose not, and let j in Fin maximize x*_j - y_j > 0 over Fin. For k in Fin, x*_j - x*_k >= y_j - y_k, and theta is nonincreasing, so theta(x*_j - x*_k) <= theta(y_j - y_k). For k not in Fin, theta(x*_j - x*_k) <= sup theta = pi - 2 eta. p is strictly increasing, so p(x*_j) > p(y_j). Hence 0 = F_j(x*) > Ft_j >= 0, a contradiction. No condition is needed at the infinite nodes. End.


## The text to referee

# Part 4. General filling

Here 1 <= n < L/2, nu = n/L, and eta in (0, 1/64] (the restriction is for convenience in the disc arguments of 4.2). Q = Q_nu is the continuum Fermi point, m(Q) = nu, r = r_Q = (I - T_Q)^(-1) a1 on [-Q, Q], and s(u) = a1(u) + int_{-Q}^{Q} T(u - v) r(v) dv is its extension to R (s = r on [-Q, Q]; this is C149's s_nu). N_0 is the cumulative function of r 1_[-Q,Q], Z_0 = p/(2 pi) + T*N_0 - q nu/2 the continuum counting function (Z_0' = s, Z_0 = N_0 - nu/2 on [-Q, Q]), t(v) = L(Z_0(v) + nu/2), so t(-Q) = 0, t(Q) = n. B(d) = int_d^inf T = q/2 - Theta(d). delta(y) = Q - |y|, rho(y) = min(1, eta + delta(y)). C denotes universal constants, not explicit.

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

Proof. V(x, h) = (1/eta) int_h^{h + eta} U(x, z) dz is positive and harmonic for 0 < h < pi/2 - eta, and eta V_h(x, h) = U(x, h + eta) - U(x, h). As h -> 0, V and V_h extend continuously to the open segment {h = 0, |x| < Q}, with V_h(x, 0) = (s(x) - f(x))/eta = 0. So the even reflection of V in h is weakly harmonic, hence harmonic, and positive, in a neighbourhood of that segment that contains every disc of radius R <= min(delta(y), 1/4) centred at (y, 0). Two standard facts for a positive harmonic W on a disc of radius R are used: W is comparable to its central value on the concentric disc of radius R/2 (Harnack), and |grad^j W| <= C_j R^(-j) W(centre) on the disc of radius R/4.

Fix y and R = min(delta(y)/2, 1/4). If R < 16 eta then rho(y) <= 33 eta, and the disc estimates for U on the disc of radius eta/2 about (y, eta), which lies in the strip, give |d_x^j s(y)| <= C eta^(-j) s(y); this also proves |s'| <= C s/eta on R. If R >= 16 eta: V(y, eta) is comparable to s(y), because U(y, z) for eta <= z <= 2 eta is comparable to U(y, eta). Harnack and the derivative estimates for the reflected V on the disc of radius R about (y, 0) give, for 0 <= h <= R/4,

    V(y, h) <= C s(y),      |d_x^j d_h V(y, h)| <= C R^(-j-1) s(y),   j = 1, 2.

Let N = floor(R/(4 eta)), so R/8 <= N eta <= R/4. The same vertical comparison gives U(y, N eta) <= C V(y, N eta) <= C s(y), and the disc estimates for U about (y, N eta) with radius N eta give |d_x^j U(y, N eta)| <= C R^(-j) s(y). Finally U(x, eta) = U(x, N eta) - eta sum_{l=1}^{N-1} V_h(x, l eta) identically in x; differentiating j times at x = y, the sum costs at most N eta C R^(-j-1) s(y) <= C R^(-j) s(y). Since R is comparable to rho(y) in this case, the theorem follows. End.

[numerical] (Fable, check_reg.py) For 1 + Delta in {1e-2, 1e-3, 1e-4} and Q in {0.3, 1, 3}: sup |r'| rho/r <= 1.19 and sup |r''| rho^2/r <= 1.71.

Consequence. Integrating the first bound, r(v)/r(y) lies between (rho(v)/rho(y))^C and its inverse for v, y on the same side of 0, and s(Q + u) >= e^(-C) r(Q) for 0 <= u <= eta.

## 4.3 The logarithm-free quadrature bound

[rigorous] (Astra, astra_B.out Section 3; checked line by line; numerically in tv_gen.py) Theorem 4.3. For eta <= 1/64, Q in (0, inf) and y in [-Q, Q], TV(D) <= C eta/rho(y), hence |S(y)| <= C_S eta/rho(y) for every L.

Proof. If delta(y) < eta, |S| <= (Theta(Q - y) + Theta(Q + y))/2 <= q/2 and eta/rho >= 1/2. Let delta(y) >= eta, r = r(y), rho = rho(y) <= 2 delta(y), and use the rescaled mass u = m/r. By Theorem 4.2 and its consequence there is a universal c > 0 such that for 0 <= u <= u_0 = c rho neither branch is capped and w(u) = w_+-(r u) satisfies w(0) = 0, w' = r/r(y +- w) comparable to 1, w'' = -r^2 a'/a^3, w''' = -r^3 a''/a^4 + 3 r^3 a'^2/a^5 with a(w) = r(y +- w), so |w''| <= C/rho, |w'''| <= C/rho^2 and

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

[rigorous] (Astra and Fable independently; checked) Proposition 4.7 (interlacing with half filling). For every eta in (0, pi/2), L >= 2 and 1 <= n <= L/2, with M = ceil(L/2), h = M - n, and X_1 < ... < X_M the half-filled roots at length 2M:

    X_j <= x*_j <= X_{j+h}   (j = 1..n),      sup_v |N_x(v) - N_0(v)| <= (1 + h)/L.

Proof. Use y_j = X_{j+h} as a trial for (L, n). Subtracting the half-filled equation for the index j + h (the quantum numbers agree, I^(M)_{j+h} = I^(n)_j + h/2),

    F_j(y) = (L - 2M) p(y_j) + sum_{k=1}^h theta(y_j - X_k) + pi h >= -(2M - L) eta + 2 h eta >= 0   for h >= 1,

because p < eta, 2M - L is 0 or 1, and theta >= -(pi - 2 eta). So x*_j <= X_{j+h} by the comparison principle, and x*_j >= X_j by reflection (x*_j = -x*_{n+1-j} >= -X_{n+1-j+h} = X_j). Counting: L N_x(v) >= #{X <= v} - h >= 2M N(v) - 1 - h >= L N_0(v) - 1 - h, by Theorem H at length 2M (valid for all eta < pi/2 by Remark (iv)) and r_Q <= r_0 (C136 (2.2), refereed). Reflection gives the upper bound. End.

[rigorous] (Fable) Corollary 4.8 (density near half filling). In the setting of Proposition 4.7, for every real u, sigma_L(u) >= r_0(u) - (h + 2) cot(eta)/(pi L).

Proof. #{x* <= v} - 2M N(v) lies in [-h - 1, 1] for every v (from X_j <= x*_j <= X_{j+h} and Theorem H), an interval of length h + 2. One integration by parts, with the constant -h/2 subtracted (int T' = 0), gives |(1/L) sum_j T(u - x*_j) - (2M/L)(T*r_0)(u)| <= (h + 2) T(0)/L. Finally a1 + (2M/L) T*r_0 >= a1 + T*r_0 = r_0. End.

