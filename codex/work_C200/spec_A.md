# Task A: a supersolution for the half-filled XXZ Bethe equations near the ferromagnet, with a uniform O(1) counting error

You are asked for a complete, rigorous proof (every step shown, explicit constants) of the statements in Parts 1 and 2 below, or, where a statement is false as written, a corrected statement with proof. Everything you need is defined here. Do not assume anything not stated. If a step can only be verified numerically, say so explicitly and separate it from the rigorous part. Plain ASCII or LaTeX-free math is fine.

## 0. Setting and notation

Parameter eta in (0, pi/4]. gamma = pi - eta, b = pi/gamma in (1, 4/3], q = 1 - 2 eta/pi.

Functions of a real variable x:

    p(x)     = 2 atan( tan(eta/2) tanh x ),            increasing, p(+-inf) = +-eta
    theta(x) = -2 atan( tanh x / tan eta ),            odd, decreasing, theta(+-inf) = -+(pi - 2 eta)
    a1(x)    = p'(x)/(2 pi) = sin(eta) / (pi (cosh 2x + cos eta)) > 0,     int_R a1 = eta/pi
    T(x)     = -theta'(x)/(2 pi) = sin(2 eta) / (pi (cosh 2x - cos 2 eta)) > 0,   even, strictly decreasing in |x|,
               int_R T = q,  T(0) = cot(eta)/pi
    Theta(x) = int_0^x T = (1/pi) atan( tanh x / tan eta ),   odd, increasing, Theta(+inf) = q/2 = 1/2 - eta/pi.
    Useful: q/2 - Theta(w) = int_w^inf T = (1/pi) [ atan( tan(eta) coth w ) - eta ]  for w > 0.

Continuum density at half filling (it solves rho = a1 + T * rho on R, where * is convolution; you may check this by Fourier transform, T^(k) = sinh((pi - 2 eta) k/2)/sinh(pi k/2), a1^(k) = sinh(eta k/2)/sinh(pi k/2), rho^(k) = 1/(2 cosh(gamma k/2))):

    rho(v) = b / (2 pi cosh(b v)),      N(v) = int_{-inf}^v rho = (1/pi) atan( e^{b v} )  in (0, 1/2).

Bethe system. L >= 2 even, n = L/2, quantum numbers I_j = j - (n+1)/2, j = 1..n. The Bethe map on R^n is

    F_j(x) = L p(x_j) - sum_{k=1}^n theta(x_j - x_k) - 2 pi I_j.

Known (proved earlier, may be cited): F has a unique real zero x* = (x*_1 < ... < x*_n), symmetric under x -> -x.

Known comparison principle (proved earlier; proof: at an index j maximizing x_j - y_j > 0, use that p is strictly increasing and theta is nonincreasing): if x, y in R^n and F_j(y) >= F_j(x) for all j then y_j >= x_j for all j.

Sawtooth: s(u) = u - floor(u) - 1/2 for real u (so |s| <= 1/2, period 1, s(u) = u - 1/2 on [0,1)).

## Part 1 (the main lemma, pure calculus)

For y in R and m >= 0 define the matched distances
    w_+(m) in (0, inf]:  N(y + w_+) = N(y) + m   if m < 1/2 - N(y),   w_+ = +inf otherwise,
    w_-(m) in (0, inf]:  N(y - w_-) = N(y) - m   if m < N(y),         w_- = +inf otherwise,
(w_+-(0) = 0), and
    D_y(m) = Theta(w_+(m)) - Theta(w_-(m)),      with Theta(inf) = q/2.
Explicitly, with alpha = pi N(y) = atan(e^{b y}) in (0, pi/2) and phi = pi m:
    w_+ = (1/b) log( tan(alpha + phi)/tan(alpha) )  for phi < pi/2 - alpha,
    w_- = (1/b) log( tan(alpha)/tan(alpha - phi) )  for phi < alpha.
D_y(0) = 0 and D_y(m) = 0 for m >= max(N(y), 1/2 - N(y)). By the symmetry y -> -y, D_{-y} = -D_y.

Claim 1. There is an explicit constant C_1 such that for every eta in (0, pi/4] and every y in R, the total variation of D_y on [0, inf) is at most 2 C_1 eta. Numerically (floating point, my own scan over eta from 1e-5 to pi/4 and alpha from 1e-9 to pi/4): for y <= 0, D_y <= 0, D_y is decreasing then increasing (one sign change of the derivative), sup |D_y| <= (b/pi) tan(eta) cos(2 alpha) roughly, and sup over y of sup_m |D_y|/eta is 0.3183 (= 1/pi) as eta -> 0 and 0.139 at eta = pi/4. So the expected truth is TV(D_y) = 2 sup|D_y| <= (2/pi) b tan(eta), i.e. C_1 about b tan(eta)/(pi eta) <= 0.55. Any explicit C_1 is acceptable, but please make it as small as you can prove; C_1 <= 2/pi would be ideal, and say what you get.

A possibly useful reduction (check it): D_y = (1/pi)[atan(tan(eta) coth w_-) - atan(tan(eta) coth w_+)] = (1/pi) atan( tan(eta) Xi / (1 + tan^2(eta) coth(w_+) coth(w_-)) ) with Xi = coth w_- - coth w_+, so TV(D_y) <= (tan(eta)/pi) (TV(Xi) + sup|Xi|) provided the factor 1/(1 + tan^2 eta coth w_+ coth w_-) is monotone in m with values in (0,1]. Xi depends on eta only through the exponent 2/b in (3/2, 2): coth w_+ = (X^(2/b) + A^(2/b))/(X^(2/b) - A^(2/b)) with A = tan alpha, X = tan(alpha + phi). In the limit b = 1, alpha -> 0, phi = alpha u: Xi = -(4 - 2u^2)/(4 - u^2) for u < 1. In general Xi(0+) = -b cos(2 alpha). If you can prove unimodality of D_y (or of Xi) directly, do that; otherwise bound the total variation by any rigorous route. The uniformity in alpha -> 0 and alpha -> pi/2 and in eta -> 0 is the point.

Claim 2 (consequence, please prove). For every y in R, every real L > 0 and every real c,
    J(y) := int_R s( L (N(v) - N(y)) ) T(y - v) dv
satisfies |J(y)| <= (1/2) TV(D_y) <= C_1 eta.
Sketch to be made rigorous: splitting v = y +- w, and using s(-u) = -s(u) for non-integer u, J = int_0^inf T(w) [ s(tau_+(w)) - s(tau_-(w)) ] dw with tau_+-(w) = +-L(N(y +- w) - N(y)) >= 0; changing variables to tau = L m on each side, J = int_0^inf s(L m) dD_y(m) as a Riemann-Stieltjes integral (D_y is continuous and of bounded variation, constant for large m; note the two sides have different finite ranges of tau, and beyond its range each Theta(w_+-) is the constant q/2), hence |J| <= (1/2) TV(D_y).

## Part 2 (the supersolution; please check every step and write the proof out rigorously)

Fix kappa0 >= 1/2 with kappa0 >= pi C_1/2. Let t(v) = L N(v) in (0, L/2) = (0, n). Define the trial configuration: for j = 1..n with t_j := j - 1/2 + kappa0 < n, let y_j be the unique real number with t(y_j) = t_j ("finite nodes"; they are increasing in j); for the remaining indices (the top ones; there are n_inf = number of j with j - 1/2 + kappa0 >= n, at least one) put y_j = +infinity. Define for finite j

    Ftilde_j = L p(y_j) - sum_{k finite} theta(y_j - y_k) - n_inf (pi - 2 eta) - 2 pi I_j,

i.e. the Bethe map with theta(y_j - y_k) replaced by its limit theta(-inf) = pi - 2 eta for the nodes at +infinity.

Claim 3 (comparison with nodes at infinity). If Ftilde_j >= 0 for every finite j, then x*_j <= y_j for every finite j. (Proof idea: at a finite index j maximizing x*_j - y_j > 0 over finite indices, theta(x*_j - x*_k) <= theta(y_j - y_k) for finite k as in the comparison principle, and theta(x*_j - x*_k) <= pi - 2 eta = sup theta for the other k; p strictly increasing gives F_j(x*) > Ftilde_j >= 0, a contradiction with F_j(x*) = 0.)

Claim 4 (residual identity). Let N_y(v) = (1/L) #{finite j : y_j <= v}. Then for finite j,
    Ftilde_j/(2 pi L) = kappa0/L - [T * (N - N_y)](y_j).
(Derivation to check: with z(u) = p(u)/(2 pi) + (1/L) sum_{k finite} Theta(u - y_k) - (n_inf/L)(q/2), Ftilde_j = 2 pi L (z(y_j) - I_j/L); integrating by parts, (1/L) sum_{k finite} Theta(u - y_k) = (T * N_y)(u) - (q/2) n_fin/L; the continuum function Z_0(u) = p(u)/(2 pi) + (T*N)(u) - q/4 equals N(u) - 1/4, which is the integrated form of rho = a1 + T*rho; and N(y_j) - 1/4 - I_j/L = kappa0/L.)

Claim 5 (clip inequality). L N_y(v) = max( floor( t(v) - kappa0 + 1/2 ), 0 ), and since kappa0 >= 1/2,
    L (N - N_y)(v) <= kappa0 + s( t(v) - kappa0 + 1/2 )   for every v,
with equality when t(v) >= kappa0 - 1/2.

Claim 6 (supersolution). Hence for finite j, using T >= 0, int T = q, and t(v) - kappa0 + 1/2 = t(v) - t(y_j) + j,
    L [T*(N - N_y)](y_j) <= q kappa0 + J(y_j) <= q kappa0 + C_1 eta <= kappa0,
because (1 - q) kappa0 = (2 eta/pi) kappa0 >= C_1 eta. So Ftilde_j >= 0 for all finite j, and by Claim 3, x*_j <= y_j for finite j.

Claim 7 (cumulative bound). Let N_x(v) = (1/L) #{j : x*_j <= v}. Then N_x >= N_y >= N - (kappa0 + 1/2)/L on R, and by the reflection symmetry of the roots and of N, |N_x(v) - N(v)| <= (kappa0 + 1/2)/L for every real v (state precisely what happens at the jump points).

Claim 8 (density bound). With sigma_L(u) = a1(u) + (1/L) sum_j T(u - x*_j),
    sigma_L - rho = T * (mu_x - rho dv),  and  |sigma_L(u) - rho(u)| <= 2 T(0) (kappa0 + 1/2)/L = (2 kappa0 + 1) cot(eta)/(pi L)  for every real u,
by one integration by parts (int |T'| = 2 T(0)). So eta L |sigma_L - rho| <= (2 kappa0 + 1) eta cot(eta)/pi <= (2 kappa0 + 1)/pi, uniformly in L and eta.

## What to return

1. A complete proof of Claim 1 with an explicit C_1 (this is the hard part; spend your effort here). State clearly which inequalities are proved analytically. If you must rely on a finite numerical check of an explicit function on a compact parameter set, isolate it as a clearly stated, explicit inequality so that it can be certified separately by interval arithmetic, and prove everything else.
2. Proofs of Claims 2 to 8, with any corrections. Flag any error in my sketches, including sign conventions and edge cases (integer points of the sawtooth, the ranges of tau, kappa0 = 1/2 exactly, small L such as L = 2, 4 where all nodes may be at infinity).
3. The final theorem with the smallest kappa0 your C_1 allows, uniform in eta in (0, pi/4] and even L >= 2.
