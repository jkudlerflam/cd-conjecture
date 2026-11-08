# C149 proof notes: the ferromagnetic strip by a discrete comparison principle

2026-11-07. Written by a Claude Opus agent, not by Codex. J = 1 throughout. Tags: [Proved] is a mathematical deduction, [Proved, cited] uses an earlier run (named, with its referee status), [Proved, conditional on (R_kappa)] is a deduction from the one hypothesis left, [Numerical] is floating point without enclosures, [Heuristic] has an unproved step, [Assumed] is a hypothesis.

Notation of C136 Section 1, C138 Section 1 and C145/C146. eta = pi - gamma in (0, pi/2], Delta = -cos eta, q = 1 - 2eta/pi, B = pi/(2eta), a = pi sin eta, b = pi/gamma. p(x) = 2 atan(tanh(x) tan(eta/2)), theta(x) = -2 atan(tanh(x) cot(eta)), a1 = p'/(2pi), T = -theta'/(2pi) >= 0, even, decreasing in |x|, with T(x) = sin(2eta)/(pi(cosh 2x - cos 2eta)), T(0) = cot(eta)/pi, ||T||_1 = q. For 1 <= n <= L/2, I_j = j - (n+1)/2, nu = n/L. The Bethe map is

    F_j(x) = L p(x_j) - sum_(k=1)^n theta(x_j - x_k) - 2 pi I_j,     x in R^n,

and the ground state roots x* are its unique real zero (C134, C136 Section 1; uniqueness is reproved in Corollary 1.2). For any configuration y in R^n, mu_y = L^(-1) sum_k delta_(y_k), N_y(v) = mu_y((-inf, v]), z_y(u) = (1/2pi)[p(u) - int theta(u - v) d mu_y(v)] is the counting function and sigma_y = z_y' = a1 + T*mu_y the counting density. For the true roots we write x, sigma_L, z_L, and Q_L for the counting boundary, z_L(Q_L) = nu/2 (C136 Section 4).

The continuum reference. Q = Q_nu is the continuum Fermi point, m(Q) = int_(-Q)^Q r_Q = nu (Q_nu = infinity at n = L/2). It is in general different from Q_L. mu_0 = r_Q 1_[-Q,Q] dx, N_0 its cumulative function, and

    Z_0(u) = (1/2pi)[p(u) - int theta(u - v) d mu_0(v)],     s_nu = Z_0' = a1 + T*mu_0 > 0 on R,

so s_nu = r_Q on [-Q, Q] (the restricted equation), Z_0(+-Q) = +-nu/2, and Z_0 is strictly increasing on R with Z_0(+-infinity) = +-[nu/2 + (eta/(2pi))(1 - 2nu)] (the continuum mass outside [-Q, Q] is int a1 + q nu - nu = (eta/pi)(1 - 2nu), split evenly). lambda_nu = eta L r_Q is the local crossover parameter of C145, but at Q_nu.

Summary.

- [Proved] Part 1. The comparison principle: F(y) >= F(x) componentwise implies y >= x, for all real configurations, all eta in (0, pi/2] and all real quantum numbers, with no ordering hypothesis. Uniqueness of the real Bethe solution is a corollary.
- [Proved] Part 2. An exact residual identity for arbitrary configurations, F_j(y)/(2pi L) = [Z_0(y_j) - I_j/L] - T*(N_0 - N_y)(y_j), and its form for the explicit trials y^(+-) = quantiles of Z_0 at the shifted levels I_j/L +- kappa0/L. The hypothesis left is (R_kappa): these trials are a super- and a subsolution.
- [Proved, conditional on (R_kappa)] Part 3. The sandwich gives the cumulative bound L sup|N_x - N_0| <= kappa0 + 1/2, and through the factorization T = T_A * T_B of the kernel into two positive unimodal kernels, log-free versions of C138's (D) and (D1) at every point of the real line: eta L |sigma_L - s_nu| <= (2kappa0 + 1) eta cot(eta)/pi and eta^2 L |sigma_L' - s_nu'| <= (2kappa0 + 1) 2 eta^2 T_A(0) T_B(0) <= (2kappa0 + 1) 8/pi^2.
- [Proved, conditional on (R_kappa) with kappa0 = 0.4] Part 4. With an energy identity at Q_nu (no density mismatch term), C138 Theorem B runs with these constants: |F_L(n) - L e_bulk| <= 136.4/L on the whole attractive side, and L(E_(2L) - 2E_L) <= 685 M. For C150's route, sigma_L >= 0.714 r_(Q_L) on its shell S' whenever Q_L <= Q_nu (observed in every sampled case).
- [Numerical] Part 5. (R_kappa) with kappa0 = 0.4 holds in every sampled case, 1 + Delta from 0.5 to 1e-6, L to 8192, all fillings including half filling. The worst margin is at half filling, about (2/pi)(kappa0 - b/4) eta/L, which is 0.093 eta/L at 1 + Delta = 1e-4. A proof of (R_kappa) is not given.

# Part 1. The comparison principle

[Proved] Theorem 1.1. Let eta in (0, pi/2], L > 0, and let I in R^n be any vector of quantum numbers (not necessarily ordered or half-integer). For x, y in R^n, if F_j(y) >= F_j(x) for every j, then y_j >= x_j for every j.

Proof. Suppose not, and let m = max_j (x_j - y_j) > 0, attained at some index j. For every k, (x_j - x_k) - (y_j - y_k) = m - (x_k - y_k) >= 0. theta is nonincreasing on R (theta' = -2pi T <= 0), so theta(x_j - x_k) <= theta(y_j - y_k) for every k (the term k = j is theta(0) = 0 on both sides). p is strictly increasing (p' = 2pi a1 > 0) and x_j > y_j. Hence

    F_j(x) - F_j(y) = L[p(x_j) - p(y_j)] - sum_k [theta(x_j - x_k) - theta(y_j - y_k)] > 0,

which contradicts F_j(y) >= F_j(x). End.

Remarks. (1) The Jacobian dF_j/dx_k = theta'(x_j - x_k) = -2pi T(x_j - x_k) <= 0 for k != j, and its row sums are 2pi L a1(x_j) > 0, so DF is a strictly diagonally dominant Z-matrix (a nonsingular M-matrix) at every configuration. F is off-diagonally antitone and strictly diagonally isotone, that is an M-function in the sense of Rheinboldt (J. Math. Anal. Appl. 32 (1970) 274; Ortega and Rheinboldt, Iterative Solution of Nonlinear Equations, Section 13.5), and Theorem 1.1 is the inverse isotonicity of M-functions. The proof above is self-contained and uses only two facts: theta depends on differences and is nonincreasing, p is strictly increasing. (2) The theta form of the map is defined and smooth on all of R^n, so no ordering of the roots and no restriction of the domain is needed. On ordered configurations it coincides with the step-free form L p(x_j) - sum_(k != j) phi(x_j - x_k) of C138 (2.2), which is the one evaluated numerically. (3) At the free point eta = pi/2, theta = 0 and the statement is trivial.

[Proved] Corollary 1.2. (i) The equation F(x) = 0 has at most one real solution, so the ground state roots x* are the unique real solution (C134 and C136 Section 1 prove existence). (ii) If F(y^+) >= 0 >= F(y^-) componentwise, then y^- <= x* <= y^+. (iii) If I' >= I componentwise and both systems have real solutions x(I), x(I'), then x(I') >= x(I).

Proof. (i) and (ii) are Theorem 1.1 applied twice. (iii) F^(I)(x(I')) = F^(I')(x(I')) + 2pi(I' - I) >= 0 = F^(I)(x(I)). End.

[Proved] Lemma 1.3 (counting functions). z_y(u) is nonincreasing in each y_k (dz_y(u)/dy_k = -T(u - y_k)/L <= 0). So y^- <= x <= y^+ gives z_(y^+) <= z_L <= z_(y^-) pointwise, and N_(y^+) <= N_x <= N_(y^-) pointwise.

# Part 2. The residual identity and the trial configurations

[Proved] Lemma 2.1 (exact residual identity). For every configuration y in R^n and every u in R,

    z_y(u) = Z_0(u) + [T*(N_y - N_0)](u),

and consequently, for every j,

    F_j(y)/(2 pi L) = [Z_0(y_j) - I_j/L] - [T*(N_0 - N_y)](y_j).                (2.1)

The same holds with mu_0 replaced by any finite measure of total mass nu with compact support.

Proof. G = N_y - N_0 is of bounded variation and vanishes outside a compact set (both cumulative functions are 0 far left and nu far right). Stieltjes integration by parts gives int theta(u - v) dG(v) = -int G(v) d_v[theta(u - v)] = -int G(v) 2pi T(u - v) dv, since d_v theta(u - v) = -theta'(u - v) dv = 2pi T(u - v) dv. Hence z_y - Z_0 = -(1/2pi) int theta(u - v) dG(v) = T*G. Finally F_j(y) = 2pi L[z_y(y_j) - I_j/L] directly from the definitions. End.

For the true roots, (2.1) reads Z_0(x_j) - I_j/L = T*(N_0 - N_x)(x_j). The first bracket of (2.1) is the counting displacement of node j relative to the continuum, the second is a T-average of the cumulative difference, which for a regular configuration is a T-smoothed sawtooth.

Definition (the trials). For kappa = kappa0/L with 0 < kappa0 < 1/2, let y^+ and y^- be defined by

    Z_0(y^(+-)_j) = I_j/L +- kappa.

They exist and are strictly increasing in j (Z_0 is a strictly increasing bijection of R onto the open interval of half-width nu/2 + (eta/(2pi))(1 - 2nu) >= nu/2, and every level satisfies |I_j/L +- kappa| <= nu/2 - (1/2 - kappa0)/L < nu/2 because kappa0 < 1/2). At even half filling Q = infinity and y^(+-)_j = (2/b) atanh(tan(pi(I_j/L +- kappa))) in closed form.

[Proved] Lemma 2.2 (residual of the trials). Let t(v) = L(Z_0(v) + nu/2). Then

    F_j(y^(+-))/(2pi L) = +- kappa - [T*G_(+-)](y^(+-)_j),
    L G_(+-)(v) = clip(t(v); 0, n) - clip(floor(t(v) -+ kappa0 + 1/2); 0, n),            (2.2)

with clip(s; 0, n) = min(max(s, 0), n). In the bulk (0 <= t <= n and 0 <= t - kappa0 + 1/2 < n + 1, for y^+), L G_+ = kappa0 + s(t - t_j) where s(u) = {u} - 1/2 is the zero-mean sawtooth and t_j = t(y^+_j) = j - 1/2 + kappa0 (so L G_+ = kappa0 - 1/2 just to the right of a node and kappa0 + 1/2 just to its left).

Proof. N_(y^+)(v) = L^(-1) #{j : I_j/L + kappa <= Z_0(v)} = L^(-1) #{j in [1, n] : j <= t(v) - kappa0 + 1/2}, and N_0(v) = clip(Z_0(v) + nu/2; 0, nu). Insert in (2.1) with Z_0(y^+_j) - I_j/L = kappa. In the bulk, t - floor(t - kappa0 + 1/2) = kappa0 - 1/2 + {t - kappa0 + 1/2} = kappa0 + s(t - t_j) since t_j - kappa0 + 1/2 = j is an integer. End.

[Assumed] Hypothesis (R_kappa). For the given kappa0, every eta in (0, pi/2], L >= 2 and 1 <= n <= L/2,

    F(y^+) >= 0 >= F(y^-)  componentwise,  equivalently  [T*G_(+-)](y^(+-)_j) <=  kappa (for +),  >= -kappa (for -).

By spin-flip and reflection symmetry (y^-_j = -y^+_(n+1-j), F_j(y^-) = -F_(n+1-j)(y^+)) the two halves are equivalent, so (R_kappa) is the single statement F(y^+) >= 0.

What (R_kappa) is. It is a statement about an explicit, perfectly regular configuration: the nodes are the quantiles of the continuum counting function at equally spaced levels. No Bethe root enters. In the variable t, (R_kappa) at node j reads

    int g(t) K_j(t) dt <= kappa0,    K_j(t) = T(y_j - V(t)) V'(t) >= 0,   int K_j = q,           (2.3)

with V the inverse of t(v) and g = L G_+ as in (2.2). In the bulk g = kappa0 + s(t - t_j), so (2.3) is the quadrature inequality

    int s(t - t_j) K_j(t) dt <= (1 - q) kappa0 + (end terms),                                 (2.4)

the first-order quadrature error of the dressing kernel against equally spaced nodes in the counting variable, compared with the uniform gain (1 - q) kappa0 = (2 eta/pi) kappa0 of a twist.

[Heuristic, confirmed numerically] Size of the left side of (2.4). Poisson summation in t, with the Lorentzian form of T on scale eta and rho = s_nu locally linear, gives

    int s(t - t_j) K_j(t) dt ~ L eta^2 rho'(y_j)/(e^(2 pi lambda_j) - 1),     lambda_j = eta L rho(y_j),

which is at most (eta/(2pi)) |rho'/rho| (sparse limit lambda -> 0) and exponentially small for lambda >> 1. It is positive where rho increases (the left half for y^+). At half filling |rho'/rho| <= b, so (2.4) needs kappa0 > b/4 in the bulk tails, and the margin at kappa0 is (2/pi)(kappa0 - b/4) eta (in units of 1/L). [Numerical] At half filling (explore_half.py, closed-form trial) the formula agrees with the computed residual to 2e-3 eta/L at 1 + Delta = 1e-2 and 2e-5 eta/L at 1 + Delta = 1e-6, and max_j L|residual|/eta -> b/(2pi). Near a Fermi edge rho'/rho is of order 1/(Q - |y|) down to the kernel scale, and the twist gain is larger there by the same profile, since (I - T_Q)1 >= T*1_out = (phi(d) - 2eta)/(2pi), about eta/(pi d) at distance d from the edge.

Where the amplification went. The spec warns that (I - T_Q)^(-1) amplifies uniform sources by B = pi/(2eta). Here the trial is built at the level of counting, not of density: a counting shift kappa buys the residual gain (1 - q)kappa = kappa/B, and the residual to beat is of order eta/L, so kappa of order B eta/L = O(1/L) suffices. The factor B is spent once, against a residual that is itself of order eta, and the outcome is an O(1) number of roots.

Why the constant shift works without a room condition. A uniform counting shift kappa on the full line costs nothing to the continuum part (it is a twist), and (2.1) is exact for every configuration, so the end effects are automatically contained in G_+. The only constraint is existence of the top node, kappa0 < 1/2 + L(Z_0(infinity) - nu/2) = 1/2 + eta(L - 2n)/(2pi), which holds for kappa0 < 1/2.

# Part 3. From the sandwich to log-free density bounds

Throughout this part assume (R_kappa) for the given (eta, L, n), so y^- <= x <= y^+ by Corollary 1.2.

[Proved, conditional on (R_kappa)] Lemma 3.1 (cumulative bound). For every v in R,

    |N_x(v) - N_0(v)| <= (kappa0 + 1/2)/L.

Proof. By Lemma 1.3, N_(y^+)(v) <= N_x(v) <= N_(y^-)(v). As in the proof of Lemma 2.2, L N_(y^(+-))(v) = clip(floor(t(v) -+ kappa0 + 1/2); 0, n), and L N_0(v) = clip(t(v); 0, n). clip(.; 0, n) is nondecreasing and 1-Lipschitz, and t - kappa0 - 1/2 <= floor(t - kappa0 + 1/2) and floor(t + kappa0 + 1/2) <= t + kappa0 + 1/2. So L N_0 - kappa0 - 1/2 <= L N_x <= L N_0 + kappa0 + 1/2. End.

This is C136's cumulative currency, with the constant kappa0 + 1/2 in place of 3B/2 (of order 1/eta), and with the continuum at Q_nu in place of Q_L. C138 Section 5.3 found numerically that the cumulative error relative to r_(Q_L) grows like eta^(-0.7) near a Fermi boundary; relative to Q_nu it stays below 0.65 in every sampled case (check_density.py), as Lemma 3.1 predicts.

[Proved] Lemma 3.2 (factorization of the kernel). For 0 < c < pi and beta > 0 the kernel k_(c,beta)(x) = sin(c)/(2 beta (cosh(pi x/beta) - cos c)) is positive, even, strictly decreasing in |x|, has Fourier transform sinh((pi - c) beta w/pi)/sinh(beta w), and k(0) = cot(c/2)/(2 beta). Consequently, with eta' = pi eta/gamma in (0, pi) for eta < pi/2,

    T = T_A * T_B,   T_A(x) = sin(eta')/(gamma (cosh(2 pi x/gamma) - cos eta')),   T_B(x) = sin(eta)/(pi (cosh 2x - cos eta)),

since T^(w) = sinh((pi - 2eta)w/2)/sinh(pi w/2) = [sinh((pi - 2eta)w/2)/sinh(gamma w/2)] [sinh(gamma w/2)/sinh(pi w/2)]. T_A(0) = cot(eta'/2)/gamma, T_B(0) = cot(eta/2)/pi, ||T_A'||_1 = 2T_A(0), ||T_B'||_1 = 2T_B(0).

Proof. The Fourier transform is the classical int e^(iwx) dx/(cosh x - cos c) = 2pi sinh((pi - c)w)/(sin(c) sinh(pi w)) after rescaling x -> pi x/beta (the same identity gives C138 (2.1) for T with (c, beta) = (2eta, pi/2) and for a1 with (pi - eta, pi/2)). T_B is (c, beta) = (eta, pi/2). T_A is (c, beta) = (eta', gamma/2): (pi - eta') gamma/(2pi) = (gamma - eta)/2 = (pi - 2eta)/2. Positivity and monotonicity are read off the formula (cosh is increasing in |x| and -cos c > -1). [Numerical] T_A * T_B = T to relative error 1.5e-11 at x in {0, eta, 5 eta} (check_density.py). End.

[Proved, conditional on (R_kappa)] Theorem 3.3 (log-free density and derivative bounds). For every u in R,

    eta L |sigma_L(u) - s_nu(u)| <= C_D^# := (2kappa0 + 1) eta cot(eta)/pi <= (2kappa0 + 1)/pi,
    eta^2 L |sigma_L'(u) - s_nu'(u)| <= C_D1^# := (2kappa0 + 1) 2 eta^2 T_A(0) T_B(0) <= (2kappa0 + 1) 8/pi^2.

On [-Q_nu, Q_nu], s_nu = r_(Q_nu) and s_nu' = r_(Q_nu)' (one-sided at the ends).

Proof. sigma_L - s_nu = T*(mu_x - mu_0) = int T(u - v) d(N_x - N_0)(v) = -int (N_x - N_0)(v) d_v T(u - v), so |sigma_L - s_nu| <= ||T'||_1 sup|N_x - N_0| = 2T(0)(kappa0 + 1/2)/L by Lemma 3.1. For the derivative, sigma_L' - s_nu' = T_A' * f with f = T_B*(mu_x - mu_0), and |f| <= 2T_B(0)(kappa0 + 1/2)/L by the same integration by parts, so |sigma_L' - s_nu'| <= ||T_A'||_1 ||f||_inf = 4T_A(0)T_B(0)(kappa0 + 1/2)/L. The suprema over eta: eta cot(eta) <= 1, and 2eta^2 T_A(0)T_B(0) = 2eta^2 cot(eta'/2) cot(eta/2)/(pi gamma) decreases from its eta -> 0 limit 8/pi^2 to 0 at eta = pi/2 (constants_C149.json; grid of 2e4 points). End.

Comparison with C142 Theorem II.6. There C_D(eta) = (4/pi)(log(1/eta) + pi/2) and C_D1(eta) of the same order: the logarithm came from charging each cell at first order through the full dressing kernel calR. Here the root positions are controlled at the level of counting (Lemma 3.1, an O(1) number of roots), and the density follows by one integration by parts against T' or T_A' * T_B', whose L^1 norms are of order 1/eta and 1/eta^2 with no logarithm. The factorization is what turns a counting bound into a pointwise derivative bound: a counting bound alone cannot control sigma' directly because T' changes sign.

[Proved, conditional on (R_kappa)] Corollary 3.4 (the counting boundary). |Z_0(Q_L) - nu/2| <= q(kappa0 + 1/2)/L and |z_L(Q_nu) - nu/2| <= q(kappa0 + 1/2)/L. Hence int between Q_nu and Q_L of sigma_L and of s_nu are both at most q(kappa0 + 1/2)/L (per side).

Proof. z_L = Z_0 + T*(N_x - N_0) (Lemma 2.1) and |T*(N_x - N_0)| <= q sup|N_x - N_0|. At u = Q_L, z_L = nu/2; at u = Q_nu, Z_0 = nu/2. The integrals are the increments of z_L and of Z_0 between the two points. End.

[Proved, conditional on (R_kappa)] Corollary 3.5 (comparability). For every u in [-Q_nu, Q_nu] with lambda_nu(u) >= Lambda,

    sigma_L(u) >= (1 - C_D^#/Lambda) r_(Q_nu)(u),

so sigma_L >= (3/4) r_(Q_nu) where lambda_nu >= 4C_D^# (<= 2.30 at kappa0 = 0.4), and sigma_L >= 0.71 r_(Q_nu) where lambda_nu >= 2. This is C145's comparability half of (G) and C150's hypothesis (H_I), but at Q_nu (see Section 4.3 for Q_L).

# Part 4. The strip theorem

## 4.1 The energy identity at Q_nu

Let eps_Q = (I - T_Q)^(-1)(e - h) with h = h(Q) = -a r_Q(Q)/Z_Q(Q) (C136 Section 5, so eps_Q(+-Q) = 0), at Q = Q_nu, and extend it to R by

    eps^(u) = e(u) - h + int_(-Q)^Q T(u - v) eps_Q(v) dv,

which equals eps_Q on [-Q, Q] and is C^1 on R.

[Proved] Lemma 4.1. For every configuration x of n points (in particular the true roots),

    F_L(n) - L e_bulk(nu) = sum_j eps^(x_j) - L int_(-Q_nu)^(Q_nu) eps_Q sigma_x dx.

At even half filling the same holds with Q = infinity, h = 0. No density-mismatch term appears, in contrast to C136 (4.3), which is written at Q_L.

Proof. On R, e - h = eps^ - T*(eps_Q 1_[-Q,Q]). Summing over the roots, sum_j (e - h)(x_j) = sum_j eps^(x_j) - L int_(-Q)^Q eps_Q (T*mu_x) = sum_j eps^(x_j) - L int eps_Q (sigma_x - a1). Also L int_(-Q)^Q (e - h) r_Q = L int [(I - T_Q) eps_Q] r_Q = L int eps_Q (I - T_Q) r_Q = L int eps_Q a1 (T_Q is symmetric). Subtract; the h terms cancel because sum_j 1 = n = L m(Q_nu). Finally F_L(n) = Delta L/4 + sum_j e(x_j) and e_bulk(nu) = Delta/4 + int_(-Q)^Q e r_Q at m(Q) = nu (C136 Section 5). End.

[Proved] Lemma 4.2 (the extension). |eps^'| <= 2a s_nu on R, and |eps^(u)| <= 2a |Z_0(u) - nu/2| for |u| >= Q_nu.

Proof. Inside, |eps_Q'| <= 2a r_Q (C136 (5.2), refereed; sharpened in C145 Lemma B.2). Outside, integrating by parts with eps_Q(+-Q) = 0, eps^' = e' + T*(eps_Q' 1_[-Q,Q]), and |e'| = a|a1'| <= 2a a1 (C145 Lemma B.2), so |eps^'| <= 2a(a1 + T*(r_Q 1)) = 2a s_nu. The second bound integrates the first from Q. End.

## 4.2 The sector theorem

[Proved, conditional on (R_kappa)] Theorem 4.3. For every eta in (0, pi/2], L >= 2 and 0 <= n <= L/2 at which (R_kappa) holds,

    |F_L(n) - L e_bulk(Delta, n/L)| <= K^#(eta)/L,
    K^# = (2a/9)[1 + C_D1^# A_core/(eta C_D^#)] + V_W(eta)/6 + 16 (a/eta) C_D^# A_tail + 16a/3 + 4aq(kappa0 + 1/2) + 4a q^2 (kappa0 + 1/2)^2,

with A_core = 2.71765, A_tail = 3.98873 (C145, closed forms in C147), V_W(eta) = a(4I(b) - 8) + 12aB (C146 Theorem 2). At kappa0 = 0.4, sup over eta of K^# is 136.39 (at eta = 0.289), and K^# -> 129.58 as eta -> 0 (constants_C149.json).

Proof. Lemma 4.1 splits the error into a quadrature part over the cells of the true roots and a boundary part:

    F_L - L e_bulk = sum_j [eps^(x_j) - L int_(X_j) eps^ sigma_L] + B_L,     B_L = L int_(-Q_L)^(Q_L) eps^ sigma_L - L int_(-Q_nu)^(Q_nu) eps_Q sigma_L,

with X_j the x-images of the counting cells, which tile [-Q_L, Q_L] (C142 Lemma II.1, refereed C147). The quadrature part is treated exactly as in C138 Theorem B (refereed in C145 Part C and C147 Section 3.4), with r = r_(Q_nu), K_c = {u in [-Q_nu, Q_nu] : lambda_nu(u) >= M}, M = 4 C_D^#, and the three-region Lemma 1 of C138:

1. Core. On K_c, Theorem 3.3 gives |sigma_L - r| <= r/4, so g = r/sigma_L is in [4/5, 4/3], and (D1) holds with C_D1^#. The C138 computation (with the proved (W2) of C145 Theorem B.1 and the optimized split A_core) gives Var_(K_c)(g) <= (32/(9M))[C_D^# + C_D1^# A_core/eta], and with |v| <= 2a and Var(v_(Q_nu)) <= V_W(eta) (C146 Theorem 2, valid for every Q) the core contributes (1/(8L))[2a Var(g) + (4/3)V_W] = [(2a/9)(1 + C_D1^# A_core/(eta C_D^#)) + V_W/6]/L.
2. Tail inside [-Q_nu, Q_nu]. On {lambda_nu < M}, |eps^'| <= 2a r and C138's tail integral with the proved kappa(d) gives 4a(M/(eta L)) A_tail = 16 (a/eta) C_D^# A_tail/L.
3. Straddling cells: two cells, each at most 8a/(3L) (r <= (4/3) sigma_L on their K_c part).
4. Cells beyond Q_nu (only if Q_L > Q_nu). By Lemma 4.2 and Corollary 3.4, int over Q_nu < |u| < Q_L of |eps^'| <= 2a x 2 q(kappa0 + 1/2)/L.
5. Boundary. Between Q_nu and Q_L (either order), |eps^| <= 2a q(kappa0 + 1/2)/L by Lemma 4.2 and Corollary 3.4, and the sigma_L-mass there is at most q(kappa0 + 1/2)/L per side (Corollary 3.4), so |B_L| <= 2L x 2aq(kappa0 + 1/2)/L x q(kappa0 + 1/2)/L.

Adding gives K^#/L. At even half filling Q_nu = Q_L = infinity, items 4 and 5 vanish and Var(v) = 2ab <= V_W. The numerical suprema are of explicit one-variable functions on a geometric grid of 2e4 points of [1e-9, pi/2], with the eta -> 0 limits in closed form. End.

[Proved, conditional on (R_kappa) at (L, n) and (2L, 2n)] Corollary 4.4 (the strip). By matched doubling (C138 (4.4), C136 Section 9, refereed), with the field, DM and twist reductions of C129,

    L(E_(2L) - 2E_L) <= C^# M,   C^# = 5 sup K^# + pi^2/4 = 684.4 <= 685     (kappa0 = 0.4),

for every bond of the U(1) class with -1 < J_z/J <= 0, every field and DM angle and every L >= 2, so eta_0 = pi/2 (the whole attractive side). With C142 on the repulsive side (4078, refereed C147) and zero defect at J_z/J = -1, the whole U(1) class has C = 4078, conditional on (R_kappa) only. Neither (W1'), (W3') nor (G) is used: the mismatch term, the only place where (W3) entered, is absent at Q_nu. For kappa0 = 0.3, 0.35, 0.45 the constant is 619.9, 652.1, 716.9.

## 4.3 The link to C150 (shell density relative to Q_L)

C150 Theorem 5.4 needs sigma_L >= c_* r_(Q_L) on S' = {2 <= eta L r_(Q_L) < 2C_D(eta)} for some c_* > 0, on 0 < eta <= 1/2.

[Proved, conditional on (R_kappa)] Proposition 4.5. (i) If Q_L <= Q_nu, then sigma_L >= c_* r_(Q_L) on S' with c_* = 1 - (2kappa0 + 1) eta cot(eta)/(2pi) >= 0.7135 (kappa0 = 0.4), hence C_I <= 3.181/0.7135 = 4.46 and L|D| <= 112.2 + 4.87 x 4.46 <= 134 by C150 Theorem 5.4. (ii) If Q_L > Q_nu and lambda_* := eta L r_(Q_nu)(Q_nu) >= Lambda_e := q(kappa0 + 1/2) e^(kappa_eta eta) (about 3.7), then delta := Q_L - Q_nu <= eta Lambda_e/lambda_* and r_(Q_L) <= e^(kappa_eta delta)(1 + C_2 delta/eta) s_nu on [-Q_L, Q_L], with kappa_eta = 2^(1/2) max(1/sin eta, 2) and C_2(eta) <= 7.51 on eta <= 1/2 (constants_C149.json, transfer). For lambda_* >= 60 this gives c_* >= 0.3.

Proof. (i) r_Q is nondecreasing in Q pointwise (C136 (2.3)), so on [-Q_L, Q_L], s_nu = r_(Q_nu) >= r_(Q_L). By Theorem 3.3, sigma_L >= r_(Q_L) - C_D^#/(eta L), and on S', 1/(eta L) <= r_(Q_L)/2. (ii) D = r_(Q_L) - s_nu solves (I - T_(Q_L))D = T*(s_nu 1_(Q_nu < |v| < Q_L)) on [-Q_L, Q_L], a nonnegative source bounded by s_nu(Q_nu) int over the shell of T(x - v) dv. For delta <= eta this is at most C_2(delta/eta) s_nu(Q_nu) T*1_(|v| > Q_L)(x), by monotonicity of T and the pointwise bound eta T(d) <= C_1 int_d^infinity T. The barrier (I - T_(Q_L))1 >= T*1_(|v| > Q_L) and positivity of the resolvent give D <= C_2(delta/eta) s_nu(Q_nu). Finally s_nu(Q_nu) <= s_nu(x) for |x| <= Q_nu and s_nu(Q_nu) <= e^(kappa_eta delta) s_nu(x) beyond (|s_nu'| <= kappa_eta s_nu, as for sigma in C146 Proposition 3.1). The bound on delta: if delta > eta then int_(Q_nu)^(Q_nu + eta) s_nu >= eta s_nu(Q_nu) e^(-kappa_eta eta) = (lambda_*/L) e^(-kappa_eta eta) > q(kappa0 + 1/2)/L, contradicting Corollary 3.4; the same argument gives delta/eta <= Lambda_e/lambda_*. End.

[Numerical] Q_L < Q_nu in every sampled case (check_density.py, 1 + Delta = 1e-2, 1e-4, 1e-6, L to 4096; for instance Q_L - Q_nu = -0.38 at n = L/2 - 1, and L(Z_0(Q_L) - nu/2) between -0.003 and -0.111). So case (i) is the one met in practice. A proof of the sign, which would make (i) unconditional in the sense of this part, is not given. The cleaner statement is Theorem 4.3, which avoids Q_L altogether. C150's own Lemma 5.1 (the sparse count, the dominant 97 of its 133) could also be rerun with Lemma 3.1 in place of C136's cumulative bound 3B/(2L), which removes its 3B term.

# Part 5. Status of (R_kappa)

## 5.1 What is proved about it

[Proved] (a) Free point eta = pi/2: T = 0, theta = 0, so F_j(y^+)/(2pi L) = Z_0(y^+_j) - I_j/L = kappa > 0 exactly, for every kappa0 > 0.
[Proved] (b) Exact reformulation (2.3): (R_kappa) is a family of inequalities for the explicit kernels K_j(t) = T(y_j - V(t)) V'(t) against the explicit staircase g, at the nodes of an explicit configuration.
[Proved] (c) A necessary condition: (R_kappa) implies, through Lemma 3.1 at x*, that L|Z_0(x*_j) - I_j/L| <= kappa0 for the true roots. [Numerical] The true maximal displacement is 0.249 at half filling and 1 + Delta = 1e-4 and 0.263 at 1e-6 (rk_test.py, maxdisp), close to b/4, so kappa0 cannot be taken below about b/4, and kappa0 < 1/2 is forced by existence of the top node at half filling. kappa0 = 0.4 sits in the middle of this window.
[Heuristic] (d) The bulk part of (2.4) is the first-order quadrature error of a T-average of a sawtooth, of size L eta^2 rho'/(e^(2pi lambda) - 1) (Section 2), and the edge part is matched by the edge profile of (I - T_Q)1. Turning this into a proof requires (i) a quadrature lemma for the sawtooth against T(y_j - V(t))V'(t) with the sharp constant b/(2pi) in the sparse limit (a crude first-order bound loses a factor 2, which is exactly the margin between b/4 and 1/2 at half filling), and (ii) control of rho = r_(Q_nu) near the Fermi edge on the kernel scale (r'/r of order 1/d down to d of order eta), for which only the reflection bound of C145 Theorem B.1 (distance >= 0.7) and the log-Lipschitz bound |r'| <= kappa_eta r are available.

## 5.2 Numerical verification

rk_test.py evaluates F(y^(+-)) exactly for the explicit trials (continuum from fftsolver.py: uniform grid h = eta/10 with Gregory end weights of order 7 and FFT Toeplitz products, m(Q) converged to 1e-12; Z_0 is the Nystrom counting function evaluated by FFT and cubic Hermite interpolation; half filling in closed form) and reports min_j F_j(y^+)/(2pi L) in units of eta/L (the subsolution value is its negative by symmetry, and both are printed). Fillings n in {L/32, L/8, L/4, 3L/8, L/2 - L/16, L/2 - 8, L/2 - 2, L/2 - 1, L/2}.

| 1 + Delta | L | min over n of L F_j(y^+)/(2 pi L eta) at kappa0 = 0.25 / 0.3 / 0.4 / 0.45 | worst n |
| --- | --- | --- | --- |
| 0.5 | 16 to 1024 | 0.154 / 0.187 / 0.246 / 0.275 | L/2 |
| 0.1 | 16 to 1024 | 0.102 / 0.141 / 0.220 / 0.260 | L/2 |
| 1e-2 | 64 to 4096 | 0.030 / 0.066 / 0.141 / 0.179 | L/2 |
| 1e-4 | 256 to 4096 | -0.008 / 0.025 / 0.093 / 0.129 | L/2 |
| 1e-6 | 16 to 8192 | -0.035 / 0.020 / 0.088 / 0.124 | L/2 (7L/16 for 0.25) |

(R_kappa) at kappa0 = 0.4 holds in all 173 cases, with the margin decreasing toward the asymptotic (2/pi)(0.4 - b/4) = 0.096 as eta -> 0 at half filling (0.088 at 1 + Delta = 1e-6, where b/4 = 0.2501 and the finite-L corrections are still visible). kappa0 = 0.25 fails at 1e-4 and 1e-6, as predicted by kappa0 > b/4. The rows at 0.5 and 0.1 also include L = 16.

## 5.3 Consistency checks

* Free point (Delta = 0) [Proved]: (R_kappa) is exact (5.1(a)); C_D^# = C_D1^# = 0 there since cot(pi/2) = 0 and cot(eta'/2) = cot(pi/2) = 0; Theorem 4.3 reduces to C138 Theorem B at T = 0.
* eta -> 0 at fixed L [Proved for the limit equations, Heuristic for uniformity]: dividing (2.1) by eta, both kappa(1 - q) = (2eta/pi)kappa and the sawtooth term are of order eta, so (R_kappa) tends to the same statement for the convex log-gas L tanh x_j = 2 sum coth(x_j - x_k) of C138 Section 2.3 (quantiles of its continuum equilibrium density, shifted by kappa0/L, against the log-gas gradient map). [Numerical] At 1 + Delta = 1e-6 and eta L = 0.023, 0.091, 0.36 (L = 16, 64, 256) the margin at kappa0 = 0.4 is at least 0.088 eta/L, the same as at large L.
* C142 Theorem II.6 [Proved, comparison]: Theorem 3.3 replaces C_D(eta) = (4/pi)(log(1/eta) + pi/2) (9.97 at 1 + Delta = 1e-6) by (2kappa0 + 1)/pi = 0.573 and C_D1(eta) (10.35) by 1.459, with r_(Q_nu) in place of r_(Q_L). [Numerical] On Bethe roots (check_density.py, 30 cases at 1 + Delta = 1e-2, 1e-4, 1e-6, L to 4096), eta L sup|sigma_L - s_nu| <= 0.317 and eta^2 L sup|sigma_L' - s_nu'| <= 0.217, and L sup|N_x - N_0| <= 0.758 against kappa0 + 1/2 = 0.9. The value 0.317 coincides with C138's and C142's largest eta L sup|sigma_L - r_(Q_L)|, attained on isolated sparse tail roots.
* Half filling [Proved]: Q_nu = Q_L = infinity, trials explicit, Theorem 4.3 has no boundary terms.
* C138 Section 5.3 [Numerical, explained]: the growth of the cumulative error and of the mismatch L(nu - m(Q_L)) near a Fermi boundary is an artifact of normalizing at Q_L, where Z_(Q_L) has a square-root edge; at Q_nu the cumulative error is bounded (Lemma 3.1) and the counting boundaries differ by at most q(kappa0 + 1/2)/L in counting units (Corollary 3.4).
