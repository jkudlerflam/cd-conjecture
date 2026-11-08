# C150 proof notes: the ferromagnetic strip by direct doubling

Written by a Claude Opus agent (2026-11-07), not by Codex. J = 1 (restore by homogeneity). Tags: [Proved] is a mathematical deduction (where a constant is the supremum of an explicit one-variable function, it was evaluated numerically on a grid, which is said), [Numerical] is floating point without enclosures, [Heuristic] has an unproved step, [Assumed] is a hypothesis. No independent referee pass has been made.

Summary.

- [Proved] An exact decomposition D = A + B of the matched defect D = F_2L(2n) - 2F_L(n) through a trial configuration of the 2L system built from the L roots, with the quantum numbers matched exactly (Section 1).
- [Proved] A homotopy lemma: for every trial configuration yt of the 2L system, |E_2L(y*) - E(yt)| <= 2a sum_k |rho_k(yt)|, where rho is the residual of yt in counting units and a = pi sin(eta). There is no linearization error, and the 1/eta amplification of (I - T)^(-1) is absorbed exactly by the row sums of the Jacobian (Section 2).
- [Proved] The residual of the refined grid is a second difference of theta(x - Y(t)) in the counting coordinate t of the L system. Per-cell Euler-Maclaurin then bounds A and B on any block of cells by block-boundary terms plus a remainder controlled by I_1 = int dy/(eta L sigma_L). On cells outside the block the doubled roots (x_j, x_j) are used, which have A = 0 and residual exactly +-1/(4L) (Sections 3 and 4).
- [Proved] Theorem 4.4, unconditional for every 0 < eta <= 1/2, L >= 2, 0 <= n <= L/2 and every central block: L|D| <= a(N_sp + (H+1)/4) + (a/eta)(9/512)(2 C_3 H + C_3e) I_1.
- [Proved] Corollary 5.3: L|D| <= 85.8 log(1/eta) + 180 on 0 < eta <= 1/2, unconditionally (Theorem II.6 of C142 supplies I_1 on its region, C136 (6.9) and C138's tail lemma count the rest).
- [Proved, conditional on (H_I)] Theorem 5.4: if int over the shell S' = {2 <= eta L r_Q < 2C_D(eta)} of dy/(eta L sigma_L) is at most C_I, then L|D| <= 112.2 + 4.87 C_I on 0 < eta <= 1/2. (H_I) follows from sigma_L >= c_* r_Q on S' with any c_* > 0, with C_I = 3.18/c_*, so C145's (G_c) (c_* = 3/4) gives L|D| <= 133. No derivative bound and no variation of r/sigma is needed, and there is no log log.
- [Assumed] (H_I) is not proved. It is the only gap. Numerically it holds with wide margin.

## 0. Conventions and inherited facts

[Proved, inherited] Notation of C136 Section 1 and C138 Section 1. eta = pi - gamma in (0, 1/2], Delta = -cos(eta), a = pi sin(eta), q = 1 - 2eta/pi, B = 1/(1-q) = pi/(2eta), b = pi/gamma in [1, 2). The kernels are, in closed form (checked symbolically in sym_check.py),

    T(s)  = sin(2eta) / (2 pi (sinh^2 s + sin^2 eta)) = -theta'(s)/(2pi),   theta(s) = -2 atan(tanh s / tan eta),
    a1(x) = sin(eta) / (2 pi (sinh^2 x + cos^2(eta/2))) = p'(x)/(2pi),
    e(x)  = -a a1(x).

The (L, n) roots x_1 < ... < x_n solve L p(x_j) - sum_k theta(x_j - x_k) = 2 pi I_j with I_j = j - (n+1)/2, at the minimizing flux L pi, and F_L(n) = Delta L/4 + sum_j e(x_j) (C136 Section 1, identification from C134). The counting function is z_L(x) = [p(x) - L^(-1) sum_k theta(x - x_k)]/(2pi), with z_L(x_j) = I_j/L and z_L' = sigma_L = a1 + L^(-1) sum_k T(. - x_k) > 0. Y denotes the inverse of z_L.

Inherited results used below, all proved and refereed. (C136 (6.9)) the cumulative function G of mu - r_Q^occ dx, mu = L^(-1) sum delta_(x_j), Q = Q_L, satisfies ||G||_inf <= 3B/(2L). (C138 Section 4.3, tail and core integrals, refereed in C145 Part C and C147 Section 3.4) for rho > 0, with K_c = {x in [-Q,Q]: r_Q(x) >= rho} and T_c its complement in [-Q,Q], int_(T_c) r_Q <= 2 rho (1 + 2/(kappa b)) and int_(K_c) dx/r_Q <= 2(1 + 1/(kappa b))/rho, with kappa the (W2) constant at distance one. (C145 Theorem B.1, refereed C147) r_Q >= (1 - 4e^(-2d)) r_0 at distance d from the Fermi point, so kappa = 1 - 4e^(-2) = 0.45866. (C142 Theorem II.6, refereed C147) eta L |sigma_L - r_Q| <= C_D(eta) = (4/pi)(log(1/eta) + pi/2) on [-Q_L, Q_L]. The roots lie in (-Q_L, Q_L), where z_L(Q_L) = nu/2 (C136 Section 4), and Q_L = infinity at even half filling.

## 1. The exact decomposition (Task 1)

### 1.1 Quantum numbers

[Proved] Put h = 1/L, t_j = I_j/L for the (L, n) system and t'_k = I'_k/(2L), I'_k = k - (2n+1)/2, k = 1, ..., 2n, for the (2L, 2n) system. Then

    t'_(2j-1) = t_j - h/4,     t'_(2j) = t_j + h/4,     I'_(2j-1) + I'_(2j) = 4 I_j.

Indeed t'_(2j) = (2j - n - 1/2)/(2L) = t_j + 1/(4L). The coarse cells C_j = [t_j - h/2, t_j + h/2] tile [-nu/2, nu/2] and the fine cells of width h/2 refine them two to one. Both systems are at their own minimizing flux (L pi and 2L pi), and both Bethe systems have the same form with these quantum numbers. Every 1 <= n <= L/2 is covered, including odd L and odd n (C136 Section 1 doubles to (2L, 2n) always).

### 1.2 The Yang-Yang structure of the doubling

[Proved] Let S_L(x) = sum_j [L P(x_j) - 2 pi I_j x_j] - sum_(j<k) Th(x_j - x_k), with P' = p and Th' = theta, so that grad S_L = 0 are the Bethe equations and Hess S_L = 2 pi L J_L. For the 2L system, S_2L(y) = sum_k [2L P(y_k) - 2pi I'_k y_k] - sum_(k<l) Th(y_k - y_l) and

    d S_2L / d y_k = 4 pi L F_k(y),     F_k(y) = z_2L[y](y_k) - t'_k,
    z_2L[y](u) = [p(u) - (2L)^(-1) sum_l theta(u - y_l)]/(2pi).

Its Hessian is 4 pi L J(y) with

    J_kk = a1(y_k) + (2L)^(-1) sum_(l != k) T(y_k - y_l),    J_kl = -(2L)^(-1) T(y_k - y_l)  (k != l).

J is symmetric, its off-diagonal entries are <= 0, and its row sums are a1(y_k) > 0, at every y in R^(2n), ordered or not. So J(y) is a strictly diagonally dominant nonsingular M-matrix, J(y)^(-1) >= 0 entrywise, and S_2L is strictly convex on all of R^(2n).

On doubled configurations (x, x), meaning y_(2j-1) = y_(2j) = x_j, one has S_2L(x, x) = 4 S_L(x) - n Th(0), E_2L(x, x) = 2 E_L(x), and grad S_2L(x*, x*) = (+pi, -pi, +pi, -pi, ...) at the L roots x*, since 2I_j - I'_(2j-1) = 1/2 and theta(0) = 0. So the doubled L ground state is a configuration of exactly twice the energy whose only defect is a uniform splitting force on each pair. This is the discrete content of "doubling each root".

### 1.3 Trials and the decomposition

[Proved] Define the refined grid xh_(2j-1) = Y(t_j - h/4), xh_(2j) = Y(t_j + h/4). Since z_L is the counting function of the L system, xh is exactly where a 2L system with counting function z_L would put its roots. Fix a set beta of consecutive coarse cells (the block, central and symmetric below) and let the mixed trial be

    yt_(2j-1), yt_(2j) = xh_(2j-1), xh_(2j)   (j in beta),       = x_j, x_j   (j not in beta).

Then, with y* the (2L, 2n) roots and E(y) = sum_k e(y_k),

    D = F_2L(2n) - 2F_L(n) = A + B,
    A = E(yt) - 2E_L(x) = sum_(j in beta) [e(xh_(2j-1)) + e(xh_(2j)) - 2 e(x_j)],
    B = E(y*) - E(yt).                                                       (1.1)

The terms Delta(2L)/4 and 2 Delta L/4 cancel. A is the quadrature part. In the coordinate t, with f = e o Y, A = sum_(j in beta) Delta2_j f where Delta2_j phi = phi(t_j + h/4) + phi(t_j - h/4) - 2 phi(t_j), that is, A = 2L (M_(h/2) - M_h) f restricted to the block, the Richardson difference of the two midpoint rules for the same function. B is the self-consistency part.

[Proved] The residual of the trial is

    rho_k = F_k(yt) = zt(yt_k) + s_k,    zt(u) = z_2L[yt](u) - z_L(u),
    s_k = 0 for members of block pairs, s_(2j-1) = +1/(4L), s_(2j) = -1/(4L) for doubled pairs.   (1.2)

For block members this is z_L(xh_k) = t'_k; for doubled members it is t_j - t'_(2j-1) = 1/(4L) and t_j - t'_(2j) = -1/(4L).

## 2. The self-consistency term (Task 2)

[Proved] Theorem 2.1. For every trial yt in R^(2n),

    |E(y*) - E(yt)| <= 2a sum_k |F_k(yt)|.

Proof. For s in [0, 1] let G_s(y) = S_2L(y) - 4 pi L (1-s) rho . y, rho = F(yt). G_0 is strictly convex with critical point yt, G_1 = S_2L is strictly convex with critical point y*. A convex function attaining its minimum on a bounded set has bounded sublevel sets, so G_0 and G_1 are coercive, hence so is G_s = (1-s) G_0 + s G_1. Thus G_s has a unique minimizer y(s), which solves F(y(s)) = (1-s) rho, with y(0) = yt and y(1) = y*. By the implicit function theorem (Hessian 4 pi L J positive definite) y(s) is C^1 and J(y(s)) y'(s) = -rho. Hence

    E(y*) - E(yt) = - int_0^1 w(s) . rho ds,    w(s) = J(y(s))^(-1) e'(y(s)),       (2.1)

using the symmetry of J. Now |e'| = a |a1'| and a1'/a1 = -2 sinh(2x)/(cosh 2x + cos eta), so |e'| <= 2a a1. The constant vector 2a.1 has J(2a.1) = 2a a1 >= |e'| componentwise, and J^(-1) >= 0, so |w_k| <= (J^(-1)|e'|)_k <= 2a for every k and s. Insert in (2.1). End of proof.

[Proved, interpretation] This answers the question of Task 2 in a sharp form. The Neumann inverse (I - T)^(-1) has norm B = pi/(2eta), and its discrete version J^(-1) has the same size on smooth data. But the energy only ever meets J^(-1) through J^(-1) e', and e' is dominated by a multiple of the row-sum vector J 1 = a1. Thus J^(-1) e' is of order a, which is the discrete statement |eps_Q'| <= 2a r_Q, v_Q <= 2a, and the 1/eta amplification is absorbed exactly, configuration by configuration, with no knowledge of the roots. The source rho is then charged at first order, but rho is itself a second difference (Section 3), so its size is that of a second-order quadrature error. [Numerical] max|w|/(2a) is 0.51 to 0.52 at eps <= 1e-2 (0.58 at eps = 0.1) both at the refined grid and at y* (the dressed velocity is about ab, half of 2a), and 0.995 at the free point, where the bound is attained in the tails.

## 3. The residual as a second difference

[Proved] Lemma 3.1. For every u,

    zt(u) = sum_(j in beta) zeta_j(u),
    zeta_j(u) = (2L)^(-1) [ int_(xh_(2j-1))^(x_j) - int_(x_j)^(xh_(2j)) ] T(u - y) dy
              = -(4 pi L)^(-1) Delta2_j g_u,     g_u(t) = theta(u - Y(t)).

Proof. The p terms cancel in zt. With Th~ = -theta/(2pi), increasing with Th~' = T, zt(u) = int Th~(u - y) d(mu~ - mu)(y), where mu~ = (2L)^(-1) sum delta_(yt_k). Doubled pairs cancel in mu~ - mu. Its cumulative function N vanishes outside the block pairs and equals +1/(2L) on (xh_(2j-1), x_j) and -1/(2L) on (x_j, xh_(2j)). Integrating by parts (N has compact support and zero total mass), zt(u) = int T(u - y) N(y) dy. In the t coordinate, T(u - Y(t)) Y'(t) = (2pi)^(-1) d/dt theta(u - Y(t)), which gives the second form. End of proof.

[Proved] Lemma 3.2 (pair width). For every root, 0 < xh_(2j) - x_j <= delta and 0 < x_j - xh_(2j-1) <= delta, with tanh(delta) = tan(eta) (eta < pi/4). Proof: 1/(4L) = int_(x_j)^(xh_(2j)) sigma_L >= L^(-1) int_0^u T = (pi L)^(-1) atan(tanh u / tan eta), so tanh u <= tan eta; the same on the left. No information on the other roots is used. [Numerical] delta/eta <= 1.007 for eta <= 0.1 and 1.226 at eta = 1/2.

[Proved] Lemma 3.3 (comparability of the trial density). For every y, (2L)^(-1) sum_k T(yt_k - y) <= H (sigma_L(y) - a1(y)) <= H sigma_L(y), with

    H(eta) = (1/2)[1 + sup_(s >= 0) (sinh^2 s + sin^2 eta)/(sinh^2 (s - delta)_+ + sin^2 eta)].

Proof: per pair, (1/2)[T(y - x_j + v) + T(y - x_j - u)] with 0 <= u, v <= delta; for s = y - x_j >= 0 the first is at most T(s) and the second at most T((s - delta)_+) since T is even and decreasing in |s|; symmetric for s < 0. Doubled pairs give exactly T(y - x_j)/L. [Proved, sup evaluated numerically] H = 1.809 as eta -> 0 (the Lorentzian value (1 + (3+sqrt5)/2)/2), 1.827 at eta = 0.1, 2.582 at eta = 1/2, increasing in eta on the grid.

[Proved] Lemma 3.4 (local derivative bounds). Put tau_1 = 2eta/sin(2eta), tau_2 = eta^2 sup_s |T''/T| (equal to 2eta^2/sin^2 eta, attained at s = 0), alpha_1 = 2, alpha_2 = 4, alpha_3 = 8, s_i = max(tau_i, eta^i alpha_i). Then |T^(i)| <= tau_i eta^(-i) T, |a1^(i)| <= alpha_i a1, |e^(i)| <= a alpha_i a1, and, since sigma_L is a positive combination of a1 and translates of T, |sigma_L^(i)| <= s_i eta^(-i) sigma_L, for i = 1, 2 (alpha up to i = 3). Proof: T'/T = -2 sinh 2s/(cosh 2s - cos 2eta) has supremum 2/sin(2eta) at cosh 2s = 1/cos 2eta; the other ratios are explicit rational functions of cosh 2s, sinh 2s (constants_C150.py), whose suprema were evaluated on a log grid from 1e-9 eta to 40 together with their limits at infinity (4 for T''/T and a1''/a1, 8 for a1'''/a1). These bounds hold for every root configuration. They are local Harnack bounds on the scale eta and contain no information on the density profile.

[Proved] Lemma 3.5 (derivatives in the counting coordinate). With D = sigma^(-1) d/dy = d/dt, sympy gives

    D^2 g_u / (2pi) = T''/sigma^3 + 3 T' sigma'/sigma^4 - T sigma''/sigma^4 + 3 T sigma'^2/sigma^5    (T^(i) at u - y),
    D^3 f = e'''/sigma^3 - 3 e'' sigma'/sigma^4 - e' sigma''/sigma^4 + 3 e' sigma'^2/sigma^5,

(the first line is d^3 g_u/dt^3 divided by 2pi), so by Lemma 3.4

    |d^3 g_u/dt^3| <= 2pi C_3 T(u - y)/(eta^2 sigma(y)^3),     C_3 = tau_2 + 3 tau_1 s_1 + s_2 + 3 s_1^2,
    |d^3 f/dt^3|   <= C_3e a a1(y)/(eta^2 sigma(y)^3),          C_3e = eta^2 alpha_3 + 3 eta alpha_2 s_1 + alpha_1 (s_2 + 3 s_1^2).

C_3 = 10.0 and C_3e = 10.0 as eta -> 0, and C_3 = 12.83, C_3e = 21.96 at eta = 1/2 (constants_C150.json, increasing in eta).

## 4. Euler-Maclaurin on the block (Task 3) and the master bound

[Proved] Lemma 4.1. Let w(tau) = (1/4 - |tau|)_+, v(tau) = int_(-1/2)^tau (w - 1/16), so that v(+-1/2) = 0, v is odd, and ||v||_inf = 9/512 (attained at tau = -+3/16). For phi in C^3 on a union beta of consecutive cells j_1, ..., j_2,

    sum_(j in beta) Delta2_j phi = (h/16)[phi'(t_(j2) + h/2) - phi'(t_(j1) - h/2)] - sum_(j in beta) int_(C_j) h^2 v((t - t_j)/h) phi'''(t) dt.

Proof. Delta2_j phi = int_(C_j) h w((t - t_j)/h) phi''(t) dt (Peano kernel of the second difference), the mean of h w over C_j is h/16, and one integration by parts against h^2 v, which vanishes at the cell edges. The first terms telescope. End of proof.

[Proved, interpretation] This is the Richardson structure of A. The leading term is (h/16)[phi'] = 2L(h^2/32)[phi'], which is 2L times the difference -(h^2/96) + (h^2/24) of the midpoint Euler-Maclaurin coefficients at steps h/2 and h; the O(1) terms cancel identically, the edge term is the only second-order term, and the remainder is third order. For the free fermion it reproduces the Casimir energy (Section 6).

[Proved] Lemma 4.2 (the self-consistency source). For any block beta of consecutive cells, with I_1 = int over the block image of dy/lambda_L(y), lambda_L = eta L sigma_L,

    sum_k |zt(yt_k)| <= H/(8L) + (9/512) C_3 H I_1/(eta L).

Proof. By Lemma 3.1 and Lemma 4.1 applied to phi = g_u, |zt(u)| <= (4pi L)^(-1)[(h/16)(|g_u'(t_b1)| + |g_u'(t_b2)|) + (9/512) h^2 int_beta |g_u'''| dt]. Here g_u'(t) = 2pi T(u - Y(t))/sigma(Y(t)), and int |g_u'''| dt = int |g_u'''| sigma dy <= 2pi C_3 int T(u - y) eta^(-2) sigma^(-2) dy. Sum over u = yt_k and use Lemma 3.3, sum_k T(yt_k - y) <= 2LH sigma_L(y). Each boundary term gives at most 2LH sigma(y_b)/(32 L^2 sigma(y_b)) = H/(16L). The remainder gives (4pi L)^(-1)(9/512) L^(-2) 2pi C_3 2LH eta^(-2) int sigma^(-1) dy, and int sigma^(-1) dy = eta L I_1. End of proof.

[Proved] Lemma 4.3 (the quadrature term). L|A| <= a/4 + (a/eta)(9/512) C_3e I_1.

Proof. Lemma 4.1 with phi = f = e o Y. The boundary terms are (h/16)|e'/sigma| <= (h/16) 2a a1/sigma <= a/(8L) each, since a1 <= sigma_L. The remainder is at most (9/512) h^2 int |f'''| sigma dy <= (9/512) L^(-2) C_3e a eta^(-2) int a1 sigma^(-2) dy <= (9/512) L^(-2) C_3e a eta^(-2) int sigma^(-1) dy. End of proof.

[Proved] Theorem 4.4 (master bound). For 0 < eta <= 1/2, every L >= 2, 0 <= n <= L/2 and every block beta of consecutive cells, with N_sp = n - |beta| the number of doubled cells,

    L |D| <= a [N_sp + (H + 1)/4] + (a/eta)(9/512)(2 C_3 H + C_3e) I_1,          (4.1)

and the term (H+1)/4 and I_1 are absent if beta is empty.

Proof. D = A + B by (1.1). By Theorem 2.1 and (1.2), L|B| <= 2aL [sum over doubled members of 1/(4L) + sum_k |zt(yt_k)|] = a N_sp + 2aL sum_k |zt(yt_k)|. Insert Lemma 4.2 and add Lemma 4.3. End of proof.

[Proved] Corollaries of (4.1) that use nothing else. (i) beta empty: |F_2L(2n) - 2F_L(n)| <= a n/L = pi sin(eta) nu, for every Delta in (-1, 0] (the restriction eta <= 1/2 is not needed when beta is empty). This is a Theorem-A-type bound for the doubling, with constant pi/2 in place of the 15pi/4 of C138 Theorem A. (ii) So L|D| <= (pi/2) eta L, which is fine for eta L bounded. The content of the rest is to choose beta so that N_sp = O(1/eta) and I_1 = O(1).

## 5. Choice of the block, counting, and the strip

[Proved] Lemma 5.1 (counting). Fix Lambda > 0 and let X_Lambda in [0, Q] be the point with eta L r_Q(X_Lambda) = Lambda (X = 0 if eta L r_Q(0) <= Lambda, X = Q if eta L r_Q(Q) >= Lambda; r_Q is even and nonincreasing on [0, Q]). Let beta be the cells whose image lies in [-X_Lambda, X_Lambda], a central set of consecutive cells. Then

    N_sp <= (2 Lambda/eta)(1 + 2/(kappa b)) + 3pi/eta + 2.

Proof. The cells not in beta are those of roots with |x_j| > X_Lambda plus at most one straddling cell per side. The number of roots in (X, Q_L) is L[G(Q_L) - G(X)] + L int_X^Q r_Q <= 3B + L int_X^Q r_Q by C136 (6.9), and likewise on the left. By C138's tail lemma with rho = Lambda/(eta L), L int_(|x| > X) r_Q <= 2 (Lambda/eta)(1 + 2/(kappa b)). End of proof.

[Proved] Lemma 5.2. If lambda_L >= c_0 lambda_Q on [-X_Lambda, X_Lambda], then I_1 <= 2(1 + 1/(kappa b))/(c_0 Lambda). Proof: C138's core integral with rho = Lambda/(eta L).

[Proved] Corollary 5.3 (unconditional, with a logarithm). For 0 < eta <= 1/2, L >= 2, 0 <= n <= L/2,

    L |F_2L(2n) - 2F_L(n)| <= 85.8 log(1/eta) + 180.

Proof. Take Lambda = 2 C_D(eta). On [-X_Lambda, X_Lambda], Theorem II.6 gives lambda_L >= lambda_Q - C_D >= lambda_Q/2, so c_0 = 1/2. Lemma 5.1 gives a N_sp <= pi[21.442 C_D + 9.425] + 2pi sin(eta), using a/eta <= pi and 2(1 + 2/(kappa b)) <= 10.721. Lemma 5.2 gives I_1 <= 6.361/C_D <= 2.21. With the constants at eta = 1/2 (all increasing in eta), a(H+1)/4 <= 1.35 and (a/eta)(9/512)(2C_3 H + C_3e) <= 4.869. Insert C_D = (4/pi) log(1/eta) + 2. The sum is 85.77 log(1/eta) + 134.72 + 29.61 + 3.01 + 1.35 + 10.75 <= 85.8 log(1/eta) + 180. End of proof.

This is length-uniform at each eta and replaces the order 5 x 247 log(1/eta) that C142's sector bound gives for the defect. It is not the target.

[Assumed] Hypothesis (H_I). With S' = {x in [-Q_L, Q_L]: 2 <= eta L r_Q(x) < 2 C_D(eta)},

    int_(S') dy / (eta L sigma_L(y)) <= C_I,

uniformly in eta in (0, eta0], L and n.

[Proved, conditional on (H_I)] Theorem 5.4. For 0 < eta <= 1/2, L >= 2, 0 <= n <= L/2,

    L |F_2L(2n) - 2F_L(n)| <= 112.2 + 4.87 C_I.

Proof. Take Lambda = 2, so beta is the set of cells inside {lambda_Q >= 2} = S' together with the deep core {lambda_Q >= 2C_D}. Lemma 5.1 gives a N_sp <= pi(4(1 + 2/(kappa b)) + 3pi) + 2pi sin(eta) <= 96.97 + 3.01. On the deep core Theorem II.6 gives lambda_L >= lambda_Q/2, so its part of I_1 is at most 2 x 2(1 + 1/(kappa b))/(2 C_D) <= 2.21; on S' it is at most C_I by (H_I). With a(H+1)/4 <= 1.35 and the coefficient 4.869 of I_1, L|D| <= 99.98 + 1.35 + 4.869 (2.21 + C_I) <= 112.2 + 4.87 C_I. End of proof.

[Proved] Sufficient conditions for (H_I). (a) If sigma_L >= c_* r_Q on S' for some c_* > 0, then C_I <= (1/c_*) int_(S') dy/lambda_Q <= 3.181/c_* by C138's core integral with rho = 2/(eta L). (b) C145's (G_c), sigma_L >= (3/4) r_Q on S = {2 <= lambda_Q < 4C_D} (S' is contained in S), gives C_I <= 4.24 and L|D| <= 133. (c) Any pointwise comparability constant works; even c_* = 1/10 gives L|D| <= 267. (H_I) needs no derivative information and no variation bound, in contrast to C145's (G), and no log log appears, in contrast to C146 Proposition 3.1, because the remainder weight 1/lambda is integrable through the shell (int_2^infinity d lambda/lambda^2 after the change of variables) rather than charged per unit length.

[Proved, bookkeeping] Composition of the constants at eta0 = 1/2. Sparse cells 96.97 (of which 67.36 from the continuum tail count and 29.61 from the cumulative error 3B of C136), straddling cells 3.01, block-boundary terms 1.35, deep core 10.75, shell 4.87 C_I. The sparse count dominates and is crude by about a factor 10 (Section 7).

## 6. Consistency checks

[Proved] Free point Delta = 0 (eta = pi/2). T = 0, so z_2L[y] = z_L = p/(2pi) for every y, the refined grid is the (2L, 2n) solution, rho = 0 and B = 0 exactly. D = A = 2L(M_(h/2) - M_h) f with f = e o Y. Here p(x) = 2 atan(tanh x), e(x) = -1/cosh(2x), and with k = p one has cos k = 1/cosh(2x), so f(t) = -cos(2 pi t) is analytic in t, including at the half-filled edges t = +-1/4 where Y runs to infinity. Lemma 4.1 with the whole range as one block gives L D = (1/16)[f'(nu/2) - f'(-nu/2)] + O(1/L^2) = (pi/4) sin(pi nu) + O(1/L^2), which is the CFT value (pi/4) c v with c = 1 and the free Fermi velocity v = sin(pi nu); at half filling it is pi/4. [Numerical] L D = 0.7898 at L = 16, against pi/4 = 0.7854.

[Numerical] CFT value at half filling. L D against pi v/4 with v = (pi/2) sin(gamma)/gamma: see RESULT.md, ratio 0.9999 at eps = 1e-2, L = 1024.

[Proved and Numerical] The limit eta -> 0 at fixed L. Both systems tend to the log-gas of C138 (2.3), whose sum rule (C138 (2.4)) gives D = eps 2n(L-n)/((2L-1)(L-1)) + O(eta^4). [Numerical] At eps = 1e-6, (L, n) = (64, 32): L D = 1.63765e-5 against 1.63820e-5, and (16, 5): 3.78488e-6 against 3.78495e-6. The bound of Corollary 4.4(i) gives L D <= pi eta n here, consistent.

[Numerical] Where the Casimir energy sits. At small eta it is in B, not in A: at eps = 1e-2, L = 1024, half filling, L A_ref = 0.0002 and L B_ref = 0.0578 for the refined trial (L D = 0.0580). The bare energy is too smooth to carry it (e = O(eta^2), smooth on scale one), the dressing carries it. At Delta = 0 it is all in A.

## 7. What is not proved, and why the shell remains

[Heuristic] The two inputs N_sp and I_1 of (4.1) are first-order, integrated quantities. N_sp is controlled by the continuum because the cumulative error of the roots is O(B/L) = O(1/(eta L)) (C136), which is exactly the accuracy needed to count O(1/eta) tail roots. I_1 asks that the true local density lambda_L does not stay of order one over an x-length of order one inside the region where the continuum density is moderately large. Pointwise, this needs sigma_L >= c r_Q on the shell; integrated, it needs the sum of squared root gaps over the shell to be O(eta). Theorem II.6 controls the density only to C_D(eta)/(eta L), which is the logarithm of Corollary 5.3. The logarithm is produced by the 1/|y| tail of T on scales between eta and 1 (the symbol 1 - T^(w) ~ eta(|w| + 2/pi) is a first-order operator whose Green function is logarithmic), so any argument that charges a discreteness source at first order through (I - T)^(-1) loses it. The doubling route removes every other loss: no comparison of sigma_L with r_Q at second order is needed, and the 1/eta amplification is absorbed exactly (Section 2).

[Proved, obstruction for this route] The comparison of 2L with L alone cannot supply (H_I). Every bound in Sections 2 to 4 holds for an arbitrary configuration x of n points with the counting function z_L (positivity of T and Lemma 3.2 are all that is used), so the shell information must come from the Bethe equations of the L system itself, through a statement about the density of its roots.

[Numerical] Margins (RESULT.md). The pointwise ratio lambda_L/lambda_Q on S' is at least 0.9993 in every sampled case (35,000 sample points at eps = 0.1 to 1e-4), against the c_* > 0 needed, and I_1 over the whole block is at most 0.97 against the 4.24 + 2.21 used. The sampled eta N_sp is at most 4.58 against the proved 30.87. The true L|D| is of order eta, at most pi v/4, while (4.1) evaluated with the actual N_sp and I_1 is at most 16.6, and the worst-case bound is 133 under (G_c).

## 8. Files

doubling.py (trials, residuals, A and B, homotopy weights, hypothesis quantities, continuum block), constants_C150.py and constants_C150.json (Lemmas 3.2 to 3.5), sym_check.py (closed forms and derivative recursions), doubling_A.json, doubling_B.json (scans), bethe_ferro.py and continuum.py (copied from work_C145).
