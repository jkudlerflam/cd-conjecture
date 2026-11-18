# Task G: (1) adversarial referee of a derivation, (2) an independent proof of one cited input

Part (1). You are the referee of the text at the end (inputs (I1) to (I4) and Corollary 4.13 with its proof, Steps 1 to 8). Find errors: signs, a term that is dropped, a case that is not covered (Q_L < Q versus Q_L > Q, K_c empty, a cell that contains all of K_c, half filling), a misuse of the midpoint-rule Peano bound, a constant that is not what the arithmetic gives, an unstated dependence on L, n, Q or eta. Re-derive every identity and inequality yourself, including (I2), which is claimed as proved in the text. For each step say "correct" or give the precise defect and the repair. Tags such as [rigorous], [cited], (Fable) are bookkeeping. End with a verdict.

Part (2). Input (I4), Var_{[-Q,Q]}(v) <= V_W with a universal V_W, comes from a run that has not been refereed. Give an independent, complete proof that sup over eta in (0, pi/4] and Q in (0, inf) of Var_{[-Q,Q]}(v_Q) is finite (any explicit or non-explicit universal constant), or, failing that, check the proof route sketched below and say exactly which step you can or cannot confirm. The sketched route (from the unrefereed run): with R = (I - T_Q)^(-1), u = eps' = -a R a1', v = u/r, G = R[T(. - Q)], H = R[T(. + Q)] (nonnegative),
   (i) for any f = src + T_Q f on [-Q, Q], (I - T_Q) f' = src' - f(Q) T(. - Q) + f(-Q) T(. + Q);
   (ii) hence v' = F_1 - (r(Q)/r)[(v(Q) - v) G + (v(Q) + v) H] with F_1 = a[(R a1')^2 - (R a1)(R a1'')]/(R a1)^2;
   (iii) H(x) = G(-x), G >= H on [0, Q], int_{-Q}^{Q} G = Z_Q(Q) - 1, where Z_Q = R 1 <= pi/(2 eta);
   (iv) 0 <= v <= 2a tanh(2Q) on [0, Q] (from |a1'| <= 2 tanh(2|x|) a1), v odd;
   (v) (F_1)_- <= a min(4, b^2 E(4 + E)/(1 - E)^2) with b = pi/gamma, E = 4 e^(-2 b d), d = Q - x (second entry when E < 1), using r_0 - r_Q <= [r_0(2Q - x) + r_0(2Q + x)] and the same for derivatives, r_0 = 1/(2 gamma cosh(b x));
   (vi) Var_{[-Q,Q]} v = 2[v(Q) + 2 int_0^Q (v')_-] <= 4a tanh(2Q) + 4a I(b) + 12 a tanh(2Q)(Z_Q(Q) - 1), and a (Z_Q(Q) - 1) <= pi^2/2.
You are free to use a different route.

## Definitions and facts to assume

eta in (0, 1/64] for the Corollary (bounds on constants are taken over eta <= pi/4), gamma = pi - eta, b = pi/gamma, q = 1 - 2 eta/pi, a = pi sin(eta).

    p(x) = 2 atan(tan(eta/2) tanh x);  theta(x) = -2 atan(tanh x/tan eta);  a1 = p'/(2 pi) = sin(eta)/(pi (cosh 2x + cos eta)) > 0;
    T = -theta'/(2 pi) = sin(2 eta)/(2 pi (sinh^2 x + sin^2 eta)) > 0, even, decreasing in |x|, int T = q, T(0) = cot(eta)/pi;
    T = T_A * T_B with T_A, T_B positive, even, unimodal kernels, T_A(0) <= 2/(pi eta), T_B(0) <= 2/(pi eta).
    (T_Q f)(x) = int_{-Q}^{Q} T(x - v) f(v) dv on [-Q, Q]; (I - T_Q)^(-1) = sum_k T_Q^k has a nonnegative kernel and is self-adjoint on L^2(-Q, Q).
    r_Q = (I - T_Q)^(-1) a1: positive, even, nonincreasing on [0, Q], r_Q <= r_0 = 1/(2 gamma cosh(b x)); m(Q) = int_{-Q}^{Q} r_Q.
Bethe roots: integers L >= 2, 1 <= n <= L/2, nu = n/L, I_j = j - (n+1)/2; x* the unique real solution of L p(x_j) - sum_k theta(x_j - x_k) = 2 pi I_j (increasing, symmetric). N_x(v) = (1/L) #{j : x*_j <= v}, mu_x = (1/L) sum_j delta_{x*_j}, sigma_L = a1 + T*mu_x.
Continuum reference: Q = Q_nu with m(Q) = nu (Q = inf at n = L/2), s = a1 + T*(r_Q 1_[-Q,Q]) on R, N_0 the cumulative function of r_Q 1_[-Q,Q], Z_0 = p/(2 pi) + T*N_0 - q nu/2 (so Z_0' = s, Z_0 = N_0 - nu/2 on [-Q, Q]).
Proved earlier, assume: (Theorem 4.11 and Theorem H) |N_x(v) - N_0(v)| <= K/L for all real v, with a universal K >= 1. (Corollary 4.12) eta L |sigma_L - s| <= C_D := 2K eta cot(eta)/pi and eta^2 L |sigma_L' - s'| <= C_D1 := 4K eta^2 T_A(0) T_B(0), on R.

## The text to referee

Inputs for the energy, all continuum or cited statements. e(x) = -a a1(x) with a = pi sin(eta), so e = -sin^2(eta)/(cosh 2x + cos eta) and |e'| <= 2a a1.
 (I1) [cited, C136 Sections 1 and 5, refereed] The sector ground energy is F_L(n) = Delta L/4 + sum_j e(x*_j), and e_bulk(Delta, nu) = Delta/4 + int_{-Q}^{Q} e r_Q at m(Q) = nu.
 (I2) [rigorous, checked] With Z_Q = (I - T_Q)^(-1) 1 >= 1 and h = h(Q) = -a r_Q(Q)/Z_Q(Q), the function eps = eps_Q := (I - T_Q)^(-1)(e - h) = -a r_Q - h Z_Q is even, C^1, vanishes at +-Q, and satisfies |eps'| <= 2a r_Q on [-Q, Q]: differentiating eps = e - h + T_Q eps and integrating by parts (the boundary terms vanish because eps(+-Q) = 0) gives eps' = (I - T_Q)^(-1) e', and |e'| <= 2a a1 with the positivity of the resolvent gives the bound. So v := eps'/r_Q has |v| <= 2a.
 (I3) [cited, C145 Theorem B.1, refereed in C147] r_Q(x) >= kappa(d) r_0(x) with kappa(d) = 1 - 4 e^(-2d) at distance d = Q - |x| from the Fermi point; kappa := kappa(1) = 0.459.
 (I4) [cited, C146 Theorem 2, not refereed; the managing agent read its proof (C146 Lemmas 2.1 to 2.3) and found no defect, without re-deriving its numerical constant] Var_{[-Q,Q]}(v) <= V_W <= 92.9 for every eta in (0, pi/2] and Q in (0, inf].

[rigorous given (I1), (I3), (I4) and Theorems H and 4.11] (Fable, following C138 Theorem B, refereed in C145 and C147, and C149 Theorem 4.3, not refereed; written out and checked) Corollary 4.13 (sector theorem on the strip). For eta in (0, 1/64], every L >= 2 and 0 <= n <= L/2,

    |F_L(n) - L e_bulk(Delta, n/L)| <= K'/L,      K' = 35 + 181 K + 9 K^2,

with K the constant of Theorem 4.11 (K >= 1).

Proof. n = 0 is exact. Let Q = Q_nu, r = r_Q, s its extension, and extend eps to R by epsh(u) = e(u) - h + int_{-Q}^{Q} T(u - v) eps(v) dv, which equals eps on [-Q, Q] (at half filling Q = inf, h = 0, eps = -a r_0, and the steps about Q_L and the boundary below are empty).

Step 1 (energy identity at Q_nu; C149 Lemma 4.1). On R, e - h = epsh - T*(eps 1_[-Q,Q]). Summing over the roots and using T*mu_x = sigma_L - a1, sum_j (e - h)(x*_j) = sum_j epsh(x*_j) - L int_{-Q}^{Q} eps (sigma_L - a1). Also L int_{-Q}^{Q} (e - h) r = L int [(I - T_Q) eps] r = L int eps a1, by the symmetry of T_Q. Subtracting, the h terms cancel because n = L m(Q), and with (I1)

    F_L(n) - L e_bulk(nu) = sum_j epsh(x*_j) - L int_{-Q}^{Q} eps sigma_L dx.

Step 2 (the extension). epsh' = e' + T*(eps' 1_[-Q,Q]) (integration by parts, eps(+-Q) = 0), so |epsh'| <= 2a (a1 + T*(r 1)) = 2a s on R. With Z_0 the continuum counting function (Z_0' = s, Z_0(+-Q) = +-nu/2), |epsh(u)| <= 2a |Z_0(u) - nu/2| for u >= Q, and |eps(u)| <= 2a (nu/2 - Z_0(u)) for 0 <= u <= Q; similarly on the left.

Step 3 (cells). z_L(u) = p(u)/(2 pi) + (T*N_x)(u) - q nu/2 is the counting function of the roots: z_L' = sigma_L > 0, z_L(x*_j) = I_j/L. Let Q_L be the point with z_L(Q_L) = nu/2. The cells X_j = z_L^(-1)([I_j/L - 1/(2L), I_j/L + 1/(2L)]) tile [-Q_L, Q_L], each has sigma_L-mass 1/L, and x*_j is its midpoint in the variable t = z_L. Since z_L - Z_0 = T*(N_x - N_0) and |N_x - N_0| <= K/L,

    |z_L - Z_0| <= q K/L on R,   hence   |Z_0(Q_L) - nu/2| <= q K/L   and   |z_L(Q) - nu/2| <= q K/L.                 (4.4)

Write F_L - L e_bulk = sum_j [epsh(x*_j) - L int_{X_j} epsh sigma_L] + B_L with B_L = L int_{-Q_L}^{Q_L} epsh sigma_L - L int_{-Q}^{Q} eps sigma_L.

Step 4 (boundary term). B_L is an integral over the two intervals between Q and Q_L (in either order) and their mirrors. There |epsh| <= 2a |Z_0(Q_L) - nu/2| <= 2a q K/L by Step 2 and the monotonicity of Z_0, and the sigma_L-mass of each interval is |z_L(Q_L) - z_L(Q)| <= q K/L by (4.4). So |B_L| <= 2 L (2a q K/L)(q K/L) = 4 a q^2 K^2/L.

Step 5 (two bounds per cell). With f(t) = epsh(z_L^(-1)(t)), the cell term is f(t_j) - L int_cell f dt over a cell of length 1/L centred at t_j. It is at most the oscillation of f on the cell, hence at most int_{X_j} |epsh'| dx. It is also at most Var_cell(f')/(8L) (Peano kernel of the midpoint rule), and f' = epsh'/sigma_L in the variable x.

Step 6 (core). Let C_D = 2K eta cot(eta)/pi and C_D1 = 4K eta^2 T_A(0) T_B(0) (Corollary 4.12), M = 4 C_D, rho_0 = M/(eta L), and K_c = {x in [-Q, Q] : r(x) >= rho_0} = [-x_c, x_c] (r is even and nonincreasing on [0, Q]). On K_c, |sigma_L - r| <= C_D/(eta L) <= r/4, so g := r/sigma_L lies in [4/5, 4/3]. Core cells are those with X_j inside K_c. On them epsh'/sigma_L = v g, so the sum of their second bounds is at most (1/(8L)) Var_{K_c}(v g) <= (1/(8L)) [2a Var_{K_c}(g) + (4/3) V_W], by (I2) and (I4). Next g' = [r'(sigma_L - r) - r(sigma_L' - r')]/sigma_L^2 and sigma_L^2 >= (9/16) r^2, so

    Var_{K_c}(g) <= (16/9) [ (C_D/(eta L)) int_{K_c} |r'|/r^2 + (C_D1/(eta^2 L)) int_{K_c} dx/r ].

int_{K_c} |r'|/r^2 = 2 [1/r(x_c) - 1/r(0)] <= 2/rho_0. For int_{K_c} dx/r, split [0, x_c] at Q - 1: the part with Q - x < 1 has length at most 1 and r >= rho_0 there; on the rest r >= kappa r_0 by (I3), int_0^y dx/r_0 = (2 gamma/b) sinh(b y) <= 1/(b r_0(y)), and r_0(x_c) >= r(x_c) >= rho_0, so that part is at most 1/(kappa b rho_0). Hence int_{K_c} dx/r <= 2 (1 + 1/(kappa b))/rho_0 and

    (2a/8) Var_{K_c}(g) <= (2a/9) [1 + C_D1 (1 + 1/(kappa b))/(eta C_D)],      a C_D1/(eta C_D) = 2 pi^2 sin(eta) tan(eta) T_A(0) T_B(0) <= 8 tan(eta)/eta.

Step 7 (the other cells). Their first bounds sum to at most int_U |epsh'|, U the union of the non-core cells. U is contained in the tail {x_c < |x| <= Q}, the K_c-parts of at most two straddling cells, and the intervals between Q and Q_L when Q_L > Q. On a straddling cell's K_c-part, r <= (4/3) sigma_L and the cell has sigma_L-mass 1/L, so int |eps'| <= 8a/(3L); two cells give 16a/(3L). On the tail r < rho_0: per side, the part within distance 1 of Q contributes at most rho_0, and on the rest r_0 < rho_0/kappa by (I3) and int_y^inf r_0 <= e^(-b y)/(gamma b) <= 2 r_0(y)/b, so int_tail r <= 2 rho_0 (1 + 2/(kappa b)) and int_tail |eps'| <= 4a rho_0 (1 + 2/(kappa b)). Beyond Q, int |epsh'| <= 2a int s = 4a |Z_0(Q_L) - nu/2| <= 4 a q K/L. If K_c is empty every cell is treated in this way and the same tail bound applies.

Step 8 (sum). Multiplying by L,

term by term, for eta <= pi/4 (tan(eta)/eta <= 4/pi, a <= pi sin(pi/4) = 2.23, b >= 1, 1/kappa = 2.18, q <= 1):
 core, g-part: (2a/9) + (2/9)(8 tan(eta)/eta)(1 + 1/kappa) <= 0.5 + 7.2;
 core, v-part: V_W/6 <= 15.5;
 tail: 4a M (1 + 2/(kappa b))/eta = 32 K cos(eta)(1 + 2/kappa) <= 171.5 K;
 straddling cells: 16a/3 <= 11.9;
 beyond Q and boundary: 4 a q K + 4 a q^2 K^2 <= 8.9 (K + K^2).
The sum is at most 35.1 + 180.4 K + 8.9 K^2. End.

Remarks. (i) The first bracket does not depend on K; the linear term comes from the first-order tail below the threshold lambda = M = 4 C_D. (ii) K is universal but not explicit, because the Harnack constants of 4.2 and the constants A, M, k of Theorem 4.11 were not tracked, so K' is finite and not explicit. Numerically the true value of L sup|N_x - N_0| is at most 0.77 (check_roots.py, check_gen_roots.py). (iii) At half filling K = 1 and the formula gives K' <= 225; Corollary 3.2, by the C151 route, gives 70.2.

