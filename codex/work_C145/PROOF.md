# C145 proof notes: the ferromagnetic strip, the shell, and the continuum bounds

2026-11-07. Written by a Claude Opus agent, not by Codex. J=1 in all analytic sections. Tags: [Proved] is a mathematical deduction, [Proved, cited] uses an earlier run (named, with its referee status), [Numerical] is floating point without enclosures, [Heuristic] has an unproved step, [Assumed] is a hypothesis.

Notation is that of C136 Section 1, C138 Section 1 and C142 Part II. eta = pi - gamma in (0, pi/2], Delta = -cos eta, q = 1 - 2eta/pi, B = pi/(2eta), a = pi sin eta <= pi min(1, eta), b = pi/gamma in [1, 2), r_0(x) = 1/(2 gamma cosh(bx)), T >= 0 even and decreasing on x > 0 with integral q, r = r_Q = (I - T_Q)^(-1) a1, Z = Z_Q = (I - T_Q)^(-1) 1, eps_Q = (I - T_Q)^(-1)(e - h(Q)), v_Q = eps_Q'/r_Q, Q = Q_L the counting boundary, sigma = sigma_L, R_L = mu - 1_[-Q,Q] sigma dx, and

    lambda(x) = eta L r_Q(x),        the local crossover parameter.

The cells X_j are the x-images of [I_j/L - 1/(2L), I_j/L + 1/(2L)] under z_L. They tile [-Q, Q] and carry sigma-mass 1/L each (C142 Lemma II.1).

Summary.

- [Proved] Part A. Route (c) of the spec, carried out. In C138 Theorem B with the C142 log bound inserted, the core term is already bounded uniformly in eta; the logarithm enters only through the first-order charge of the shell S = {M_0 <= lambda < 4C_D(eta)}, which holds of order log(1/eta)/eta roots, so the shell cannot be charged first order and cannot be counted away. Theorem A.2 replaces (D), (D1) by a single shell hypothesis (G) (comparability and an integrated derivative bound on S only), and proves the remaining core piece log-free.
- [Proved] Part B. (W2) holds uniformly in eta and Q with kappa(d) = 1 - 4 exp(-2d) at distance d from the Fermi point (reflection supersolution). (W1) and (W3) are reduced to two one-line continuum statements, monotonicity of v_Q and an edge lower bound on Z_Q(Q), with proved implications and explicit constants. These two are not proved.
- [Proved] Part C. Referee report on C138 Theorem B. Correct as stated; two remarks.
- [Assumed] Part D. Why the shell is not closed here, and the exact inequality left, with its numerical margin.
- Part E. Assembly.

# Part A. What Theorem B needs from the shell

## A.1 Where the logarithm sits

[Proved, bookkeeping on C138 (4.2)] With (D), (D1) available only in the C142 form (constants C_D(eta), C_D1(eta), C142 Theorem II.6, valid at every x in [-Q,Q]) one must take the core threshold M = M_1(eta) := 4C_D(eta), since comparability |sigma - r| <= r/4 then holds exactly on K_1 = {lambda >= M_1}. The five terms of C138 (4.2) then behave as follows.

1. Core: (8a/(9M_1))[C_D(eta) + C_D1(eta)(1 + 1/(kappa b))/eta] <= (2pi/9)[1 + (C_D1(eta)/C_D(eta))(1 + 1/kappa)] using a/eta <= pi. The ratio C_D1(eta)/C_D(eta) is bounded on (0, pi/2): both grow like (const) log(1/eta) as eta -> 0 (C142 Theorem II.6), and at the free end both vanish linearly in cos eta. So the core term is bounded uniformly. Section A.3 makes this explicit.
2. Var(v_Q) term C_W eta/6 and 3. straddling cells 16a/3: uniform.
4. Mismatch C_e eta^2 (B-1)^2/2 <= pi^2 C_e/8: uniform.
5. Tail: 4aM_1(1 + 2/(kappa b))/eta, of order 16pi C_D(eta), that is (64) log(1/eta).

So the logarithm comes only from the first-order charge of the cells in the tail T_1 = {lambda < M_1}, and within it only from the shell S = {M_0 <= lambda < M_1}, because the part {lambda < M_0} costs 4aM_0(...)/eta with M_0 fixed.

[Proved for the continuum density; the root count differs from it by the cumulative error] The shell cannot be counted away. At half filling, int_y^inf r_0 = atan(e^(-by))/pi (C142 Lemma I.1) and r_0(y) = e^(-by)/(gamma(1 + e^(-2by))), so int_y^inf r_0 = (gamma/pi) r_0(y)(1 + O(r_0(y)^2)). The continuum mass of S is therefore 2 gamma (M_1 - M_0)/(pi eta L) to leading order, that is about 2(M_1 - M_0)/eta roots, which is of order log(1/eta)/eta. Each is charged a first-order cell error of size up to 2a/L in C138's tail bound. Route (c) as phrased in the spec ("the shell contains few roots") therefore fails; Theorem A cannot help either, since it is an O(eta) statement about the whole sum. What route (c) does deliver is that only second-order control on S is needed, not (D) there, and not (D) anywhere with a log-free constant.

## A.2 The conditional theorem with a shell hypothesis

[Assumed] (G) Shell hypothesis, with constants M_0 > 0 and C_GS. On S = {x in [-Q,Q] : M_0 <= lambda(x) < M_1(eta)},

    sigma_L >= (3/4) r_Q     and     eta Var_S(r_Q/sigma_L) <= C_GS,

uniformly in eta in (0, pi/2), L, n (S is empty when lambda < M_0 everywhere or M_1 <= M_0).

(G) is weaker than C138's (D), (D1) restricted to S: it asks for comparability and for the total variation of g = r/sigma over S, not for a pointwise derivative bound, and nothing at all on {lambda >= M_1}.

[Proved, conditional on (W1), (W3), (G); (W2) is proved in Part B] Theorem A.2. For every -1 < Delta < 0, L >= 2, 0 <= n <= L/2,

    |F_L(n) - L e_bulk(Delta, n/L)| <= K_G/L,
    K_G = (pi/4)(C_G1 + C_GS) + pi C_W/12 + 4 pi M_0 A_tail + 16pi/3 + pi^2 C_e/8,

with A_tail = min_d [d + 2/(1 - 4e^(-2d))] = 3.9887 (at d = 1.575) and C_G1 = sup_eta (8/(9C_D(eta)))[eta C_D(eta) + C_D1(eta) A_core], A_core = min_d [d + 1/(1 - 4e^(-2d))] = 2.7177 (at d = 1.352). C_G1 <= 13.19 analytically, 10.11 numerically (Section A.3).

Proof. Use the exact identity C136 (4.3) and C138 Lemma 1 (re-derived in Part C) with core cells those whose x-image lies in K_0 = {lambda >= M_0} = K_1 union S. K_0 and K_1 are superlevel sets of the even function r, nonincreasing on [0,Q] (C136 Section 2), so they are centred intervals [-x_0, x_0] containing [-x_1, x_1].

Comparability on K_0. On K_1, C142 Theorem II.6 gives |sigma - r| <= C_D(eta)/(eta L) <= r/4. On S, (G) gives sigma >= 3r/4. So g = r/sigma <= 4/3 on K_0, and on K_1 also g >= 4/5.

Core. As in C138, Var_(K_0)(eps'/sigma) <= sup|v| Var_(K_0)(g) + sup_(K_0) g Var(v) <= 2a Var_(K_0)(g) + (4/3) C_W eta, with |v| <= 2a (C136 (5.2)). g is continuous, so Var_(K_0)(g) = Var_(K_1)(g) + Var_S(g). On K_1 the computation of C138 applies verbatim with C_D, C_D1 replaced by C_D(eta), C_D1(eta) and M = M_1, except that the split of int_(K_1) dx/r is made at distance d from the Fermi point and uses (W2) in the proved form r >= kappa(d) r_0 with kappa(d) = 1 - 4e^(-2d) (Theorem B.1): per side int dx/r <= (d + 1/(kappa(d) b))/rho_1 <= A_core/rho_1, rho_1 = M_1/(eta L). Hence

    eta Var_(K_1)(g) <= (16 eta/9)[(C_D(eta)/(eta L))(2/rho_1) + (C_D1(eta)/(eta^2 L)) 2A_core/rho_1]
                     = (8/(9C_D(eta)))[eta C_D(eta) + C_D1(eta) A_core] <= C_G1.

With (G), the core contributes at most (1/(8L))[2a(C_G1 + C_GS)/eta + (4/3) C_W eta] <= [(pi/4)(C_G1 + C_GS) + pi C_W/12]/L.

Straddling cells. At most one cell straddles each endpoint of K_0; on its K_0-part r <= (4/3) sigma, so the integral of |eps'| <= 2ar there is at most 8a/(3L). Two cells: 16a/(3L) <= (16pi/3)/L.

Tail {lambda < M_0}. As in C138 with rho_0 = M_0/(eta L) and the split at distance d: per side, the part within d of the Fermi point has r < rho_0 and length at most d, and the rest has r_0 <= r/kappa(d) < rho_0/kappa(d), so its r-integral is at most int_(r_0 < rho_0/kappa) r_0 <= 2rho_0/(kappa(d) b). With |eps'| <= 2ar and two sides, the tail is at most 4a rho_0 (d + 2/(kappa(d) b)) <= 4 pi M_0 A_tail/L.

Mismatch. C138 verbatim with (W3): at most pi^2 C_e/(8L).

Add the four pieces. At half filling (Q = infinity) the same holds with r = r_0, kappa = 1, no mismatch, and Var(v) = 2ab <= 8 eta. End.

[Numerical, constants_C145.json] With M_0 = 2, C_W = 8, C_e = 13/10 (C138's sampled bounds) and C_GS = 1: K_G = 129.4 and C = 5K_G + pi^2/4 = 649.6. With C_GS = 5: K_G = 132.6, C = 665.3. The constant is dominated by the tail term 4pi M_0 A_tail = 100.2.

## A.3 The proved core piece is log-free

[Proved] eta Var_(K_1)(r/sigma_L) <= (8/(9C_D(eta)))[eta C_D(eta) + C_D1(eta) A_core] = (8/9)[eta + (C_D1/C_D) A_core], with C_D(eta) = 2(1+q) eta calR(0) and C_D1(eta) = eta^2[4M_1^T + 4(2+q)T(0) calR(0)] from C142 Theorem II.6. The ratio is

    C_D1/C_D = 2 eta M_1^T/((1+q) calR(0)) + 2(2+q) eta T(0)/(1+q).

The second term is at most 4/pi, since (2+q)/(1+q) <= 2 and eta T(0) = eta cot(eta)/pi <= 1/pi. For the first, calR = T + T*calR >= T pointwise, so calR(0) >= T(0) = cot(eta)/pi, and with M_1^T = (2cos eta/pi)(0.32476/sin^2 eta + 0.25/sin eta) (C142 Lemma II.5),

    2 eta M_1^T/((1+q) calR(0)) <= 4 eta sin(eta)(0.32476/sin^2 eta + 0.25/sin eta) = 4(0.32476 eta/sin eta + 0.25 eta) <= 4 * 0.9026 = 3.61,

using eta/sin eta <= pi/2 and eta <= pi/2. So C_D1/C_D <= 3.61 + 4/pi = 4.88 on (0, pi/2), and

    C_G1 <= (8/9)[pi/2 + 4.88 * 2.7177] = 13.19   (analytic).

[Numerical, constants_C145.json] The exact expression on 400 log-spaced points of [1e-7, pi/2 - 1e-3] (nine points where the quadrature of calR(0) returned NaN were dropped) is 2.40 at eta = 1e-7 and rises to 10.11 at the free end. Theorem A.2 uses the analytic 13.19 in the stated constant and 10.11 only in the numerical evaluation.

The point is that the variation of g over the deep core is controlled with a log-free constant even though the pointwise (D) there is proved only with the factor log(1/eta). The log in the pointwise bound is compensated by the threshold M_1 = 4C_D(eta) in the denominator.

# Part B. The continuum bounds (W1)-(W3) of C138 Section 4.2

## B.1 (W2), proved

[Proved] Theorem B.1 (reflection bound). For every eta in (0, pi/2], Q in (0, infinity) and x in [-Q, Q],

    0 <= r_0(x) - r_Q(x) <= r_0(2Q - x) + r_0(2Q + x),

and consequently, with d = Q - |x|,

    r_Q(x) >= (1 - 4 e^(-2bd)) r_0(x) >= (1 - 4 e^(-2d)) r_0(x).

In particular (W2) holds with kappa = 1 - 4e^(-2) = 0.4587 at d >= 1, and with kappa(d) = 1 - 4e^(-2d) at every d > (log 4)/2 = 0.693, uniformly in eta and Q.

Proof. Let h = r_0 - r_Q on [-Q, Q]. Restricting the full-line equation r_0 = a1 + T*r_0 to [-Q, Q] gives r_0 = a1 + T_Q r_0 + g with g(x) = int_(|u|>Q) T(x - u) r_0(u) du, and r_Q = a1 + T_Q r_Q. Subtracting, (I - T_Q) h = g, so h = (I - T_Q)^(-1) g >= 0 (positive Neumann series, ||T_Q|| <= q < 1). Split g = g_+ + g_-, the parts with u > Q and u < -Q.

Let psi_+(x) = r_0(2Q - x), a reflection of a translate of r_0, so psi_+ solves the full-line equation psi_+ = a1(2Q - .) + T*psi_+ (T is even). On [-Q, Q],

    (I - T_Q) psi_+ = a1(2Q - x) + int_(|u|>Q) T(x - u) psi_+(u) du >= int_(u>Q) T(x - u) r_0(2Q - u) du.

For u > Q >= 0 one has |2Q - u| <= |u| (equivalent to 4Q^2 - 4Qu <= 0), and r_0 is even and decreasing in |x|, so r_0(2Q - u) >= r_0(u). Since T >= 0, (I - T_Q) psi_+ >= g_+ on [-Q, Q]. In the same way psi_-(x) = r_0(2Q + x) gives (I - T_Q) psi_- >= g_-. Hence (I - T_Q)(psi_+ + psi_- - h) >= 0, and positivity of (I - T_Q)^(-1) gives h <= psi_+ + psi_-.

For the ratio, cosh(bx) <= e^(b|x|) and cosh(bs) >= e^(bs)/2 for s >= 0. With |x| <= Q, both 2Q - x and 2Q + x are at least 2Q - |x| >= |x|, so r_0(2Q -+ x)/r_0(x) <= 2 e^(b|x| - b(2Q - |x|)) = 2 e^(-2b(Q - |x|)). Adding the two terms and using b >= 1 gives the second display. End.

Consistency checks. [Proved] At Q -> infinity the right side vanishes and r_Q -> r_0, as it must. At the free point T = 0 one has r_Q = r_0 = a1 on [-Q,Q], so h = 0, consistent. The bound is vacuous within distance 0.69 of the Fermi point, where C138 found r_Q(Q)/r_0(Q) of order eta^(1/2) numerically, so no bound of this form can hold up to the edge uniformly.

[Numerical, check_W2.json, check_W2_refine.json] 30 cases, eps = 1 + Delta in {0.5, 0.1, 1e-2, 1e-4, 1e-6}, Q in {0.3, 1, 1.5, 2, 4, 6}, Nystrom mesh N = 900. Both inequalities hold with margin, except at Q >= 4 for eps <= 1e-4, where the computed r_Q/r_0 also exceeds 1 (contradicting the proved r_Q <= r_0). Under refinement N = 900, 1800, 2700 these apparent violations shrink as 0.0105, 0.0025, 0.0007 (eps = 1e-6, Q = 4) and 0.0020, 0.00024, 0.00005 (eps = 1e-4, Q = 6), so they are mesh artifacts. min r_Q/r_0 at distance >= 1 is 0.762 (Q = 1, eps = 1e-6) against the proved 0.459, in agreement with C138's 0.762.

This also answers C142's remark: the positivity argument through the full-line kernel calR loses a factor 1/eta because calR has mass of order 1/eta, while the reflected source psi_+ has the right size by construction. No Wiener-Hopf factorization is needed for (W2).

## B.2 Reduction of (W1) and (W3)

[Proved] Lemma B.2. For x in [-Q, Q], |eps_Q'(x)| <= 2a tanh(2Q) r_Q(x), so |v_Q| <= 2a tanh(2Q) <= 2a.

Proof. a1(x) = sin(eta)/(pi(cosh 2x + cos eta)), so |a1'(x)| = 2 sinh(2|x|) a1(x)/(cosh 2x + cos eta) <= 2 tanh(2|x|) a1(x), since cos eta >= 0. By C136 (5.1), (I - T_Q) eps_Q' = e' = -a a1', so |eps_Q'| <= (I - T_Q)^(-1) |a a1'| <= 2a (I - T_Q)^(-1)[tanh(2|.|) a1] <= 2a tanh(2Q) (I - T_Q)^(-1) a1 = 2a tanh(2Q) r_Q, by positivity of the inverse and tanh(2|y|) <= tanh(2Q) on [-Q, Q]. End.

This sharpens C136 (5.2) by the factor tanh(2Q), which matters only in the dilute regime.

[Proved] Proposition B.3 (two sufficient conditions).

(i) If v_Q is nondecreasing on [0, Q], then Var_[-Q,Q](v_Q) = 2v_Q(Q) <= 4a tanh(2Q) <= 4pi eta, that is (W1) with C_W = 4pi.

(ii) e_bulk''(nu) <= a tanh(2Q)/Z_Q(Q)^2 at nu = m(Q). Hence, if eta Z_Q(Q)^2 >= c_Z tanh(2Q) for all Q, then (W3) holds with C_e = pi/c_Z.

Proof. (i) v_Q is odd (eps_Q even, r_Q even), so v_Q(0) = 0, and v_Q >= 0 on (0, Q] (C136 Section 5, eps_Q' > 0 there). If it is nondecreasing there, Var over [-Q,Q] is 2(v_Q(Q) - v_Q(0)) = 2v_Q(Q), and Lemma B.2 bounds it, with a <= pi eta. (ii) C136 (5.3) gives e_bulk'' = eps_Q'(Q)/(2 r_Q(Q) Z_Q(Q)^2) = v_Q(Q)/(2Z_Q(Q)^2), and Lemma B.2 gives v_Q(Q) <= 2a tanh(2Q). With a <= pi eta the last claim follows. End.

[Proved, structural remark] Both conditions are statements about the edge. The derivative identities (I - T_Q) r' = a1' - r(Q)J_0 and (I - T_Q) Z' = -Z(Q) J_0, J_0(x) = T(x - Q) - T(x + Q), obtained by differentiating the restricted equations and integrating by parts, show that r_Q'/r_Q and Z_Q'/Z_Q carry the same edge-singular part r(Q)J/r and Z(Q)J/Z, J = (I - T_Q)^(-1) J_0, which cancels exactly at x = Q in v_Q(Q) = a(Z'/Z - r'/r)(Q) (this is C136 (5.1) seen from the other side). The remaining edge statement in (ii) is the size of Z_Q(Q): it is of order 1 for Q << eta and of order (min(Q,1)/eta)^(1/2) otherwise [Heuristic, Wiener-Hopf square-root edge], which is what c_Z > 0 expresses. Known lower bounds are much weaker: Z_Q(Q) >= 1/(1 - q/2) from Z(Q) = 1 + int T(Q - y) Z(y) dy >= 1 + Z(Q) int_0^(2Q) T, and Z_Q >= r_Q/a1(0) from (I - T_Q)(r_Q/a1(0)) = a1/a1(0) <= 1, which loses cosh(bQ).

[Assumed] (W1') v_Q is nondecreasing on [0, Q], for every eta and Q. [Assumed] (W3') eta Z_Q(Q)^2 >= c_Z tanh(2Q) for every eta and Q, with c_Z > 0 independent of both.

[Numerical, check_W13.json] 56 cases, 1 + Delta in {0.5, 0.1, 1e-2, ..., 1e-6}, Q in {0.003, 0.01, 0.03, 0.1, 0.3, 1, 2, 3}, mesh N = 900. (W1'): v_Q is nondecreasing on the mesh in 55 cases; the one exception (eps = 1e-6, Q = 3) has a relative dip of 3.4e-6, at the mesh-noise level of that underresolved case. Var(v_Q)/eta <= 7.81 (proved under (W1'): 4pi = 12.57). v_Q(Q)/(2a tanh 2Q) <= 0.750 (proved <= 1). (W3'): the minimum of eta Z_Q(Q)^2/tanh(2Q) over the Q grid is 1.570, 1.459, 1.149, 0.965, 0.895, 0.833, 0.808 at eps = 0.5, 0.1, 1e-2, ..., 1e-6, attained at Q/eta between 20 and 70. In the dilute window the values depend on Q/eta only (for example 1.080, 1.067, 1.080 at Q/eta = 2.1 to 2.2 for eps = 1e-6, 1e-5, 1e-4), consistent with a limit function of Q/eta as eta -> 0 with infimum near 0.8. Below we use c_Z = 3/4, so C_e = 4pi/3 = 4.19, against the sampled e_bulk''/eta^2 <= 1.24.


# Part C. Referee report on C138 Theorem B (C138 PROOF Sections 4.1 to 4.3)

I re-derived every line as I used it in Part A. Verdict: correct as stated, conditional on (D), (D1), (W1)-(W3) as formulated in C138 Section 4.2. No error, no hidden dependence on L, Q or eta. Two remarks and one improvement follow.

C.1 Lemma 1 (three-region quadrature). The midpoint Peano kernel for f with f' of bounded variation on a cell of length h = 1/L is k(s) = (h/2 - |s - t_j|)^2/2 in [0, h^2/8], and int_cell f - h f(t_j) = -int k df'. So |f(t_j) - L int_cell f| <= Var_cell(f')/(8L). Variation is invariant under the monotone change of variable t = z_L(x), and df/dt = eps_Q'/sigma_L. For the other cells, |f(t_j) - L int_cell f| <= osc_cell f <= int_cell |f'| dt = int_(X_j) |eps_Q'| dx. Correct. The cells tile [-Q, Q] because z_L(+-Q_L) = +-nu/2 (C136 Section 4), and at Q = infinity they tile R.

C.2 Core. (D) and r >= rho = 4C_D/(eta L) give |sigma - r| <= r/4, hence g = r/sigma in [4/5, 4/3]. Product rule Var(vg) <= sup|v| Var(g) + sup|g| Var(v), with sup|v| <= 2a from C136 (5.2) (refereed, C136 REFEREE Section 5). g' = [r'(sigma - r) - r(sigma' - r')]/sigma^2 is an exact identity, and sigma^2 >= 9r^2/16. int_(K_c) |r'|/r^2 = 2[1/r(x_c) - 1/r(0)] <= 2/rho uses that r is even and nonincreasing on [0, Q] (C136 Section 2, refereed). The split of int_(K_c) dx/r and the bound int_0^y dx/r_0 = (2 gamma/b) sinh(by) <= 1/(b r_0(y)) are correct. Arithmetic: Var(g) <= (32/(9M))[C_D + C_D1(1 + 1/(kappa b))/eta] and the core term (8a/(9M))[...] + C_W eta/6 follow.

C.3 Tail. Straddling cells: r <= (4/3) sigma on the K_c-part and cell mass 1/L give 8a/(3L) each. The tail integral uses r < rho on T_c, r >= kappa r_0 on {Q - |x| >= 1}, r <= r_0 (C136 (2.2)), and int_y^inf r_0 <= e^(-by)/(gamma b) <= 2r_0(y)/b (since e^(-by) <= 2 gamma r_0(y)). If rho/kappa exceeds r_0(0) the whole half line is used, and int_0^inf r_0 = 1/4 <= 1/pi = 2r_0(0)/b, so the bound still holds. Result 4a rho(1 + 2/(kappa b)) per both sides, times L. Correct.

C.4 Mismatch. F_L - L e_bulk(nu) = L int eps_Q dR_L - L[e_bulk(nu) - e_bulk(m) - e_bulk'(m)(nu - m)] uses E(Q) = e_bulk(m(Q)) and h(Q) = e_bulk'(m) (C136 Section 5, refereed). (W3) and |nu - m| <= (B - 1)/L (C136 (4.2), refereed) give C_e eta^2 (B-1)^2/(2L). Correct.

C.5 (4.3) and the numerical value. a <= pi min(1, eta), b >= 1, (B - 1) eta <= pi/2, eta <= pi/2 give (4.3) term by term. With C_D = 1/2, C_D1 = 2, kappa = 7/10, C_W = 8, C_e = 13/10 the five terms are 7.48, 2.09, 96.9, 16.76, 1.60, sum 124.87 <= 125, and 5 * 125 + pi^2/4 = 627.5 <= 628 (constants_C145.json, KB_C138). Correct.

C.6 Half filling. At Q = infinity, v = ab tanh(bx) and Var(v) = 2ab = 2 pi^2 sin(eta)/(pi - eta). Remark: this equals C_W eta with C_W = 8 exactly at the free point eta = pi/2 and is below it for eta < pi/2 (sin(eta)/(eta(pi - eta)) increases from 1/pi to 4/pi^2), so C_W = 8 is sharp there and not merely a sampled bound. (W2) is not needed (r = r_0).

C.7 Corollary (4.4). Matched doubling: F_(2L)(2n) - 2F_L(n) <= K/(2L) + 2K/L at the common density nu, so L[...] <= 5K/2 in units of J, that is 5K M_0 with M_0 = J/2. The field, DM angle and twist reductions (pi^2 M/4) are those of C129 and C136 Section 9, refereed in C136 REFEREE (Assembly). Correct.

Remark 1. (D1) is used at every point of K_c, and K_c reaches the Fermi point whenever r_Q(Q) >= rho. C138's numerics for (D1) stopped a quarter kernel width short of Q. In Part A this does not matter: on K_1 the C142 Theorem II.6 derivative bound is proved at every x in [-Q, Q], and on the shell only the integrated (G) is used.

Remark 2. The dominant term 16pi C_D(1 + 2/kappa) is the first-order tail charge. Any improvement of Theorem B should go after the threshold M, not after the core. This is the observation behind Part A.

Improvement. With (W2) proved in the form kappa(d) = 1 - 4e^(-2d) (Theorem B.1), the factors (1 + 1/(kappa b)) and (1 + 2/(kappa b)) of C138 become A_core = 2.7177 and A_tail = 3.9887 after optimizing the split distance d. With C138's other sampled constants this gives K_* = 128.99 and C = 647.4 (constants_C145.json, KB_with_proved_W2), slightly above C138's 125 and 628 because the proved kappa is smaller than the sampled 0.7.

# Part D. The shell: what was tried and what is left

## D.1 Why the three routes of the spec do not close it here

Route (c) is Part A: it isolates (G) and shows that nothing cheaper than second-order control on S will do (the shell holds of order log(1/eta)/eta roots).

Route (a), multiscale bootstrap. [Proved, the computation behind the obstruction] Write sigma - r = int calR_Q(x, y) dR_L(y) (C142 Lemma II.2) and split the y-cells by their distance u = |y - x|. Charging a cell first order costs (1/L) osc_cell calR_Q(x, .), and for the full-line kernel calR'(u) behaves like -1/(pi eta u) for eta < |u| < 1 (C142 II.4), so the first-order charge of all cells in a set U, times eta L, is (1/pi) int_U du/|u| up to O(1). A second-order (trapezoid) charge of a cell is smaller by the factor (1/8)[1/(L|u|sigma) + |sigma'|/(L sigma^2)], so it wins wherever the cell is narrow compared with u and sigma varies little across the cell. With the a priori input of C142 Theorem II.6 the second factor is |sigma'|/(L sigma^2) <= b/(Lr) + C_D1(eta)/(eta^2 L^2 r^2) up to constants, and its du/|u| integral over eta < |u| < 1 gives C_D1(eta) log(1/eta)/lambda^2. This is bounded only for lambda >= C log(1/eta), which is the C142 threshold. Feeding back a log-free C_D1 does not help: the term is then C_D1 log(1/eta)/lambda^2, still unbounded on the shell. So the bootstrap needs a derivative bound that is better than (D1) away from the edge, of the form |sigma' - r'| <= C/(eta L) times (1 + 1/dist(y, sparse set)), which is a different statement and was not proved. [Heuristic] Physically, the sparse spikes act as sources of size 1/(eta L) localized on scale eta, and calR spreads each into (1/(pi eta)) log(1/|u|), whose derivative varies on the scale u, not eta. That is what the true sigma - r does numerically.

A continuity argument in eta (along the analytic branch of C142 Theorem I.4 and C136 Section 1) with the hypothesis "sigma >= r/2 on {lambda >= M}" also fails, for a geometric reason. [Proved, computation] At a point x with lambda(x) slightly above M, the cells just outside {lambda >= M} are not covered by the hypothesis and must be charged first order; their charge is (1/pi) log(l/max(D, eta)) with D the distance from x to the complement and l the length of the adjacent uncovered stretch within distance one. As lambda(x) decreases to M, D tends to 0 and the charge tends to (1/pi) log(1/eta), so the improved bound does not return the hypothesis at its boundary.

Route (b), direct second order with the exact sawtooth. [Proved] The identity int phi dR_L = -int psi(t) S(t) dt = int psi'(t) P(t) dt, with t = z_L(x), S the exact zero-mean sawtooth of period 1/L, P its antiderivative (0 <= -P <= 1/(8L^2)) and psi = d phi(y(t))/dt = phi'(y)/sigma(y), is exact and uses only monotonicity of z_L and the exact cell masses. It gives |int phi dR_L| <= Var_y(phi'/sigma)/(8L^2), which is C138 Lemma 1 in another form. The variation of phi'/sigma still needs a lower bound on sigma and control of sigma' on the unit neighbourhood of x, so route (b) meets the same obstruction as (a).

## D.2 The inequality left

[Assumed] (G) on S = {M_0 <= eta L r_Q < 4C_D(eta)}, with M_0 = 2:

    sigma_L >= (3/4) r_Q     and     eta Var_S(r_Q/sigma_L) <= C_GS.

Together with (W1') and (W3') of Part B, this is everything that remains for the strip; every other ingredient of Theorem A.2 is proved. (G) is implied by C138's (D), (D1) restricted to S, and it implies nothing on the deep core (where it is proved log-free in A.3) or on the sparse tail (which is charged first order at cost 4pi M_0 A_tail). The minimal energy-level statement behind it is |sum over cells in S of [eps_Q(x_j) - L int_(X_j) eps_Q sigma_L]| <= C/L.

[Numerical, check_shell_*.json, check_shell_big_*.json] Margins, eps = 1 + Delta in {1e-2, 1e-4, 1e-5, 1e-6}, L from 64 to 16384, n in {L/8, L/4, 3L/8, L/2-2, L/2}, with sigma_L sampled on a grid of spacing min(eta, 1/(L sigma_max))/25 on K_0 = {lambda >= 2}:

| 1+Delta | cases with K_0 nonempty | L max | min sigma/r on K_0 (needed >= 3/4) | max eta Var_(K_0)(r/sigma) | max eta Var_S(r/sigma) |
| --- | --- | --- | --- | --- | --- |
| 1e-2 | 15 | 2048 | 1.0000 | 0.0018 | 0.0015 |
| 1e-4 | 9 | 2048 | 1.0000 | 0.00018 | 0.00018 |
| 1e-5 | 8 | 8192 | 0.9950 | 0.00024 | 0.00024 |
| 1e-6 | 4 | 16384 | 1.0000 | 0.00011 | 0.00011 |

(The 0.995 is at n = L/2 - 2, L = 8192, where Q_L is large and the restricted continuum solve is mesh-limited, as in C138 Section 5; every other case has sigma >= r on K_0.) So (G) holds with C_GS of order 1e-3 and comparability far above 3/4, in every sampled case. K_0 is empty when eta L r(0) < 2, which at 1 + Delta = 1e-6 means L below about 8900; there the whole occupied interval is first-order tail, already covered by Theorem A.2. At 1 + Delta <= 1e-4 every sampled K_0 lies inside the shell (lambda(0) < M_1(eta) = 28 to 40), so these cases test (G) on the whole shell.

[Proved] Relation to C138's (D). If (D) holds with C_D <= 1/2 on {lambda >= 2}, then the comparability half of (G) holds. The numerical value of eta L sup|sigma - r| is at most 0.32 (C138, C142), so the comparability half of (G) has the same numerical margin as (D) at M_0 = 2, namely 0.32 against 0.5.

# Part E. Assembly

## E.1 The strip

[Proved, conditional on (G) with M_0 = 2, C_GS = 1, (W1') and (W3') with c_Z = 3/4] For every -1 < Delta <= 0 (so delta_0 = 1 is allowed: Theorem A.2 is uniform on the whole attractive side), every L >= 2 and 0 <= n <= L/2,

    |F_L(n) - L e_bulk(Delta, n/L)| <= 136.61/L,

with the five terms of Theorem A.2 equal to 11.15 (core, C_G1 = 13.19 analytic), 3.29 (Var v_Q, C_W = 4pi from (W1')), 100.25 (first-order tail {lambda < 2}), 16.76 (straddling cells), 5.17 (mismatch, C_e = 4pi/3 from (W3')). By matched doubling (C138 (4.4), C136 Section 9),

    L(E_(2L) - 2E_L) <= C_strip M,    C_strip = 5 * 136.61 + pi^2/4 = 685.5 <= 686,

for every bond of the U(1) class with -1 < J_z/J <= 0, every field and DM angle, every L >= 2 (constants_C145.json, C_strip_final). With the sampled values in place of the conservative ones (C_GS = 0.002, C_W = 8, c_Z = 0.808, C_G1 = 10.11) the same formula gives 661.6.

[Proved, unconditional] On the strip itself nothing new is unconditional beyond C142 Theorem II.6 and C138 Theorem A. The unconditional pieces of the strip theorem are the core (Section A.3), the tail and straddling cells, and (W2) (Theorem B.1).

## E.2 The whole class J_z/J in [-1, 1]

| region | constant C in L(E_2L - 2E_L) <= C M | source | status |
| --- | --- | --- | --- |
| 0 <= J_z/J <= 1 | 4078 | C142 Corollary I.10 | proved, not refereed |
| -1 + delta <= J_z/J <= 0 | 6 * 10^15 delta^(-25/3) | C136 | proved, refereed |
| J_z/J = -1 | 0 | C129 Section 8 | proved |
| -1 < J_z/J <= 0 | 686 | Theorem A.2 here (C138 Theorem B refereed in Part C) | conditional on (G), (W1'), (W3') |

Seams. At J_z/J = 0 the repulsive theorem (closed interval [0,1]) and either attractive statement (closed interval containing 0) use the same C129/C136 reductions (DM angle to flux, field cancellation in matched sectors, flux cost pi^2 M/4) and the same M >= M_0 = J/2 on |Delta| <= 1, so the constants combine by taking the maximum. At J_z/J = 1 the bound is reached by finite-matrix continuity at fixed L (C142). At J_z/J = -1 the defect is zero. So:

[Proved, conditional on (G), (W1'), (W3')] For the whole U(1) spin-1/2 class with -1 <= J_z/J <= 1, every field and DM angle and every L >= 2,

    L(E_(2L) - 2E_L) <= 4078 M,

the constant coming from the repulsive side (max(4078, 686) = 4078).

[Proved, unconditional] On [-1 + delta, 1] for every delta in (0, 1], the constant max(4078, 6 * 10^15 delta^(-25/3)) = 6 * 10^15 delta^(-25/3), and 0 at J_z/J = -1. The strip (-1, -1 + delta) is not covered unconditionally for any fixed delta.

## E.3 Consistency checks

* Free point Delta = 0 [Proved]: T = 0, C_D = 0, M_1 = 0, K_1 = [-Q, Q], g = 1, and Theorem A.2 reduces to the free midpoint estimate; (W2) is trivial (r_Q = r_0 = a1), v_Q is monotone, Var(v)/eta = 8 at half filling (C138's equality case).
* Half filling [Proved]: Q = infinity, r = r_0, kappa = 1, no mismatch; Theorem B.1 is vacuous and not needed.
* Dilute limit [Proved]: Lemma B.2 gives |v_Q| <= 4aQ, so Var(v_Q) and the bound of Proposition B.3(ii) vanish linearly in Q as Q -> 0 [Heuristic: for Q << eta, T_Q is negligible and e_bulk'' is about pi^2 nu] (sampled e_bulk''/eta^2 = 0.010 at eps = 0.5, Q = 0.003).
* Delta -> -1 [Numerical]: e_bulk''/eta^2 -> 1.00 at Q >= 1 (check_W13.json), matching C138; min r_Q/r_0 at distance one converges to 0.762.
* (D) and (G) [Numerical]: sigma_L >= r_Q on K_0 in all sampled cases, which is stronger than the needed 3/4. No proof of the sign is offered.
