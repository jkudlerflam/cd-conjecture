# Task B: a supersolution for the XXZ Bethe equations near the ferromagnet at general filling (Fermi edge present)

Goal: a rigorous proof of Theorem G below, or the sharpest rigorous partial result you can reach (for example the theorem on a range of fillings or under one explicitly stated, purely continuum hypothesis on the function r_Q, with everything else proved). Every step must be shown. Say plainly what is proved and what is not. Everything you need is defined here; assume nothing else. If you find an error in the reductions I state in Section 3, say so.

## 1. Definitions

eta in (0, pi/4], gamma = pi - eta, b = pi/gamma, q = 1 - 2 eta/pi.

    p(x)     = 2 atan( tan(eta/2) tanh x ),   p(+-inf) = +-eta,  strictly increasing
    theta(x) = -2 atan( tanh x / tan eta ),   odd, decreasing, theta(+-inf) = -+(pi - 2 eta)
    a1(x)    = p'/(2 pi) = sin(eta)/(pi (cosh 2x + cos eta)) > 0,   int a1 = eta/pi,   |a1'| <= 2 a1
    T(x)     = -theta'/(2 pi) = sin(2 eta)/(pi (cosh 2x - cos 2 eta)) = sin(2 eta)/(2 pi (sinh^2 x + sin^2 eta)) > 0,
               even, strictly decreasing in |x|, int T = q < 1, T(0) = cot(eta)/pi
    Theta(x) = int_0^x T = (1/pi) atan( tanh x / tan eta ),  Theta(inf) = q/2;   int_w^inf T = (1/pi)[atan(tan(eta) coth w) - eta].
    Fourier: T^(k) = sinh((pi - 2 eta) k/2)/sinh(pi k/2),  a1^(k) = sinh(eta k/2)/sinh(pi k/2).
    T = (1/(2 pi i)) [coth(x - i eta) - coth(x + i eta)], analytic in |Im x| < eta, Re T > 0 there.

Restricted continuum problem. For Q in (0, inf], (T_Q f)(x) = int_{-Q}^{Q} T(x - v) f(v) dv on [-Q, Q], and r_Q = (I - T_Q)^{-1} a1 (Neumann series, positive kernel, norm <= q). m(Q) = int_{-Q}^{Q} r_Q is continuous, strictly increasing from 0 to 1/2. For Q = inf: r_inf(v) = r_0(v) = b/(2 pi cosh(b v)).

Facts about r_Q that are proved and may be cited:
 (F1) r_Q > 0, even, nonincreasing on [0, Q]; r_Q <= r_0 pointwise; r_Q(x) is nondecreasing in Q for fixed x.
 (F2) 0 <= r_0(x) - r_Q(x) <= r_0(2Q - x) + r_0(2Q + x) on [-Q, Q] (reflection supersolution).
 (F3) The resolvent (I - T_Q)^{-1} = I + R_Q has a nonnegative kernel; (I - T_Q) 1 >= int_{|v|>Q} T(x - v) dv.
 (F4) The extension s(u) := a1(u) + int_{-Q}^{Q} T(u - v) r_Q(v) dv, u in R, equals r_Q on [-Q, Q], is positive, even, unimodal on R (convolution of even unimodal functions), analytic in |Im u| < eta with Re s > 0 there, and |s'| <= C s/eta on R for a universal C.
 Nothing else about r_Q is known rigorously (in particular no Wiener-Hopf edge asymptotics, no derivative bound on the unit scale, no convexity in Q). Heuristically r_Q vanishes like a square root at the edge on scales between eta and 1, i.e. r_Q(Q - u) ~ c sqrt(max(u, eta)).

Bethe system. Integers L >= 2, 1 <= n <= L/2, nu = n/L, I_j = j - (n+1)/2. F_j(x) = L p(x_j) - sum_{k=1}^n theta(x_j - x_k) - 2 pi I_j on R^n. It has a unique real zero x* (increasing, symmetric under x -> -x). Comparison principle (proved): F_j(y) >= F_j(x) for all j implies y_j >= x_j for all j. Variant with nodes at +infinity (proved): if some trial nodes are at +infinity and for the finite ones Ftilde_j := L p(y_j) - sum_{k finite} theta(y_j - y_k) - n_inf (pi - 2 eta) - 2 pi I_j >= 0, then x*_j <= y_j for all finite j.

Continuum reference at filling nu: Q = Q_nu with m(Q) = nu (Q = inf iff nu = 1/2). N_0(v) = int_{-Q}^{min(v,Q)} r_Q for v >= -Q, 0 for v < -Q (the cumulative function of r_Q 1_[-Q,Q], with N_0 = nu for v >= Q). The counting function Z_0(u) = p(u)/(2 pi) + (T*N_0)(u) - q nu/2 satisfies Z_0' = s (F4), Z_0 = N_0 - nu/2 on [-Q, Q], and is strictly increasing on R with range (-nu/2 - x0, nu/2 + x0), x0 = (eta/(2 pi))(1 - 2 nu). Put t(v) = L (Z_0(v) + nu/2), so t = L N_0 on [-Q, Q], t(-Q) = 0, t(Q) = n, t' = L s.

Sawtooth s#(u) = u - floor(u) - 1/2 (I write s# to distinguish it from the density s).

## 2. Target

Theorem G (wanted). There are universal constants K and eta_0 > 0 such that for all eta in (0, eta_0], all L >= 2 and 1 <= n <= L/2, the root counting function N_x(v) = (1/L) #{j : x*_j <= v} satisfies |N_x(v) - N_0(v)| <= K/L for every real v. (Any K is fine. By the reflection symmetry it suffices to prove N_x >= N_0 - K/L.)

At half filling this is proved with K = 1 (statement and proof sketch in Section 4). The task is nu < 1/2.

## 3. The reduction I have (please check it, then prove the remaining inequality or replace the route)

Trial configuration with shift kappa0 = K0 + 1/2, K0 a nonnegative integer: nodes y_j defined by t(y_j) = t_j := j + K0, j = 1..n, when t_j is in the range of t; otherwise y_j = +infinity. So nodes j <= n - K0 lie in [-Q, Q] (node j = n - K0 is at Q) and the K0 top nodes lie beyond Q or at infinity, at the points where t = n + i, i = 1..K0. Let N_y count the finite nodes. Exact residual identity for finite j:

    Ftilde_j/(2 pi L) = kappa0/L - [T*(N_0 - N_y)](y_j),

and L(N_0 - N_y)(v) equals:  0 for v < -Q;  kappa0 + s#(t(v) - t_j) - (a nonnegative integer deficit, present only where 0 <= t(v) < K0) for -Q <= v <= Q;  K0 - min(floor(t(v) - n), K0) for v > Q.
So, for a node y = y_j in [-Q, Q], Ftilde_j >= 0 follows from the inequality

  (**)   int_{-Q}^{Q} T(y - v) s#(t(v) - t(y)) dv
           <=  (2 eta/pi) kappa0  +  kappa0 int_{v < -Q} T(y - v) dv  +  int_{v > Q} T(y - v) [ 1/2 + min(floor(t(v) - n), K0) ] dv.

For the finite nodes beyond Q the inequality Ftilde_j >= 0 holds automatically (pairing v = y +- w: on the right the integrand of L T*(N_0 - N_y) minus kappa0 T is at most -T(w)/2, on the left at most +T(w)/2).
If (**) holds at all nodes in [-Q, Q], then x*_j <= y_j for finite j, hence N_x >= N_y >= N_0 - (kappa0 + 1/2)/L, which is Theorem G with K = K0 + 1.

Rewriting of the left side of (**). For y in [-Q, Q] let d = Q - y, d' = Q + y, and for m >= 0 let w_+(m), w_-(m) be the matched distances, N_0(y + w_+) - N_0(y) = m = N_0(y) - N_0(y - w_-), capped at d and d' respectively once the edge is reached. With Dt_y(m) = Theta(w_+(m)) - Theta(w_-(m)) one has
    LHS of (**) = int_0^inf s#(L m) dDt_y(m)      (Riemann-Stieltjes),
so |LHS| <= (1/2) TV(Dt_y), and on any stretch where Dt_y is smooth the contribution is a trapezoid-rule error (sum over integer points of L m minus integral). TV(Dt_y) does not depend on L.

Numerics (floating point, FFT solver for r_Q, eta from 0.14 down to 0.0014, L up to 1024, eight fillings): with X := [LHS of (**) - (1/2) int_{|v|>Q} T(y - v) dv]/eta evaluated at every node in [-Q, Q], max X = 0.13, attained close to half filling in the sparse left tail; near the Fermi edge X is negative and large (-0.2 to -10); in the dense interior at distance d of order one from the edge X is about -(1/(2 eta)) int_{|v|>Q} T(y - v) dv plus an exponentially small term. So (**) appears to hold already with kappa0 = 1/2 (K0 = 0), where it needs X <= 1/pi = 0.318, and with a wide margin for larger kappa0. Since any K is acceptable, kappa0 may be taken as large as convenient, and constants may be lossy everywhere.

Where the difficulty is. (i) Interior: bounding the trapezoid error of Dt_y by C eta needs some regularity of r_Q on scales between the lattice spacing and 1 (at half filling an exact structure gives unimodality of D_y, see Section 4). A crude first-order bound produces log(1/eta) and is not acceptable. (ii) Fermi edge, when the edge is neither dense nor sparse: the density s drops across Q on the scale eta (outside, s(Q + u) is about r_Q(Q) eta/(pi u) for eta < u < 1), so the lattice is strongly asymmetric there and the estimate must be sign-aware; the available gain on the right is (1/2 + number of trial nodes between Q and v) T, and on the left kappa0 T.

## 4. The half-filled case (proved; you may use it and imitate it)

nu = 1/2, Q = inf, N(v) = (1/pi) atan(e^{b v}), alpha = pi N(y), phi = pi m, a = 2/b in [3/2, 2). Then coth(b w_+/2) = z + d0 and coth(b w_-/2) = z - d0 with z = sin(2 alpha) cot(phi) and d0 = cos(2 alpha), so -D_y(m) = int_{z - d0}^{z + d0} h(x) dx with h(x) = a T(a u) sinh^2 u = (a sin(2 eta)/(2 pi)) sinh^2 u/(sinh^2(a u) + sin^2 eta) for x = coth u > 1 and h = 0 for x <= 1. h is unimodal (because sinh^2(a u) is convex as a function of sinh^2 u), so D_y is unimodal, TV(D_y) = 2 sup|D_y| <= 2 |d0| sup h, and sinh^2(a u) + sin^2 eta >= 2 a sinh^2 u gives sup h <= sin(2 eta)/(4 pi). Hence |int s#(L m) dD_y| <= eta/pi for every real y and L, the trial with kappa0 = 1/2 (top node at infinity) is a supersolution, and |N_x - N| <= 1/L.

## 5. What to return

Either a full proof of Theorem G for nu < 1/2, or a theorem of the form "Theorem G holds provided (H)", where (H) is an explicit statement about the continuum function r_Q alone (no L, no Bethe roots, no trial lattice), as weak and as checkable as you can make it, together with a complete proof of the implication and an honest account of why (H) is not proved. If you can prove (H) from (F1)-(F4), do so. If you can prove Theorem G on a sub-range (for instance lambda_Q := eta L r_Q(Q) >= Lambda_0, or n <= c L, or L - 2n <= C), state and prove that. You are free to choose a different trial configuration or a different reference density (the residual identity holds for any reference), to use larger kappa0, non-constant shifts, nodes at infinity, or any other device. Identify exactly which inequality fails if a route fails.
