# Task D: the sparse right Fermi edge (last missing case of a counting bound for XXZ Bethe roots near the ferromagnet)

Goal: a rigorous proof of Theorem G below in the one regime that is still open (Section 4), with every step shown; constants may be universal and non-explicit, and may be lossy. If you cannot close it, prove as much as you can (for example the many-node part of Section 5) and isolate the remaining inequality as sharply as possible, saying exactly what input is missing. Everything needed is stated here. Sections 2 and 3 are proved and may be used without proof. If you find an error in a formula of Sections 1 to 3, say so.

## 1. Definitions

eta in (0, 1/64], gamma = pi - eta, q = 1 - 2 eta/pi.

    p(x)     = 2 atan( tan(eta/2) tanh x ),   strictly increasing, p(+-inf) = +-eta
    theta(x) = -2 atan( tanh x / tan eta ),   odd, decreasing, theta(+-inf) = -+(pi - 2 eta)
    a1(x)    = p'/(2 pi) = sin(eta)/(pi (cosh 2x + cos eta)) > 0
    T(x)     = -theta'/(2 pi) = sin(2 eta)/(2 pi (sinh^2 x + sin^2 eta)) > 0, even, decreasing in |x|, int T = q, T(0) = cot(eta)/pi,
               T(w) <= eta/(pi w^2), |T'(w)| <= C T(w)/(w + eta), |T''(w)| <= C T(w)/(w + eta)^2 for 0 <= w <= 1
    Theta(x) = int_0^x T = (1/pi) atan(tanh x/tan eta), Theta(inf) = q/2,   B(d) := int_d^inf T = q/2 - Theta(d) (about eta/(pi d) for eta << d << 1).

Restricted continuum problem: for 0 < Q < inf, r_Q is the solution on [-Q, Q] of r_Q(x) = a1(x) + int_{-Q}^{Q} T(x - v) r_Q(v) dv. Its extension s(u) = a1(u) + int_{-Q}^{Q} T(u - v) r_Q(v) dv is defined on R and equals r_Q on [-Q, Q]. m(Q) = int_{-Q}^{Q} r_Q increases continuously from 0 to 1/2.

Bethe system: integers L >= 2, 1 <= n < L/2, I_j = j - (n+1)/2, F_j(x) = L p(x_j) - sum_{k=1}^n theta(x_j - x_k) - 2 pi I_j on R^n. It has a unique real zero x* (increasing, symmetric under x -> -x). The Jacobian of F is, at every configuration, a symmetric matrix with nonpositive off-diagonal entries dF_j/dx_k = -2 pi T(x_j - x_k) and row sums 2 pi L a1(x_j) > 0 (an M-matrix).

Comparison principle (proved): if F_j(y) >= F_j(x) for all j then y_j >= x_j for all j. Variant (proved): trial nodes may be placed at +infinity; if Ftilde_j := L p(y_j) - sum_{k finite} theta(y_j - y_k) - n_inf (pi - 2 eta) - 2 pi I_j >= 0 at every finite node, then x*_j <= y_j at every finite node. Consequence (proved): F is an M-function, so for any subset A of indices and any fixed values of the other coordinates, the subsystem {F_j = 0, j in A} has a unique solution, which is nondecreasing in the fixed coordinates and obeys the same comparison principle.

Continuum reference: nu = n/L, Q = Q_nu with m(Q) = nu, N_0(v) = int_{-Q}^{v} r_Q clipped to [0, nu] outside [-Q, Q], t(v) = L int_{-Q}^{v} s (so t = L N_0 on [-Q, Q], t(-Q) = 0, t(Q) = n, and t continues to increase beyond Q with t' = L s). Sawtooth s#(u) = u - floor(u) - 1/2, fractional part {u} = u - floor(u). N_x(v) = (1/L) #{j : x*_j <= v}.

Theorem G (wanted). There are universal K, eta_0 > 0 with |N_x(v) - N_0(v)| <= K/L for all real v, all eta <= eta_0, L, n. By reflection it is enough to prove N_x >= N_0 - K/L, i.e. to exhibit a supersolution y (in the sense of the variant above) with N_y >= N_0 - K/L.

## 2. Proved facts about r_Q and s

 (F1) r_Q > 0, even, nonincreasing on [0, Q]; s is even and unimodal on R; r_Q <= r_0 := 1/(2 gamma cosh(pi x/gamma)).
 (F2) Harmonic representation. With P_h(x) = sin(2h)/(pi (cosh 2x - cos 2h)) (so P_eta = T), f = r_Q 1_[-Q,Q] and g(x) = 1/(2 eta cosh(pi x/eta)), the function U(x, h) = (P_h * f)(x) + (P_{pi/2 - h} * g)(x) is positive and harmonic in the strip 0 < h < pi/2, with U(x, 0) = f(x), U(x, eta) = s(x). Hence U(x, 0) = U(x, eta) for |x| < Q, and V(x, h) = (1/eta) int_h^{h+eta} U(x, z) dz is positive harmonic with zero normal derivative on the segment {h = 0, |x| < Q}, so its even reflection in h is positive harmonic in a neighbourhood of that segment.
 (F3) Regularity on the distance-to-edge scale (from (F2), Harnack and derivative estimates on discs). With delta(y) = Q - |y| and rho(y) = min(1, eta + delta(y)):
        |r_Q'(y)| <= C r_Q(y)/rho(y),   |r_Q''(y)| <= C r_Q(y)/rho(y)^2   for |y| <= Q,   and |s'| <= C s/eta on R.
      In particular r_Q(v)/r_Q(y) is between (rho(v)/rho(y))^C and its inverse when v, y are on the same side, and s(Q + u) >= e^{-C} r_Q(Q) for 0 <= u <= eta.
 Not known: any explicit constant in (F3); the edge exponent (heuristically r_Q(Q - u) is like sqrt(max(u, eta)) for u < 1); convexity or concavity of r_Q; monotonicity of r_Q'/r_Q.

## 3. The trial, the supersolution inequality, and what is proved

Standard trial with shift kappa = K0 + 1/2 (K0 >= 0 integer, n > K0): nodes y_j with t(y_j) = j + K0 when this level is attained, otherwise y_j = +infinity. Node j = n - K0 is at Q; the K0 top nodes are beyond Q at the points where t = n + i, or at infinity.

For a node y = y_j in [-Q, Q] put delta = Q - y, delta' = Q + y, tau_+(w) = t(y + w) - t(y) for 0 <= w <= delta, tau_-(w) = t(y) - t(y - w) for 0 <= w <= delta' (both are integers at the caps, because t(y), t(Q) = n and t(-Q) = 0 are integers), and

    S(y) := int_{-Q}^{Q} T(y - v) s#(t(v) - t(y)) dv = int_0^delta T(w) {tau_+(w)} dw - int_0^{delta'} T(w) {tau_-(w)} dw + (Theta(delta') - Theta(delta))/2.

Exact criterion: Ftilde_j >= 0 at the node y_j in [-Q, Q] if

    (**)  S(y) <= (2 eta/pi) kappa + kappa B(delta') + int_{v > Q} T(y - v) [1/2 + min(floor(t(v) - n), K0)] dv,

and the last integral is at least B(delta)/2. Finite nodes beyond Q satisfy Ftilde_j >= 0 automatically. If (**) holds at all nodes in [-Q, Q] then N_x >= N_y >= N_0 - (K0 + 1)/L.

Equivalent form of (**) with only the guaranteed gain on the right:

    (*)   int_0^delta T(w) {tau_+(w)} dw  <=  int_0^{delta'} T(w) {tau_-(w)} dw + (2 eta/pi) kappa + (kappa + 1/2) B(delta').

Proved:
 (P1) Log-free quadrature bound: |S(y)| <= C_S eta/rho(y) for all y in [-Q, Q]. Proof idea: S(y) = int s#(L m) dD(m) with D(m) = Theta(w_+(m)) - Theta(w_-(m)), w_+- the matched distances (equal r_Q-mass m to the right and to the left of y, capped at the edges); in the rescaled mass u = m/r_Q(y), for u <= u_0 = c rho(y), w_+-(u) = u + A_+- u^2 + O(u^3/rho^2) with A_+- = -+ r_Q'(y)/(2 r_Q(y)), and T(w) w' = T(u) + A (u^2 T(u))' + O(u^2 T(u)/rho^2); the function u^2 T(u) = (sin(2 eta)/(2 pi)) u^2/(sinh^2 u + sin^2 eta) has total variation at most sin(2 eta)/pi; the rest of each branch has variation at most its remaining kernel mass C eta/u_0.
 (P2) Theorem G holds if lambda_Q := eta L r_Q(Q) >= Lambda (dense edge), for universal K0, Lambda: then t(Q + eta) - n >= K0, the bracket in (**) equals kappa on v >= Q + eta, and (P1) is absorbed by kappa [2 eta/pi + B(delta') + B(delta + eta)] since eta/rho <= C[eta + B(delta + eta)].
 (P3) For y <= 0 (left half) (**) holds for kappa large: there rho(y) is governed by the left edge and (P1) is absorbed by kappa [2 eta/pi + B(delta')].
 (P4) For y > 0 with delta >= 1, (**) holds for kappa large.

## 4. What is open

The case lambda_Q < Lambda (the Fermi edge is not dense on the kernel scale) and nodes 0 < y < Q with delta = Q - y < 1. There (P1) gives C_S eta/delta while the guaranteed gain is B(delta)/2, about eta/(2 pi delta), with a fixed coefficient that does not grow with kappa (in this regime the K0 trial nodes beyond Q are far from Q, because s is small there, and cannot be moved closer: the top node's own inequality forbids it). So the sign of S must be used. Numerically (floating point, eta from 0.14 to 0.0014, L up to 1024, eight fillings) S(y) - B(delta)/2 is negative at every node near the right edge, between -0.2 eta and -10 eta, and the criterion holds already with kappa = 1/2.

## 5. The plan I would like you to carry out or replace

(a) Sign of the main term. In (P1) the first-order term of S is (A_+ - A_-) Mfun, with A_+ - A_- = -r_Q'(y)/r_Q(y) >= 0 for y > 0 and
        Mfun = int_0^{u_0} s#(u/h) (u^2 T(u))' du,   h = 1/(L r_Q(y)) the local lattice spacing in u.
    Heuristically (Poisson summation with the Lorentzian part of u^2 T) Mfun is about -eta lambda/(e^{2 pi lambda} - 1) <= 0 with lambda = eta/h, up to O(eta h^2) and O(eta h/u_0)-type corrections. So for y > 0 the main term has the favourable sign. Prove a rigorous one-sided bound of the form Mfun <= C eta (h/u_0 + h^2) or better.

(b) Many-node regime: at least M0 cells of the lattice between y and Q (tau_+(delta) >= M0). Then the lattice spacing near y is at most C delta/M0 by (F3). Split w into (0, c delta), (c delta, delta), (delta, delta'). On (0, c delta) use (a) and the second-order remainder, which is at most C c eta/delta. On (c delta, delta) both branches are oscillatory integrals of the sawtooth against T(w) times a density that varies on the scale of the distance to the edge, so a second mean value or summation by parts argument should give C eta/(c^2 delta M0)-type smallness, except within a few cells of the cap, which must be treated with the sign ({tau_+} is nonnegative, and near the cap T is about T(delta)). On (delta, delta') only the left branch remains, and in the form (*) it enters with the favourable sign: int_delta^{delta'} T {tau_-} is at least a fixed fraction of B(delta) - B(2 delta) when there are many left cells in (delta, 2 delta), which unimodality guarantees (the r_Q-mass of [y - 2 delta, y - delta] is at least that of [y, Q]). Conclude (*) for c small and M0 large.

(c) Few-node regime: fewer than M0 cells between y and Q, so at most M0 nodes are concerned. Idea: do not use the standard positions for these nodes. Let A be the set of the top M0 nodes inside (and the K0 nodes beyond Q). Replace their positions by the exact solution of the subsystem {F_j = 0, j in A} with the other nodes frozen at the standard positions (or by any explicit subsystem supersolution). Then the inequalities at j in A hold by construction. The other nodes lose at most 2 pi sum_{k in A} int_{y_k}^{ytilde_k} T(y_j - v) dv where ytilde_k is the new position, so one needs an upper bound on how far the nodes of A move to the right. Diagonal dominance should give it: in counting units a single node's deficit is repaired by shifting that node by a comparable amount when the lattice is locally sparse (own-cell kernel mass close to q), and by at most C(Lambda, M0) times the deficit when it is not (then the node is within C eta of the edge by (F3), and lambda there is at most C Lambda). The frozen nodes at distance D feel at most (number of moved nodes) x (shift in x) x T(D). Show that the many-node margin of (b) absorbs this, possibly with two thresholds M0 << M1. The price is that N_y >= N_0 - K/L holds with K increased by the total counting shift of the moved nodes, which is fine if it is universal.

You may replace any part of this plan. What matters is a correct proof of Theorem G for lambda_Q < Lambda, or a precise statement of the one inequality that remains, with everything else proved.
