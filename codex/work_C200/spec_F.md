# Task F: adversarial referee of a proof (counting bound for XXZ Bethe roots near the ferromagnet at every filling)

You are the referee. Below are definitions, facts to assume, and a text (Sections 4.5 and the first corollary of 4.6) with statements and proofs. The text claims to close a case that an earlier attempt regarded as an obstruction, so be suspicious. Find errors: wrong signs, wrong index bookkeeping (which node has which quantum number, which nodes lie to the left or right), a hidden dependence of a "universal" constant on eta, Q, L, n or on each other in a circular way (the order of choice is A, then M, then k), missing cases (n <= k, Q small, the last cell longer than Q/2, the node at Q, nodes with y <= 0 but Q - y < 1, Case 1 versus Case 2 boundary), illegitimate use of monotonicity of r, and any step that does not follow. Re-derive each identity and inequality yourself, in particular: the residual identity with nodes at infinity, inequality (4.2), the identity L p(z)/(2 pi) = E_i + (1 - q) n/2 + int B(z - v) dmu(v), the constant terms summing to q I_j, the right-endpoint comparison, the logarithmic lower bound on E_i, and the absorption in (d) and in Case 2. Also test the construction against simple sanity scenarios (for example n = k + 1, or a very dilute sea with ell comparable to Q) and say whether the proof covers them. For each lettered step say "correct" or give the precise defect and, if you can, the repair. Finish with a verdict per item: correct as stated, correct after listed repairs, or wrong. Tags such as [rigorous], (Astra), (Fable), (checked), [numerical] are bookkeeping and can be ignored; numerical remarks are not part of the proof.

## Definitions and facts to assume

eta in (0, 1/64], gamma = pi - eta, q = 1 - 2 eta/pi.

    p(x)     = 2 atan(tan(eta/2) tanh x),   strictly increasing, odd, p(+-inf) = +-eta
    theta(x) = -2 atan(tanh x/tan eta),     odd, strictly decreasing, theta(+-inf) = -+(pi - 2 eta)
    a1       = p'/(2 pi) = sin(eta)/(pi (cosh 2x + cos eta)) > 0
    T        = -theta'/(2 pi) = sin(2 eta)/(2 pi (sinh^2 x + sin^2 eta)) > 0, even, decreasing in |x|, int T = q, T(0) = cot(eta)/pi
    Theta(x) = int_0^x T (odd, Theta(inf) = q/2),   B(d) = int_d^inf T = q/2 - Theta(d) for d >= 0.

Integers L >= 2, 1 <= n < L/2, nu = n/L, I_j = j - (n+1)/2. Bethe map F_j(x) = L p(x_j) - sum_{k=1}^n theta(x_j - x_k) - 2 pi I_j on R^n.
Fact 1: F has a unique real zero x* (increasing), symmetric under x -> -x. N_x(v) = (1/L) #{j : x*_j <= v}; sigma_L = a1 + (1/L) sum_j T(. - x*_j).
Fact 2 (continuum): Q = Q_nu in (0, inf) and a continuous r = r_Q on [-Q, Q] with r = a1 + int_{-Q}^{Q} T(. - v) r(v) dv on [-Q, Q], int_{-Q}^{Q} r = nu; r > 0, even, nonincreasing on [0, Q]. s(u) := a1(u) + int_{-Q}^{Q} T(u - v) r(v) dv on R (s = r on [-Q, Q]). N_0 = cumulative function of r 1_[-Q,Q]. t(v) := L int_{-Q}^{v} s (so t = L N_0 on [-Q, Q], t(-Q) = 0, t(Q) = n, t(0) = n/2, t' = L s > 0 on R).
Fact 3 (Lemma 1.1, comparison with nodes at +infinity): if a trial configuration has finite nodes y_j (j in Fin) and the other n_inf nodes at +infinity, and Ft_j := L p(y_j) - sum_{k in Fin} theta(y_j - y_k) - n_inf (pi - 2 eta) - 2 pi I_j >= 0 for every j in Fin, then x*_j <= y_j for every j in Fin. (The indices of the finite nodes are the lowest ones, in increasing order of position.)
Fact 4 (Theorem 4.3): with delta(y) = Q - |y|, rho(y) = min(1, eta + delta(y)), sawtooth s#(u) = u - floor(u) - 1/2, and S(y) := int_{-Q}^{Q} T(y - v) s#(t(v) - t(y)) dv, one has |S(y)| <= C_S eta/rho(y) for all y in [-Q, Q], with C_S universal (independent of eta, Q, L, n, y).
Section 4.1 of the earlier text (assume): the trial with interior nodes t(y_j) = j + K0; S as above.

## The text to referee

## 4.5 Every filling: exterior nodes on the scale of the last cell

Theorem 4.5 used the trial nodes beyond Q at the quantile positions of t, which are far from Q when the edge is sparse, so that the right exterior contributes only B(Q - y)/2. The managing agent took this for an obstruction (spec_D.md, Section 4) on the strength of a miscount in the top node's own inequality. It is not one: the exterior nodes may be placed at explicit positions on the scale of the last continuum cell, and they are supersolution components there because the continuum counting function grows logarithmically beyond Q.

[rigorous] (Astra, astra_D.out Section 1; checked) Lemma 4.10 (arbitrary exterior positions). Let k = K0 >= 0, kappa = k + 1/2, n > k. Keep the interior nodes t(y_j) = j + k, 1 <= j <= n - k, place m <= k nodes at points Q < z_1 < ... < z_m, and the other k - m nodes at +infinity. Then at every interior node y = y_j,

    Ft_j/(2 pi) >= (1 - q) kappa + kappa B(Q + y) + (1/2) B(Q - y) + sum_{i=1}^m B(z_i - y) - S(y).              (4.2)

Proof. The residual identity of 4.1 holds for any finite positions: Ft_j/(2 pi) = t(y_j) - n/2 - I_j - int T(y_j - v) D(v) dv with D = L N_0 - #{finite nodes <= v} (it follows from L p(u)/(2 pi) + int Theta(u - v) L r(v) 1_[-Q,Q] dv = t(u) - n/2, both sides having derivative L s and vanishing at 0, and Theta(-inf) = -q/2 for the nodes at infinity). At an interior node the first three terms give kappa. D = 0 below -Q, D <= kappa + s#(t(v)) on [-Q, Q], and D = k - sum_i 1{v >= z_i} above Q, with int_{v > Q} T(y - v) 1{v >= z_i} dv = B(z_i - y). Substituting, and using that t(y_j) is an integer, gives (4.2). End.

So each exterior node at a finite point contributes its full tail B(z_i - y).

[rigorous] (Astra and checked) Kernel bounds, with universal b, c > 0:

    b eta/x <= B(x) <= eta/(pi x)  (eta <= x <= 2; the upper bound for all x > 0),      (1 - q) + B(x) >= c eta/min(1, eta + x)  (x >= 0).      (4.3)

(T(w) <= eta/(pi w^2) for all w, T(w) >= c_1 eta/w^2 for eta <= w <= 4, integrate over [x, 2x]; for x < eta use B(x) >= B(eta) >= b, for x >= 1 use 1 - q = 2 eta/pi.)

[rigorous] (Astra, astra_D.out Sections 2 to 5; checked line by line; numerically in ext_trial.py) Theorem 4.11 (counting bound at every filling). There is a universal K such that for every eta in (0, 1/64], L >= 2 and 1 <= n < L/2,

    |N_x(v) - N_0(v)| <= K/L   for all real v.

Proof. Universal constants A, M, k are chosen in this order, with k >= M and k >= 2; C_S is the constant of Theorem 4.3. If n <= k the bound is trivial with K = k. Let n > k, so n >= 3 and the last continuum cell lies in the right half: define ell by L int_{Q - ell}^{Q} r = 1, and d = max(eta, ell), H = A d, R = min(1, Q).

Case 1: M H <= R/4. Use m = M exterior nodes at z_i = Q + i H, and k - M nodes at infinity.

(a) Density at depth. r is nonincreasing on [0, Q], so r(Q - ell) is at least the average of r over the last cell, 1/(L ell), and L r(Q - u) >= 1/ell >= 1/d for ell <= u <= Q.

(b) Logarithmic growth of the counting function beyond Q. Fix i, x = i H <= R/4 <= Q/4, E_i = t(Q + x) - n = L int_0^x s(Q + w) dw. Since a1 > 0 and s(Q + w) >= int_d^x T(u + w) r(Q - u) du, and for d <= u <= x, 0 <= w <= u one has eta <= u + w <= 2u <= 1/2 and T(u + w) >= c_1 eta/(4 u^2),

    E_i >= int_d^x L r(Q - u) int_0^u T(u + w) dw du >= (c_1 eta/(4 d)) int_d^x du/u = c_2 (eta/d) log(A i).

(c) The exterior nodes are supersolution components. The node z = z_i has index j = n - k + i. Writing Theta(w) = q/2 - B(w) at positive separations, the constant parts of the interaction sum to q I_j and

    Ft_j/(2 pi) = L p(z)/(2 pi) - (1 - q) I_j - sum_{a < j} B(z - y_a) + sum_{a > j, finite} B(y_a - z).

By the identity in the proof of Lemma 4.10 at u = z > Q, L p(z)/(2 pi) = E_i + (1 - q) n/2 + int B(z - v) L r(v) 1_[-Q,Q] dv. The function v -> B(z - v) is increasing, so on each of the n unit-mass cells of [-Q, Q] its value at the right endpoint exceeds its cell integral by at most its increment over the cell; summing, the sum over all n right endpoints exceeds the integral by at most B(z - Q) = B(i H). The interior nodes are the right endpoints of the cells k + 1, ..., n, so their sum is no larger. Dropping the nonnegative terms from nodes to the right, and using n/2 - I_j = k - i + 1/2 > 0 and B(l H) <= eta/(pi l A d),

    Ft_j/(2 pi) >= E_i - B(i H) - sum_{l=1}^{i-1} B(l H) >= (eta/d) [c_2 log(A i) - (1 + log i)/(pi A)].

Choose A >= e with 1/(pi A) <= c_2. The bracket is at least c_2 (1 + log i) - c_2 (1 + log i) = 0 for every i >= 1, independently of M.

(d) Interior nodes. Let 0 < y < Q be an interior node, delta = Q - y. Since t(Q) - t(y) is a positive integer, delta >= ell. If delta < 1, then rho(y) = min(1, eta + delta) >= delta and rho(y) >= d (d <= 1 by the case assumption), and eta <= z_i - y = delta + i A d <= 5/4, so by (4.3)

    sum_{i=1}^M B(z_i - y) >= b eta sum_{i=1}^M 1/(delta + i A d) >= (b eta/rho(y)) sum_{i=1}^M 1/(1 + i A).

The harmonic sum diverges, so M can be chosen (A being fixed) with b sum_{i <= M} 1/(1 + i A) >= C_S. Then S(y) <= C_S eta/rho(y) is covered in (4.2). If delta >= 1, rho(y) = 1 and (1 - q) kappa >= C_S eta for kappa >= pi C_S/2. If y <= 0, rho(y) = min(1, eta + Q + y) and kappa [(1 - q) + B(Q + y)] >= kappa c eta/rho(y) >= C_S eta/rho(y) for kappa >= C_S/c, by (4.3). At y = Q, S(Q) <= (1/2) int_{-Q}^{Q} T(Q - v) dv <= B(0)/2, which is the third term of (4.2).

Case 2: R < 4 A M d. Put all k exterior nodes at infinity (m = 0). The arguments for y <= 0, for delta >= 1 and for y = Q are unchanged. For 0 < y < Q with delta < 1, again delta >= ell, so rho(y) >= min(1, d), while min(1, eta + Q + y) <= min(1, eta + 2R) <= (1 + 8 A M) min(1, d) <= (1 + 8 A M) rho(y). By (4.3), kappa [(1 - q) + B(Q + y)] >= (c kappa/(1 + 8 A M)) eta/rho(y), which covers C_S eta/rho(y) once c kappa >= (1 + 8 A M) C_S. (When the Fermi sea is short on the scale of its last cell, the left exterior alone suffices.)

Choose k >= M with all the stated lower bounds on kappa. Then Ft_j >= 0 at every finite trial node, Lemma 1.1 gives x*_j <= y_j there, and N_x >= N_y >= N_0 - (k + 1)/L (on [-Q, Q] the trial count is max(0, floor(t) - k), above Q it is at least n - k). Reflection, with left limits at the jumps, gives the upper bound. K = k + 1. End.

[numerical] (Fable, ext_trial.py) With k = M = 4 and A = 2, exterior nodes at Q + 2 i d, i = 1..4, for 1 + Delta in {1e-3, 1e-4, 1e-5}, L in {64, 256, 1024} and five fillings (lambda_Q from 1.2 down to 0.0001, last cell from 0.7 eta to 240 eta): min over exterior nodes of Ft/(2 pi eta) is between 0.945 and 46.6, and min over interior nodes between 2.80 and 4.43. So the construction is a supersolution in every sampled case, with small explicit constants.

Together with Theorem H (K = 1, even L, n = L/2), the counting bound holds for every L and every 1 <= n <= L/2.

## 4.6 Consequences: logarithm-free (D) at every filling, and the strip

[rigorous] (Fable; this is C149 Theorem 3.3 with kappa0 + 1/2 replaced by K, re-derived) Corollary 4.12. For eta in (0, 1/64], every L and 1 <= n <= L/2, and every real u,

    eta L |sigma_L(u) - s(u)| <= 2 K eta cot(eta)/pi <= 2K/pi,      eta^2 L |sigma_L'(u) - s'(u)| <= 4 K eta^2 T_A(0) T_B(0) <= 16 K/pi^2,

where s = r_{Q_nu} on [-Q_nu, Q_nu]. Proof: sigma_L - s = T*(mu_x - mu_0) = int T'(u - v)(N_x - N_0)(v) dv and int |T'| = 2 T(0). For the derivative, T = T_A * T_B with T_A, T_B positive, even and unimodal (C149 Lemma 3.2: the Fourier transform sinh((pi - 2 eta) w/2)/sinh(pi w/2) factors as [sinh((pi - 2 eta) w/2)/sinh(gamma w/2)] [sinh(gamma w/2)/sinh(pi w/2)], each factor the transform of a kernel of the form sin(c)/(2 beta (cosh(pi x/beta) - cos c)); checked), so sigma_L' - s' = T_A' * [T_B * (mu_x - mu_0)], and each convolution costs twice the central value. This is the logarithm-free (D) and (D1) of C138, relative to the continuum Fermi point Q_nu, with sigma_L >= (1 - 2K/(pi lambda)) r_{Q_nu} wherever lambda = eta L r_{Q_nu} > 2K/pi.

