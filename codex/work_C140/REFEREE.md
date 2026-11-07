# Referee report on C139 Theorem A + C140 Theorem 4.1 (repulsive U(1) class, 0 <= Delta <= 63/64)

2026-11-07. Independent referee pass by a Claude Opus agent. Read-only on C139/C140 except this file and the two scripts `referee_thm41.py`, `referee_ed.py` (outputs `referee_thm41.json`, `referee_ed.json`). J=1.

## Verdict

The headline claim is correct as stated, up to minor presentational repairs. That claim is |F_L(n) - L e_bulk| <= (1086/delta)/L, and hence L(E_2L - 2E_L) <= 345600 M on 0 <= J_z/J <= 63/64 for every field, DM angle and L >= 2. I found no fatal error and no gap in the logic of C140 Theorem 4.1, and none in the C139 replacements I rechecked line by line. The secondary claim (19935 + pi^2/4) M on the same range is not established to the same standard. It is a degree-counting rescaling that the authors say they did not redo, yet RESULT.md lists it under a [Proved] heading (issue 5 below). It is plausible, and it does not affect 345600.

## 1. C140 Theorem 4.1, line by line

(F1) Dual identity. Correct. sigma + K*sigma = a1 + K*(sigma dx - mu), and (I+K*)^(-1) = I - J* on bounded integrable functions, with (I+K*)^(-1)a1 = r0. I checked the factorization of Jhat = Khat/(1+Khat) into the k1 and r0 transforms directly (sum-to-product on 1 + Khat). It is right, so J >= 0 is even and unimodal, with ||J||_1 = (pi-2gamma)/(2(pi-gamma)).

(F2), (F3). Correct. z_L(x_(j+1)) - z_L(x_j) = 1/L. With p(inf) = pi - gamma and theta(inf) = pi - 2gamma, c = (pi-gamma)/(2pi) - nu(pi-2gamma)/(2pi), 2c - nu = (1-gamma/pi)(1-2nu), and z_L(x_n) = nu/2 - 1/(2L). Numerically the tail mass matches c - nu/2 + 1/(2L) to 5e-8 (quadrature error) and the gap residuals are 4e-15.

Section 2, kernel constants. Correct. sinh(Aw)/sinh(Cw) <= e^(-(C-A)w) reduces to C >= A, with C - A = gamma, and u/(1+u) is increasing. The identity sum_(m>=2) t^(m-1)(1+t)^(-m) = t/(1+t) checks. The step J^(*m)(0) = sup J^(*m) uses Jhat >= 0 (positive definiteness), which holds but should be stated. Numerically J(0) = 1.2412 against the bound 1.2465 at Delta = 63/64, and sum_(m>=2) J^(*m)(0) = 0.538 against 0.552.

Lemma 3.1, the two-neighbour lemma. Correct. For x in g_k, J(x-y) is nondecreasing in y on (-inf, x] and nonincreasing on [x, inf). So pairing g_j with x_j (j < k) and with x_(j+1) (j > k) gives nonnegative pairs under (P), and g_k with its two endpoints contributes at least min/L - (min+max)/L = -max/L. In (b) the two pairings give the stated lower and upper bounds, and the mirror half of out is inside J_out v.

Corollary 3.2 and (3.1), (3.2). Correct. ||J_out|| <= q_J < 1 on L^inf and on L^1, the resolvent is positive, and v is bounded, so the inequalities propagate through (I - J_out)^(-1). ||psi||_1 = q_J, so ||phit||_1 <= q_J/(1-q_J) <= 1. The majorization J_out^k psi <= J^(*(k+1))(.-Q) + J^(*(k+1))(.+Q) gives (3.2). Numerically the upper comparison sigma <= s_(x_n) beyond x_n holds in every case (max violation -2e-11, i.e. none).

Lemma 3.3, outermost-root location. Correct. 1/(2L) <= int_(x_n)^inf sigma <= int_Q^inf s_Q = h_Q/2 <= int_Q^inf r0/(1-q_J) <= 2 r0(Q)/(b(1-q_J)), with 1 - q_J = pi/(2(pi-gamma)). Numerically L r0(x_n) exceeds the proved bound by a factor 2.07 to 3.34, and 2L int_Q^inf s_Q >= 1.09 against the needed 1.

Theorem 3.4. Correct. The loss is (ln2 + 2(1-ln2))/(pi gamma L) = (2-ln2)/(pi gamma L), and the conversion to theta r0 uses r0(x) >= r0(x_n) on [x_1, x_n] (symmetry) and Lemma 3.3. Numerically sigma - [r0 - (2-ln2)/(pi gamma L)] > 0 in every case, and min sigma/r0 on [x_1,x_n] is 0.974 at Delta = 63/64, against theta = 3.9e-4.

The 63/64 cutoff. Correct, but tight. Positivity of theta requires gamma > gamma_c = pi - pi^3/(8(2-ln2)) = 0.1758535, i.e. Delta < 0.9845776. arccos(63/64) = 0.1770077, so theta(63/64) = 3.89e-4 > 0. The continuity argument only needs theta > 0, so the small margin is harmless.

Section 4, continuity argument. Valid. It does not need the ground state to stay on the branch a priori, and a change of branch cannot occur. The branch is defined by continuation of a nondegenerate critical point, and the identification is made afterwards by analyticity (below). Openness of G follows from the implicit function theorem with the positive definite Hessian 2pi L diag(sigma(x_j)) + Theta (I rechecked the Hessian formula, including the diagonal theta'(0) and the Bochner positivity of Theta). Closedness uses the a priori ball of C139 Section 2.1 (I rechecked the recession bound (pi-gamma)(L-2n+2) sum|x_j| - C), distinctness of coordinates, and non-strict (P) in the limit, which is all Lemmas 3.1 and 3.3 need. Theorem 3.4 then upgrades (P) to strict. Reflection symmetry has to hold along the branch for Theorem 3.4 (Q = x_n = -x_1). It does, because the reflected branch solves the same equations and the coincidence set is open and closed. The n = 1 case checks: a1/K lies between 1/(2cos gamma) and (1+cos gamma)/cos gamma, both > 1/2 >= 1/L.

Numerical continuation from gamma = pi/2 (seed x_j = artanh tan(pi I_j/L)) to arccos(63/64), in 400 steps, for (L,n) = (7,3), (8,4), (12,6), (13,6), (33,16), (64,32). It lands on the authors' direct solution to 2e-15, with L sigma(x_j) >= 0.97 along the path.

## 2. Does Theorem 4.1 give (I) at both parities?

Yes. (N) follows from strict (P) at the roots through the Hessian. The identification then runs through C139 Section 2.3 and the C134 argument, which needs no Yang-Yang convexity. I checked its ingredients for this setting. The sector Hamiltonian at flux pi per bond (total L pi) has off-diagonal entries -J/2 in the S^z basis (check_c129.matrix), so it is stoquastic and irreducible at every Delta, and F_L(n)(gamma) is simple and real-analytic. At gamma = pi/2 the Bethe vector is the free-fermion Slater determinant, which is nonzero and is the ground state. Hence det(H_n(gamma) - E(gamma)) vanishes identically along the analytic branch, E = F near pi/2, and E = F on [gamma_0, pi/2] by the identity theorem. This holds for odd L as well, since the centered I_j encode flux L pi at either parity.

The ED cross-check confirms all of this (Section 5).

## 3. C139 sign audit, spot checks

- Lemma 4.1, ||D|| <= B/(2L) from ordering alone. Correct. z_L - z_nu = -K*D by parts, and z_nu clipped equals Cum(1_Q r) - nu/2 because z_nu(-Q_nu) = -c + h_Q/2 = -nu/2 (h_Q = 2c - nu). At roots D(x_j^+) = 1/(2L) + P(z_L(x_j)) - P(z_nu(x_j)), D is nonincreasing between roots, and D(+-inf) = 0.
- Lemma 4.3, the clipping step. Correct. |Cum R| <= (1+2q)B/(2L) <= 3B/(2L), and the mismatch bound uses Var(Z) <= q_J. I rechecked that bound by the double integral 2 int_0^inf u|J'| = q_J.
- Restricted inverse <= 3. Correct. F = g1_Q - K*(1_Q f) equals f on [-Q,Q] and satisfies F = (I-J*)(g1_Q) + J*(1_out F), so ||F|| <= (1+q_J)/(1-q_J) ||g||. Invertibility comes from L^2 positivity of K_Q plus the same inequality with g = 0.
- r0 <= r_Q <= (1+4/pi) r0. Correct. int_out r0 = (2/pi) atan(e^(-bQ)), U1 <= J(0) h_Q <= (1/(2gamma))(4/pi)e^(-bQ), and r0 >= e^(-b|x|)/(2gamma).
- Quadrature (5.2) and ends (5.3). Correct. These cover the trapezoid Peano bound Var(f')/(8L), |(r/sigma)'| <= 4(bD'+D1')/(Lr), int dx/(L r0) <= 1/(bD'), the end identity -L int eps' Cum R (using Cum R(x_j2^+) = 1/(2L) and eps(Q_L) = 0), and int_(x_j2)^(Q_L) r0 <= [3B/2 + 4D'/b]/L. The printed bound has an extra factor c_r, which is slack.
- epsilon_Q' >= 0 by the odd reduction J(x-y) - J(x+y) >= 0. Correct, and the boundary terms vanish because eps(+-Q) = 0.

## 4. Constants

345600 = (5400/delta) M with delta = 1/64, from C139 Section 6. Using constants.py on a 20001-point grid of [arccos(63/64), pi/2], sup K_rep = 65190 <= 1086*64 = 69504, and sup max(4 L_cut, 5 K_large) = 325952 <= 5373*64 = 343872. Adding pi^2/4 gives 325954 <= 345600. The majorants K_large <= 47.1 + 2055.3/gamma^2 and K_rep <= 2171.5/gamma^2 hold pointwise, and gamma_0 = 0.17701 >= sqrt(2/64) = 0.17678. The constant follows.

Sawtooth bookkeeping (C140 5.2). The bound sup|Cum R| <= 3/(2L) is correct, from z_L increasing on [x_1,x_n] and (3.1). For the rescaling of C139 8.2 I checked the degrees. D, D1, ||G||, |nu-m| and L_cut are linear in the sawtooth constant, the ends and curvature terms are quadratic, and V_f is invariant. A factor 9 is therefore a valid (loose) majorant, but only if every step of 8.2 survives with 1/2 replaced by 3/2. I did not redo that, and neither did the authors.

## 5. ED cross-check (referee_ed.py)

- Bethe vs ED at flux pi per bond, L = 2..12, all 1 <= n <= L/2, Delta in {0.25, 0.5, 0.9, 0.98}, 36 sectors per Delta, both parities. The maximum error is 4.8e-14.
- Flux pi per bond is the sector minimizer. Over 11 other fluxes pi + phi/L (L <= 11), every ED energy is higher, by at least 1.1e-3.
- L|F_L(n) - L e_bulk| / K_rep(gamma) <= 1.0e-3 (largest at Delta = 0.25).
- End-to-end L(E_2L - 2E_L)/M for L = 2..6, DM in {0, 0.3, 1.0}, fields 0 to 2.2J (23 values), Delta in {0.25, 0.5, 0.9, 0.98}. The maximum is 3.97 (Delta = 0.98, L = 2), far below 345600 and consistent with the authors' supremum 4 at Delta = 1.
- Theorem 4.1 numerics (referee_thm41.py), Delta in {0.25, 0.5, 0.9, 0.98, 63/64}, L = 2..12 and 16 to 257, up to 6 fillings per L (2 <= n <= L/2), both parities. Every lemma inequality holds, the minimum Hessian eigenvalue is 6.8, and no negative mass beyond x_n was found.

## 6. Issues

1. C140 Section 4, closedness of G. Minor. "By local uniqueness the branch reaches gamma_*" is terse. The limit is taken along a subsequence, and one needs the open-closed argument on (gamma_*, gamma_*+eps) to show that the existing branch coincides with the IFT branch through the limit point. Repair: add two sentences.
2. C140 Theorem 4.1, identification. Minor. It cites C139 Section 2.3, which is itself a one-paragraph reference to C134. Repair: state once the facts used (stoquastic and irreducible at flux L pi with entries -J/2, nonzero Slater determinant at gamma = pi/2, identity theorem for det(H_n - E(gamma))), and note that they hold at odd L.
3. C140 Section 2. Minor. J^(*m)(0) = sup J^(*m) needs Jhat >= 0. Say so.
4. C140 Section 3. Minor. The cutoff should be recorded as Delta < 0.98458 (gamma > 0.175853), with theta(63/64) = 3.9e-4. The text says "theta = 0.008 (0.98)" but does not give the 63/64 value. The bracket "proved value is 0.125 to 0.235" should read up to 0.25 at Delta = 0. The symbol theta(gamma) clashes with C139's phase theta(x). Rename it.
5. C140 RESULT bottom line and 5.2, the constant 19935 + pi^2/4. Gap, fixable. It appears under [Proved] although the 8.2 rescaling was "not redone line by line". If correct it supersedes 345600 on the same range, so the headline is internally inconsistent. Repair: redo C139 8.2 with sawtooth 3/(2L), or retag it as [Proved, bookkeeping sketch] and keep 345600 as the proved headline.
6. C140 5.3 and RESULT, (L*). Minor, outside the refereed claim. (L*) must be strict (at lambda_* the margin theta is zero and openness fails). The phrase "(L*) then propagates" is circular as written. Repair: assume L r0(x_n) > lambda_* b for every symmetric ordered critical point satisfying (P), and run the same continuity argument.
7. C140 RESULT. Minor. "ED at L <= 12 not rerun" is now covered by this report at four new Delta values.

No issue is fatal.
