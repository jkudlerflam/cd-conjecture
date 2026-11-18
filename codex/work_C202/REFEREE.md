This referee report was written by an independent Claude agent (Fable 5.2, 2026-11-18), not by Codex. One adversarial line read of the proof was delegated to Astra (model gpt-6.2-astra, effort high, transcript 20261118T110136-57918), and I checked each of its findings myself. No call to gpt-6-astra was made.

# C202 referee report on C148 (the easy-axis side, J_z/J > 1)

J = 1 throughout, Delta = cosh(lambda), M_0 = (1 + Delta)/4. Read: work_C148 (PROOF, RESULT, pos_interval.py, pos_rigorous.py, pos_check.py, pos_window.py, constants.py and their outputs), REPORT Section 9.106, work_C147/REFEREE.md, work_C142/PROOF.md Part I (for the lemmas C148 cites as "verbatim"), and work_C129 (RESULT, and PROOF Sections 3 to 6 for the reductions). I did not reuse any solver from those folders. Every script in work_C202 was written from scratch (ref148.py holds the ED, the circle Bethe solver and the restricted continuum).

## Verdicts in one place

| item | verdict |
| --- | --- |
| 1. Positivity of sigma_L (Lemmas 4.1 to 4.3, Theorem 4.4, the interval certificate, the lambda >= 1.36 bound) | correct after minor repairs (one unstated subcase in Lemma 4.1(a), two decimals rounded the wrong way) |
| 2. Ground-state identification at both parities and all fillings (Theorem 4.5) | correct after a minor repair (the roots are analytic through Delta = 1 in u, not in x). Completeness of the Bethe ansatz is not used |
| 3. Uniformity of K in the Fermi point, half filling and the dilute regime near h_c1 (Sections 3, 5, 6) | correct after a minor repair (one case of the end estimate is not written, and it is realized at L = 4, n = 2) |
| 4. Constant chain, (5/2)K <= 10778 M_0, plus pi^2/4, and M >= M_0 | correct after a minor repair (displayed intermediate decimals do not sum to 218.4, the unrounded majorant is 218.381 and gives 4310.68 <= 4311) |
| 5. Independent numerics | no violation of any inequality, all of C148's quoted numbers that I recomputed are reproduced |

I found no error that changes the theorem or the constant 10781. Every problem below has a short repair that stays inside the stated budget.

## 1. Positivity of the counting density

### 1.1 The pairing argument on the circle (Lemma 4.1, Corollary 4.2)

I rederived sigma_L = r0 + J (*) (sigma_L dx) - J (*) mu from (I - J(*)) K(*) = J(*) (multipliers (1 - t/(1+t)) t = t/(1+t)), and went through the pairing gap by gap.

(a), case 1 (|x| + Q <= pi/2). Circular distance from x equals |x - y| on all of in, so J(x - .) is unimodal about x there. Gaps to the left pay for their left endpoints, gaps to the right for their right endpoints, and the gap containing x pays for the smaller of its two endpoint values. One root term of size at most J(0)/L is left. Correct.

(a), case 2 (x > 0, x + Q > pi/2). With x* = x - pi/2, the distance from x is pi - x + y on [-Q, x*] (increasing in y) and x - y on [x*, x]. If x* is in g_m and x in g_k with m < k, gaps i < m pay for roots 2..m, g_m is dropped, gaps m < i < k pay for roots m+1..k-1, g_k pays for the smaller of roots k, k+1, and gaps i > k pay for roots k+2..n. Root 1 and one term of size at most J(0)/L are left, and no root is used twice. For y in out, 0 < y - x <= pi - Q - x = d(x, x_1) < pi/2, so J(x - y) >= J(x - x_1). Correct.

Minor repair. If x and x* lie in the same gap (m = k, which needs a gap longer than pi/2, for instance n = 2 with Q > pi/4), the text both drops that gap and pairs it with two endpoints. Drop it. Then gaps i < m pay for roots 2..m and gaps i > m for roots m+2..n, and the two unpaid roots are x_1 and x_(m+1), with J(x - x_(m+1)) <= J(0). The stated inequality is unchanged.

(b). For y in out with y* in a gap g_m, J(y - .) decreases toward y* from both sides. Farther endpoints (lower bound) cover roots 2..n-1 once and leave x_1 and x_n, which is psi/L. Nearer endpoints (upper bound), with g_m taking its endpoint of larger J, leave one root with a minus sign. If y* lies in out, J(y - .) is monotone on all of in and the same assignment leaves one end root. Correct. Corollary 4.2 follows from positivity of the Neumann series of J_out (norm at most omega_Q <= 1/2 on L^inf and L^1), with ||psi||_1 <= 1, phit even, and sum_(m>=2) J^((*)m)(0) = K(0) - J(0). Correct. It does not need v >= 0 on out.

### 1.2 The outermost root (Lemma 4.3)

T(Q) = int_Q^(pi/2) r0 = (pi/2 - am(u_Q))/(2 pi) (the substitution u = 2 K_e x/pi gives r0 dx = dn du/(2 pi)). From (F3) and v <= s_Q, 1/L <= h_Q <= 2 T(Q)/(1 - omega_Q), and L >= 2n >= 4. In units of K_e/pi^2, with phi = am(u_Q), the bound is dn(phi) max(4, pi/(pi - 2 phi)). For phi <= 3 pi/8, 4 dn >= 4 cos(3 pi/8) = 4 sin(pi/8). For psi_a = pi/2 - phi <= pi/8, dn >= k sin(psi_a) and sin(psi_a)/psi_a >= sin(pi/8)/(pi/8). Correct. I checked k K_e = (pi/2) theta_2(q)^2 and theta_2(e^(-lambda)) >= sqrt(pi/lambda)(1 - 2 e^(-pi^2/lambda)) in 40 digits at five values of lambda (check_constants.json). The second is the alternating series of the Jacobi imaginary transformation.

### 1.3 Theorem 4.4 and the certificate

Both cases of Lemma 4.1(a) give sigma_L(x) >= r0(x_n) - (2K(0) - J(0))/L on [x_1, x_n] (in case 2 the out term is at least J(x - x_1) h_Q with h_Q >= 1/L, which pays for the left-over root). With K(0) = coth(lambda)/pi this is the stated

    Phi(lambda) <= lambda [2 coth(lambda) - pi J(0)] / [2 pi sin(pi/8) (1 - 2 e^(-pi^2/lambda))^2].

pi J(0) >= ln2/lambda: the positivity of the Poisson images that the proof invokes holds (the line kernel is a convolution of a Lorentzian and a sech). Astra gave a shorter argument that avoids it: f(u) = 1/(1 + e^(2 lambda u)) is convex on u >= 0, so its trapezoid sum dominates its integral. Either is fine.

Monotonicity directions, each checked by hand and on a grid.
- On (0, 1]: lambda coth(lambda) is increasing, 1 - 2e^(-pi^2/lambda) is positive and decreasing, so the majorant is increasing. Its value at lambda = 1 is 0.80405, and its limit at 0 is 0.54351.
- On a cell [l0, l1] of [1, 1.36]: lambda <= l1, coth(lambda) <= coth(l0) (coth decreasing), 1/2 + 2 sum_(j<=3) t^j/(1 + t^j) is decreasing in lambda (t = e^(-2 lambda) decreasing, t^j/(1+t^j) increasing in t), so its value at l1 is a lower bound, and the denominator at l1 is a lower bound. The bracket 2 coth(l0) - (...)(l1) is positive on every cell, which is needed to multiply the upper bounds. All four directions are the ones the script uses.

Reruns. pos_interval.py gives worst cell 0.938945433861 < 1, direct margin 0.020999136, f_small(1) = 0.804053542. pos_rigorous.py gives the same in floating point. My own interval version (check_constants.py, mpmath iv at 40 digits) gives 0.9389454 with 36 cells, 0.9498 with 9 cells and 0.9355 with 3600 cells. As a check that does not depend on any monotonicity direction, I also evaluated the expression with each whole cell inserted as an interval. The 36-cell upper enclosure is then 0.9636 < 1. So the certificate is sound.

The lambda >= 1.36 claim. min a1 = a1(pi/2) = tanh(lambda/2)/pi and max K = K(0) = coth(lambda)/pi, nu <= 1/2, so sigma_L >= [tanh(lambda/2) - coth(lambda)/2]/pi for any n <= L/2 points, and the bracket is increasing. It vanishes at lambda = 1.3170 and equals 0.0209991 at 1.36. Correct. The two windows overlap: the pairing majorant with the exact pi J(0) is 0.9355 at 1.36 and passes 1 only near 1.48 (it is 1.0097 at 1.5), and the direct bound is positive from 1.317.

Minor repairs. (i) "0.021/pi" should be 0.02099/pi, since 0.0209991 < 0.021. (ii) "Phi <= 0.804 for lambda <= 1" should be 0.805, since the value is 0.80405. (iii) The same-gap subcase of Lemma 4.1(a) above. (iv) Lemma 2.2(i) reuses the symbol q for ||K||_1 = 1 while q = e^(-lambda) is the nome.

n = 1: min a1/K = Delta/(1 + Delta) > 1/2 >= 1/L (the ratio (cosh 2 lambda - c)/(cosh lambda - c) is increasing in c = cos 2x). Correct.

Numerically (scan_bethe.jsonl, 3172 rings with n >= 2): min of sigma_L/r0(x_n) over [x_1, x_n] is 0.974 (proved 1 - Phi >= 0.061 for lambda <= 1.36), L r0(x_n) is at least 1.81 times the Lemma 4.3 bound, sigma_L is positive on the out arc in every case (the minimum is 4.5e-10 at Delta = 1.0001, where r0 itself is exponentially small near pi/2, and the negative mass is zero in every case), and the smallest Hessian eigenvalue is at least 2 pi L r0(x_n) for every ring with L <= 40 at six values of Delta.

## 2. Ground-state identification

The logic of Theorem 4.5 is: (1) C142 Theorem I.4 gives the branch on 0 < gamma <= pi/2 with sigma >= 0.168 r0 on [x_1, x_n], cosh(pi u_n) <= L, and E = F_L(n) at both parities. (2) In XXX rapidities u the equations and the energy are analytic in s = gamma^2 = -lambda^2 near s = 0, a subsequential limit u* at s = 0 has distinct coordinates and positive definite Hessian, and the implicit function theorem with the open-closed coincidence argument continues the branch to small s < 0. (3) On lambda in (0, infinity) the set where the branch exists with (P) strict and positive definite Hessian is open, and closed by compactness of T, (F1) (no ordinary collision), (F4) (no collision across +-pi/2, since symmetry would force x_1 + pi = x_n, which needs n = (L+1)/2) and Theorem 4.4 (non-strict (P) at the limit becomes strict). (4) E(Delta) and F_L(n)(Delta) are real-analytic and agree on (1 - eps, 1), so they agree for all Delta > 1.

I checked each step. The Jacobian of the Bethe map is 2 pi L diag(sigma_L(x_j)) + Theta with Theta_jk = theta'(x_j - x_k) positive semidefinite (Khat_k = t^|k| > 0), and sigma_L(x_j) includes the k = j term, as it should. F_L(n) is a simple eigenvalue for every real Delta (off-diagonal entries -1/2, or -1 at L = 2 where the two oriented bonds add, irreducible in the sector), so it is analytic on R. Step (4) uses only that E is an analytic function equal to F_L(n) on an interval. Neither completeness of the Bethe ansatz nor nonvanishing of the Bethe vector on Delta > 1 is used, and the proof never needs to know that E(Delta) is an eigenvalue there. Parity enters only through C142 (flux pi per bond, refereed in C147) and (F4). There is no missing case at n = L/2 or n = (L-1)/2.

Minor repair. Theorem 4.5 as stated asserts a real-analytic family Delta -> x(Delta) of ordered solutions on (1 - eps, infinity). This is false at Delta = 1 in the x coordinate: x = lambda u collapses to zero there and x_j ~ sqrt(2(Delta - 1)) u*_j has no finite derivative (L = 4, n = 2: u* = 1/(2 sqrt 3), which my solver reproduces, x_2/lambda = 0.288664, 0.288675, 0.2886751 at lambda = 0.3, 0.1, 0.03). The proof itself works in u through Delta = 1 and that is what is true. State the theorem as: u(Delta) is analytic and ordered on (1 - eps, infinity), x is analytic and ordered on each open side, and E(Delta) is analytic through 1. The identification is unaffected.

Numerical (check_ed_bethe.py). My Bethe solver reaches the roots by a different route from the proof, Newton continuation from the Ising-side configuration x_j = pi I_j/(L - n) down in lambda. Its energies equal the ED sector minima at flux pi per bond in all 588 sectors 1 <= n <= L/2, L = 2 to 14 (both parities), Delta in {1.0001, 1.001, 1.01, 1.1, 1.5, 2, 2.0764, 3, 5, 10, 20, 100}, maximum error 2.9e-12. Every solution is ordered and reflection symmetric.

## 3. Uniformity in the Fermi point

Section 2. I rechecked Lemma 2.1 (the overlap derivative 1_A(c + beta) - 1_A(c - beta) and |c - beta| <= min(c + beta, pi - c - beta)), r0 = (K_e/pi^2) dn from DLMF 22.11.3, J = a1 (*) r0, sum_m J^((*)m) = K, ||K'||_1 = 4/(pi sinh 2 lambda), M_1/||K'||_1 <= coth(lambda)/2, Var(ell) <= 4b, |ell'| <= 2b^2, and Lemma 2.3: (i) the translates of the out arc sit inside {|y| >= A}, (ii) c_p = 2 + 4/(e^(pi^2/lambda) - 1) from the two nearest images plus two geometric tails, (iii) both entries of kappa_A, (iv) int dn = am. All correct.

Section 3. Every constant comes from omega_Q <= int_T J = 1/2 and from the unimodality of r0 and J, and none divides by the length or mass of out. I rederived the dual equations, r0 <= r_Q <= 2 r0, m = (1 - h_Q)/2, 1/2 <= Ztilde <= 1, Var(Z) <= int |u||J'| = 1/2 - pi J(pi/2) (a translate of an arc by u loses at most d(u, 0) of its length), the inverse bound 3, |epstilde'| <= ab s_Q with vanishing boundary terms at both ends of the out arc, eps' >= 0 by the odd reduction with kernel J(x - y) - J(x + y) >= 0, and the four terms 1 + 2 + 1 + 2 in |w_v'| <= 6 a b^2 w. At Q = pi/2, out is empty, r = r0, Z = 1/2, m = 1/2, and the statements hold. Correct, uniform on (0, pi/2].

Half filling. e_bulk is used only on [0, 1/2], where it is C^2 with 0 <= e_bulk'' <= 2ab up to the endpoint (one-sided there). Both nu and m(Q_L) lie in [0, 1/2], so the curvature step does not cross the cusp at nu = 1/2. Note also that the doubling comparison uses the same function e_bulk at the same density on both rings, so the theorem does not depend on e_bulk being the true bulk energy. On even rings at nu = 1/2 either Q_L = pi/2 (then nu = m and the identity holds with h = -2a r0(pi/2)) or the first crossing is earlier, which the general case covers. On odd rings n = (L-1)/2 is an ordinary sector with out mass 2/L. The positive odd-ring defect described in PROOF Section 9 is a numerical remark with a heuristic explanation and nothing rests on it.

Near h_c1. The field enters only through the choice of the minimizing sector, and the matched comparison (L, n*) against (2L, 2n*) holds sector by sector, so no dilute-hole estimate is needed. This is right, since the sector theorem is uniform in n.

Sections 5 and 6. I rederived (I + K_Q) nu' = R, the sawtooth (Cum R runs from 1/(2L) to -1/(2L) on each gap, and on (x_n, Q] it is nu/2 - z_L in [0, 1/(2L) + negative mass]), |nu - m| <= s/(2L), the equation G + K_Q G = Cum R + c_0 - (nu - m) psi_G with c_0 = (K_Q G)(-Q) and psi_G in [0, k(2Q)] inside [0, 1], g = 13.5 s, Lemma 5.3 with ||K''||_1 = 4 M_1, the exact identity, the Peano bound Var(f')/(8L), |g_r'| <= 4(bD + D1)/(L r), and the cancellation of the half weights eps(x_j2)/(2L) between the trapezoid sum and the end. sigma_L being negative somewhere on (x_n, Q] does no harm: only the core quadrature inverts z_L, and there sigma_L >= r/2 is proved independently from Lemma 5.3.

Minor repair (a case not written). The ends paragraph assumes that a nonempty core [-A, A] contains a root, and the "empty core" paragraph covers only r0(0) < 2D/L. A nonempty [-A, A] with no root in it occurs: at L = 4, n = 2, lambda = 1.8 (s = 1/2) the roots are +-0.49587, r0(0) = 0.28242 > 2D/L = 0.27703, and A = 0.11474 (found by Astra, reproduced by me). Then int_(-A)^A r = G(-A) - G(A) <= 2g/L, each tail is at most 4dD/L with d = min(lambda/2, pi/2), and the total is s ab (2g + 8dD)/L, below the budget 2 s ab (2g + 4dD)/L. Also, when A >= Q the end integral should stop at Q, int_(x_n)^Q r = G(x_n^+) - G(Q) <= 2g/L. Neither changes K.

Numerical (check_intermediate.py, 156 rings, six values of Delta, L to 129, fillings 2, L/4, L/2 - 2, L/2 - 1, L/2). L sup|Cum R| <= 0.5000 (proved s), L|nu - m(Q_L)| <= 0.0115 (proved s/2), L sup|sigma_L - r_Q| is at most 0.0124 D, no violation.

## 4. The constant chain and M >= M_0

Section 7, recomputed in 40 digits (check_constants.json). With s = 3/2, g = 20.25, c = 2 + 4/(e^2 - 1) = 2.6260706:

| piece | value in units of ab |
| --- | --- |
| d D <= g/pi + (3/8) coth(pi) | 6.8221784 |
| ends, 2s(2g + 4 d D) | 203.36614 |
| kappa_A D1/(bD) <= (4 + 1/27) c (pi/4) coth(pi^2/2) | 8.3272947 |
| kappa_w <= 1 + ln(2c) | 2.6586358 |
| V_f | 115.62072 |
| core, V_f/8 | 14.452590 |
| curvature | 0.5625 |
| K/(ab) | 218.38123 |

So K <= 218.3813 ab and, with ab/M_0 = 4 pi^2 tanh(lambda/2)/lambda <= 2 pi^2, K <= 4310.68 M_0 <= 4311 M_0. Then (5/2)K <= 10776.7 M_0 <= 10778 M_0, and with pi^2/4 the total is 10779.2 <= 10781. The stated constants hold.

Minor repair. The displayed decimals are rounded in the wrong direction in three places: 203.4 + 14.5 + 0.5625 = 218.46 > 218.4, 4.037 * 2.07 = 8.357 > 8.35, and 218.4 * 2 pi^2 = 4311.04 > 4311. Each stated bound is true (the unrounded values are 218.381, 8.327 and 4310.68), but the chain as printed does not prove them. Print the unrounded values.

The grid evaluation of the exact lambda-dependent formula (not a proof input) reproduces: sup K/M_0 = 4153.38 as lambda -> 0 with s = 3/2, and 460.9 on my grid from lambda = 1.365 with s = 1/2 (Astra gets 461.75 at 1.36 itself), consistent with the quoted 4153.4 and 462. K/M_0 <= 218.4 * 4 pi^2/lambda does tend to zero.

Matched doubling. L[F_(2L)(2n) - 2F_L(n)] <= L[K/(2L) + 2K/L] = 5K/2, valid because both rings have the same lambda and density. Sectors n > L/2 follow by spin flip.

Reductions (C129, which I checked by hand rather than taking as cited). (i) The sector energy at any flux is at least the one in the stoquastic gauge (compare <psi|H|psi> with <|psi||H_stoq||psi|>), so E_L >= min_n [F_L(n) - B(L/2 - n)]. (ii) On the doubled ring the flux differs from the minimizing one by at most pi modulo 2 pi, so the twist per bond is at most pi/(2L) and the trial-state cost is at most (2L/2)(1 - cos(pi/(2L))) <= pi^2/(8L) <= pi^2 M/(4L), using M >= 1/2. (iii) The field term is linear in n and cancels in matched sectors. (iv) M >= M_0: a gauge term A x 1 - 1 x A has zero expectation in any two-site state with equal marginals. The states |up up>, |down down> and (|up down> +- e^(i phi)|down up>)/sqrt 2 have equal marginals and energies J_z/4 -+ B/2 and -J_z/4 +- J/2, so no gauge has a spread below that of h itself, and for Delta >= 1 the half-spread is at least (J_z/2 + J/2)/2 = M_0. Correct on this side.

Glue at Delta = 1. At fixed L, field and DM angle, E_L, E_(2L) and M are continuous in Delta and M(Delta = 1) >= 1/2, so either bound passes to Delta = 1. Correct.

## 5. Independent numerics

All with my own code (ref148.py). No solver failure in the final runs. (My first continuation, with a fixed step in lambda and a fixed Newton tolerance, lost 15 half-filled rings at Delta <= 1.2 and L >= 128. An adaptive step with a tangent predictor and a tolerance scaled with L removed the problem. It was in my solver, not in the equations.)

- Bethe against ED: 588 sectors, maximum error 2.9e-12 (Section 2 above).
- ED sector defect L(F_(2L)(2n) - 2F_L(n))/M_0, L = 2 to 10, all 1 <= n <= L/2, nine values of Delta (225 cases): maximum 3.99993 at Delta = 1.0001, (L, n) = (2, 1), decreasing with Delta (3.9993, 3.9933, 3.934, 3.678, 3.381, 2.877, 2.170, 0.724 at Delta = 1.001, 1.01, 1.1, 1.5, 2, 3, 5, 20).
- End-to-end ED of L(E_(2L) - 2E_L)/M with the full bond, L = 2 to 6, DM angle in {0, 0.3, 1, pi/2, 2.5, pi}, 24 fields in [-0.37, 1.1] h_c2, nine values of Delta (6480 samples): maximum 3.99993 (Delta = 1.0001, L = 2, B = 0), minimum -17.57 (Delta = 20, L = 5, antiferromagnetic transverse sign on an odd ring). The maximum at every Delta equals the sector maximum, as the reductions predict.
- Bethe doubling to 2L = 2048 and sector error against my own Nystrom continuum (scan_bethe.jsonl): 1885 pairs (L, 2L) with no solver failure, 13 values of Delta from 1.0001 to 100 (including Delta = cosh 1.36 = 2.0764), L in {2 to 9, 12, 13, 16, 17, 32, 33, 64, 65, 128, 129, 256, 257, 512, 513, 1024}, fillings n in {1, 2, L/8, L/4, L/3, L/2 - 3, L/2 - 2, L/2 - 1, floor(L/2)}. The defect L(F_(2L)(2n) - 2F_L(n))/M_0 has maximum 3.99993 at (L, n) = (2, 1), Delta = 1.0001, in agreement with ED. For L >= 128 its maximum is 2.4714 (Delta = 1.0001, L = 128, half filling, against C148's 2.4711 at Delta = 1.001), 1.55 at Delta = 2, 0.224 at Delta = 20 and 0.047 at Delta = 100. It is never below -4e-10, which is rounding. The sector error L|F_L(n) - L e_bulk(n/L)|/M_0, on 3504 rings (all fillings for lambda >= 0.044, half filling only at Delta = 1.0001, where the Nystrom grid would be too large), has maximum 2.455 at L = 2 and is at most 1.65 for L >= 128, against the proved 4311. It is stable under doubling (for example -1.5230 and -1.5226 at Delta = 1.2, n = 253 on L = 513 and its double), which is the K/L law. Near h_c1 at Delta = 2 the values at L = 1024 and n = L/2 - 1, L/2 - 2, L/2 - 3 are 0.091, 0.177, 0.258, the same on the doubled ring to three digits, and zero to rounding at n = L/2. At Delta = 5 they are 0.014, 0.027, 0.041. So the dilute regime shows no growth with L.
- Half-filling parity (my solver): even L, the defect over M_0 is 1.84, 0.399, 7.4e-5, 0 to rounding at Delta = 2 and L = 8, 32, 130, 512. Odd L, n = (L-1)/2: 1.821, 1.384, 0.508, 0.137 at Delta = 2 and L = 9, 33, 129, 513, and 0.796, 0.292, 0.081, 0.021 at Delta = 5. These reproduce C148's numbers and extend them by one size. h_c1 = 2a r0(pi/2) = 0.389802 at Delta = 2 against F_L(L/2 - 1) - F_L(L/2) = 0.3900 at L = 1024, and 3.120998 against 3.121058 at Delta = 5.

One correction to a quoted number. RESULT and REPORT give the largest sector error as L|F - L e_bulk| = 1.82 M_0. That is the largest over C148's sample, which starts at L = 4. At L = 2, n = 1 it is 2.45 M_0 as Delta -> 1+ (F_2(1) = -1 - Delta/2, e_bulk(1/2) -> 1/4 - ln 2). This is immaterial against 4311 M_0.

## 6. Repairs, most significant first

1. Section 6: add the case of a nonempty core with no root in it (realized at L = 4, n = 2, lambda = 1.8), and stop the end integral at Q when A >= Q. Same constant.
2. Theorem 4.5: state analyticity through Delta = 1 for u(Delta) and E(Delta), not for x(Delta).
3. Section 7: print unrounded intermediate values (218.381 ab, 8.327, 4310.68 M_0).
4. Lemma 4.1(a): write the subcase in which x and its antipode lie in the same gap.
5. Decimals: 0.02099/pi in place of 0.021/pi, 0.805 in place of 0.804.
6. Lemma 2.2(i): remove the second use of the symbol q.
7. RESULT/REPORT: the largest computed sector error is 2.45 M_0 (L = 2), not 1.82 M_0.

## 7. What I checked and what I only read

Checked by hand, line by line: PROOF Sections 1 to 8, including every lemma that is cited "verbatim" from C139/C140/C142 wherever the circle changes something (the pairing, the two-sided out arc, the tails, the restricted inverse, the sawtooth, G, the ends). The C129 reductions used on this side (stoquastic gauge, twist cost, field cancellation, M >= M_0) by my own short arguments.

Checked by Astra (gpt-6.2-astra) and then by me: the same sections in one adversarial pass with C142 Part I attached. Its findings were repairs 1 to 6 above. I had found 3 and the content of 1 independently, and I verified 1, 2, 4 and 5 numerically or by hand before listing them. It reported no error that changes a theorem.

Checked numerically with my own code: Section 5 above, the interval certificate, every Section 7 number, and the Sections 4 to 6 inequalities on Bethe roots.

Only read, not rederived: C142 Theorem I.4 (the branch and identification on 0 <= Delta < 1, refereed in C147) and C142 Lemma I.2 (cosh(pi u_n) <= L), which Theorem 4.5 uses as its starting point. The grid suprema 4153.4 and 462 are numerical and are not proof inputs.

## 8. Overall verdict

The statement as written, "for every bond of the U(1) class with J_z/J > 1, every field, DM angle and L >= 2, L(E_(2L) - 2E_L) <= 10781 M", is correct after minor repairs. The sector theorem |F_L(n) - L e_bulk(n/L)| <= K/L holds with K <= 4310.68 M_0 <= 4311 M_0, uniformly in lambda > 0, L >= 2 and 0 <= n <= L/2, and the assembly gives 10779.2 <= 10781. The positivity step is rigorous: the interval certificate is sound, its monotonicity directions are the right ones, and it also passes without assuming them. The identification of the ground state needs no completeness. The repairs are an unwritten case in the end estimate, the coordinate in which Theorem 4.5 is stated, and several decimals rounded the wrong way. None touches a constant. The result rests on the earlier refereed C142 Theorem I.4 for the branch below Delta = 1. The status in REPORT Section 9.106 can change from "proved, not refereed" to "refereed, correct after minor repairs".

## Files

ref148.py (ED, circle Bethe solver, continuum), check_ed_bethe.py (check_ed_bethe.json), check_full_ed.py (check_full_ed.json, log_full_ed.txt), scan_bethe.py (scan_bethe.jsonl, log_scan_bethe.txt), check_intermediate.py (check_intermediate.json), check_constants.py (check_constants.json). All ran locally, each in under two minutes. astra62_read_spec.md and astra62_read.out are the spec and the raw output of the gpt-6.2-astra read.
