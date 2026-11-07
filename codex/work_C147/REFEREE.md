This referee report was written by an independent Claude Opus agent (2026-11-07), not by Codex.

# C147 referee report on C142 (Parts I and II) and C145

J=1 throughout. Read: specs/C147_referee_C142_C145.md, specs/_Cd_target.md, specs/_bound_background.md, work_C142 (RESULT, PROGRESS, PROOF, constants_rep.py, check_repulsive.py, check_attractive.py, ed_check.py and their outputs), work_C145 (RESULT, PROGRESS, PROOF, constants_C145.py and outputs), and for the citations work_C140 (PROOF, REFEREE), work_C139 (PROOF), work_C138 (PROOF), work_C136 (PROOF Sections 1 and 2). I did not reuse any solver from those folders. Every script here was written from scratch (ref_core.py holds the ED and both Bethe solvers).

## Verdicts in one place

| item | verdict |
| --- | --- |
| C142 Lemma I.1 and its insertion into C140 Lemma 3.3 (Lemma I.2, Theorem I.3) | correct |
| C142 Theorem I.4 (continuation on all of 0 < gamma <= pi/2, theta_1 >= 0.168, negative mass) | correct after minor repairs (wording of the gamma_* step, n = 1 clause) |
| C142 Lemmas I.5 to I.8, Theorem I.9 | correct after minor repairs (Q = infinity case of I.5, one misstated location of a supremum) |
| C142 Corollary I.10, 5K + pi^2/4, seam 4L <= 112, Delta = 1 by continuity | correct |
| C142 Part II, Lemmas II.1 to II.5 and Theorem II.6 | correct (one constant rounded down in the third digit) |
| C145 Theorem A.2 and the reduction to (G), (W1'), (W3') | correct as a conditional theorem (one constant rounded down in the third digit) |
| C145 Theorem B.1 (W2 by reflection) | correct |
| C145 Lemma B.2, Proposition B.3 | correct |
| C145 claim that C138 Theorem B is correct, arithmetic 124.87, 627.5, 647.4, 686 | correct, every number reproduced |

I found no error that changes any stated theorem or constant. The problems are presentational or in the third significant digit of a majorant, and each has a one-line repair given below.

## 1. C142 Part I, the repulsive side

### 1.1 Lemma I.1 and its use

Checked by hand. With t = exp(-bQ), int_Q^inf sech(bx) dx = (2/b) atan(t), so int_Q^inf r0 = atan(t)/(gamma b) and r0(Q) = t/(gamma(1+t^2)). The ratio is f(t)/b with f(t) = atan(t)(1+t^2)/t, and f'(t) = [t - (1-t^2) atan t]/t^2 >= 0 on (0,1], f(1) = pi/2, f(0+) = 1. Correct. Numerically (check_constants.py) f is increasing on a 1e5-point grid, and the ratio of the two sides is 1 at Q = 0 and 2/pi in the far tail, as stated.

Insertion into C140 Lemma 3.3 (Lemma I.2). The chain 1/(2L) <= int_(x_n)^inf sigma <= int_Q^inf s_Q = h_Q/2 <= int_Q^inf r0/(1-q_J) (Q = x_n) is C140's refereed chain, and only its last inequality changes. With 1 - q_J = pi/(2(pi-gamma)) this gives L r0(x_n) >= b/(2(pi-gamma)) >= b/(2pi). The ratio to lambda_* b = (2-ln2) b/pi^2 is pi^2/(2(pi-gamma)(2-ln2)) >= 1.20197 (recomputed on a 4e5-point grid). Correct.

Theorem I.3. From C140 Theorem 3.4, sigma >= r0 - (2-ln2)/(pi gamma L) on [x_1, x_n], and 1/L <= 2(pi-gamma) r0(x_n)/b with pi gamma b = pi^2. So theta_1 = 1 - 2(2-ln2)(pi-gamma)/pi^2, with minimum 0.168032 as gamma -> 0. Correct.

### 1.2 Theorem I.4, the continuation argument

The structure is sound. (P) is closed along the branch by joint continuity, and open because Theorem I.3 upgrades (P) to the strict margin theta_1 r0 > 0 on a compact set. The limit point at gamma_* has a positive definite Hessian because sigma(x_j(gamma_k)) >= theta_1 b/(2(pi-gamma_k)L) passes to the limit, and C139 Lemma 2.2 then applies. Step 3 is the open-closed argument that C140's referee asked for, and it is correct. The identification by analyticity of the simple stoquastic ground eigenvalue and of the Bethe energy, equal near gamma = pi/2, needs no nonvanishing of the Bethe vector away from the free point. Correct. Part (ii): phit is even because psi and J_out commute with reflection, and ||psi||_1 = q_J gives int_(y > x_n) phit <= q_J/(2(1-q_J)) = (pi-2gamma)/(2pi). Correct.

Minor repairs.
(a) gamma_* is defined as an infimum over gamma_1 in [gamma_0, pi/2], so gamma_* >= gamma_0 always. The sentence "This contradicts the definition of gamma_* unless gamma_* < gamma_0" should read: Step 3 shows the branch extends below gamma_*, so gamma_* = gamma_0 and the branch covers [gamma_0, pi/2] (Step 3 applied at gamma_0 itself gives the closed endpoint).
(b) For n = 1 the proof cites C140 Theorem 4.1, which gives sigma > 0 on R but not the clause (i) sigma(0) >= theta_1 r0(0). Clause (i) is not used downstream at n = 1 (Lemma I.5 needs only sigma > 0, Lemma I.8 does not use (i)), so either state (i) for n >= 2 only, or note that the proof of Theorem I.3 goes through at the single point x_1 = 0.

Numerical (mine, check_ed.py, check_rep.py). Continuation by Newton from the free solution at gamma = pi/2 to the target gamma, for every 1 <= n <= L/2, L = 2 to 12, Delta in {0.3, 0.7, 0.9, 0.99, 0.999}. The ordering is never lost, and the minimum Hessian eigenvalue divided by L along every path is at least 0.519. Bethe energies equal ED sector energies at flux pi per bond to 5.2e-13 (largest at Delta = 0.999, from the stiff continuation; 1e-14 elsewhere). Twisted ED (check_flux.py, L = 3 to 10, all sectors, Delta in {0.5, 0.99, 0.999, 1, -0.5, -0.99}) confirms that flux pi per bond is the sector minimizer.

### 1.3 Lemmas I.5 to I.8 and Theorem I.9, term by term

Lemma I.5. On the gaps, Cum R goes from 1/(2L) down to -1/(2L) because z_L increases. On (x_n, Q_L], Cum R = nu/2 - z_L lies in [0, (1/2 + (pi-2gamma)/(2pi))/L] = [0, s/L]. Correct. Minor repair: at even half filling with Q = infinity, the sentence "the same argument applies" hides a sign change. There Cum R(x) = int_x^inf sigma, which can be negative, and the bound is |Cum R| <= max(s, (pi-2gamma)/(2pi))/L = s/L, using the negative mass bound once more. The conclusion is unchanged, but the sentence should be written out.

Lemma I.6. Mismatch: nu - m = int Z dR = -int Z' Cum R, so |nu - m| <= q_J s/L (C139 Lemma 3.4). Cumulative error: I rederived the integration by parts, G + K_Q G = S + c_0 - (nu-m)psi with psi(x) = k(x-Q) + k(2Q) in [0, k(2Q)], and then G = g_1 + c_0 Z with |c_0| Z <= (1+q)||g_1|| from 1/(1+q) <= Z <= 1. So ||G|| <= (2+q)[(1+q_J)/(1-q_J)] s (1 + q q_J/2)/L. Each factor decreases in gamma (checked on a grid), so g <= 45/4 at gamma -> 0 (grid sup 11.249999). Correct.

Lemma I.7. sigma - r = -K*nu' = -int K'(x-u) G(u) du on R (boundary terms vanish since K -> 0). ||K'||_1 = 2K(0) = 2 cot(gamma)/pi <= 2/(pi gamma). ||K''||_1 = 4 max|K'| and the bound M_1 = (2/pi)(0.325/sin^2 gamma + 0.25/sin gamma) (I rederived 9/(16 sqrt3) = 0.32476 and 1/4). Correct.

Lemma I.8. Core: A exists iff L >= 4 gamma D = 8g/pi <= 28.65, so L >= 29 suffices. The trapezoid Peano bound Var(f')/(8L), |g'| <= 4(bD + D1)/(L r), and int_(-A)^A dx/r0 = (4gamma/b) sinh(bA) <= L/(bD) are correct. 4ab D1/(bD) = 8 pi^2 gamma sin(gamma) M_1 and gamma sin(gamma) M_1 <= 0.575 (grid max 0.575 at gamma = pi/2). So V_f <= 8pi^2 + 80pi + 45.40 = 375.69 and V_f/8 <= 46.96. Ends: I rederived the end identity (sum - L int eps sigma + eps(x_j2)/2 = -L int eps' Cum R), int_(x_j2)^A r = G(x_j2^+) - G(A) <= 2g/L when no root lies in (x_j2, A], and int_A^Q r <= c_r (pi/(2b)) r0(A) = c_r pi D/(bL) with D/b = 2g/pi^2. Both ends together are at most 2 s ab g (2c_r/pi + 2)/L. The empty-core bound s ab [2g + 4c_r g/pi]/L is below this. Curvature: a2 (q_J s)^2/2 <= 2pi^2/8 = pi^2/4. All correct.

Theorem I.9, majorants recomputed (check_constants.py). Core 46.96, ends 2 pi^2 (45/4)(2(1+4/pi)/pi + 2) = 765.50, curvature 2.467, sum 814.93 <= 815. The exact gamma-dependent K_large has grid supremum 811.30 at gamma -> 0. For L < 29 the bound |F - L e_bulk| <= L (bond spectral diameter 1) gives L|F - L e_bulk| <= 784 <= 841. Correct. One misstatement: the text says the supremum is approached as gamma -> 0 "where V_f/8 -> 43.42". The gamma -> 0 limit of V_f/8 is 43.33, and 43.42 is its maximum, attained near gamma = 0.12. This does not affect any bound.

Numerical (check_rep.py, mine, 128 cases: Delta in {0.5, 0.9, 0.99, 0.999}, L in {6, 13, 32, 64, 129, 256}, six fillings each, restricted continuum by Gauss-Legendre Nystrom). Every inequality holds. L r0(x_n)/b >= 0.465 (proved >= 0.159), min sigma/r0 on [x_1, x_n] >= 0.979 (proved >= 0.168), negative mass beyond x_n zero, L sup|Cum R| <= 0.5000 (proved s), L|nu - m| <= 0.012, L||G|| <= 0.521 (proved g <= 11.25), L sup|sigma - r_Q| <= 0.81, L|F - L e_bulk| <= 0.86, and L|F - L e_bulk|/K_large <= 0.0051. These agree with C142's own table.

### 1.4 Corollary I.10 and the Heisenberg end

Matched doubling: F_(2L)(2n) - 2F_L(n) <= K/(2L) + 2K/L, so L[...] <= 5K/2 = 5K M_0 for L >= 29 (both L and 2L use K_large <= 815). For L <= 28, seam gluing costs at most twice the bond diameter, 2 = 4M_0 per doubling, so L[...] <= 4L M_0 <= 112 M_0. Adding the flux cost pi^2 M/4 and M >= M_0 gives 5 * 815 + pi^2/4 = 4077.47 <= 4078. The grid version gives 4058.96. At Delta = 1, E_L, E_(2L) and M are continuous functions of a finite matrix at fixed L, and M(Delta = 1) = J/2 > 0, so the bound passes to the limit at every fixed field and DM angle. Correct. The field, DM and flux reductions are cited from the refereed C129/C136 assembly and I did not re-referee them.

ED (check_ed.py, mine). Sector defect L(F_(2L)(2n) - 2F_L(n))/M_0 for L <= 6: maximum 3.864, 3.987, 3.9987, 4.000 at Delta = 0.9, 0.99, 0.999, 1, all at (L, n) = (2, 1), and 2.3431 = 8 - 4 sqrt2 at Delta = 0 (the exact free value of C129).

## 2. C142 Part II, Theorem II.6

Lemma II.1. The midpoint cells tile [-nu/2, nu/2] = z_L([-Q_L, Q_L]) (I_1/L - 1/(2L) = -nu/2), each X_j carries sigma-mass 1/L and contains x_j, so int f dR_L = sum_j int_(X_j) [f(x_j) - f(x)] sigma dx. Correct.

Lemma II.2. (I - T_Q)(sigma - r) = T_Q R_L on [-Q, Q], using that the roots lie in (-Q_L, Q_L). Correct.

Lemma II.3. calR = sum_(k>=1) T^(*k) is even and nonincreasing on x > 0. Wintner's theorem (convolutions of symmetric unimodal densities are symmetric unimodal) applies to every term, and the monotone limit preserves the property. I rederived That/(1 - That) = sinh(A)/(2cosh(A + eta w/2) sinh(eta w/2)), the two pointwise bounds 1/(e^(eta w) - 1) and q/(1-q), the split at w = 1, and log(1/(1-e^(-eta))) <= log(1/eta) + eta/2 from 2 sinh(eta/2) >= eta. Correct. Numerically (check_att.py, my own overflow-safe quadrature) calR(0) = 0.2387, 1.3464, 7.1496, 31.05 at 1 + Delta = 0.5, 0.1, 1e-2, 1e-3, against the bound 0.463, 1.671, 7.930, 33.29, and calR(0) >= T(0) in every case.

Lemma II.4. The outside-channel comparison w = U - u = (I - T_Q)^(-1) T*(1_out U_y) >= 0 and the identity int_out K_out(x, u) du = 1 - (1-q) Z_Q(x) <= q are correct. The variation of a nonnegative superposition of translates of calR is at most the superposition of the variations. Correct, uniform in Q including Q = infinity. Numerically (Nystrom, x in {0, Q/2, Q - eta, Q}, Q in {0.3, 1, 3}) Var_y calR_Q(x, .) <= 0.47, 2.63, 14.1, 60.7 against 2(1+q) calR(0) = 0.64, 4.61, 27.3, 122.4.

Lemma II.5. calR'' = T'' + T'*calR', ||T''||_1 = 4 max|T'| (T has the form of C139's K with gamma replaced by eta, so the unimodality of T' carries over), and int_out |d_x K_out(x,u)| du <= 2T(0)(1+q). Correct. Numerically the variation of d_x calR_Q(x, .) is at most 0.70, 5.19, 55.0, 553 against the bound 1.33, 14.8, 231.5, 3053.

Theorem II.6. The arithmetic checks, with eta T(0) <= 1/pi and 2 + q <= 3. One rounding: the constant written 3.61 is sup (8/pi)(0.32476 (eta/sin eta)^2 + 0.25 eta^2/sin eta) = 3.6113 at eta = pi/2. Write 3.62.

Numerical (check_att.py, my own step-free attractive Bethe solver, checked against ED to 1e-14 at L <= 12 for 1 + Delta in {0.5, 0.1, 1e-2, 1e-4}). 72 cases, 1 + Delta in {0.5, 0.1, 1e-2, 1e-3}, L in {16, 64, 256}, six fillings including half filling. Here sigma_L' - r_Q' is evaluated up to and including the Fermi point, using the exact derivative of the Nystrom interpolant, which avoids the mesh issue that made C138 and C142 stop eta/4 short of Q_L. Results: max eta L sup|sigma - r_Q| = 0.037, 0.123, 0.233, 0.288 and max eta^2 L sup|sigma' - r_Q'| = 0.085, 0.183, 0.216, 0.214, against C_D(eta) = 0.67, 2.08, 3.87, 5.48 and C_D1(eta) = 1.46, 3.02, 4.64, 6.11. The largest ratio to the proved bound is 0.061. So (D1) holds numerically at the Fermi point as well, which settles C145 Remark 1 at the sampled values.

Part II.4, the conditional bound. Inserting Theorem II.6 into C138 Theorem B with M = 4C_D(eta) is legitimate (Theorem II.6 holds on all of [-Q, Q]). The tail coefficient 16 pi (1 + 2/0.7) = 193.88 is correct.

## 3. C145

### 3.1 Theorem A.2

I rederived the proof line by line against C138 Theorem B.
- Comparability on K_1 = {lambda >= 4C_D(eta)} follows from Theorem II.6 at every point, and on S from (G). So g <= 4/3 on K_0.
- eta Var_(K_1)(g) <= (16 eta/9)[(C_D/(eta L))(2/rho_1) + (C_D1/(eta^2 L)) 2A_core/rho_1] = (8/9)[eta + (C_D1/C_D) A_core]. Correct.
- C_D1/C_D = 2 eta M_1^T/((1+q) calR(0)) + 2(2+q) eta T(0)/(1+q), with calR(0) >= T(0). The first term is at most 4(0.32476 eta/sin eta + 0.25 eta), whose supremum is 3.6113 (the text writes 4 * 0.9026 = 3.61, and the bracket is 0.9028), and the second is at most 4/pi. So C_D1/C_D <= 4.8846 and C_G1 <= 13.1958, not 13.19. Repair: write C_G1 <= 13.20. Effect: K_G = 136.610 and C_strip = 685.52 <= 686, so the stated 686 stands.
- A_core = 2.71765 at d = (1/2) log(4/(2 - sqrt3)) = 1.3516, and A_tail = 3.98873 at d = (1/2) log(4/(3 - 2sqrt2)) = 1.5745 (closed forms; C145 used a grid). Both correct.
- The tail {lambda < M_0}, the straddling cells (8a/(3L) each, using sigma >= 3r/4 on K_0), the mismatch with (W3), the half-filled case (Var(v) = 2ab <= 8 eta), and the degenerate case M_1 <= M_0 (S empty, K_0 inside K_1) are all correct.
- Final terms recomputed: 11.149, 3.290, 100.248, 16.755, 5.168, sum 136.61, and 5 * 136.61 + pi^2/4 = 685.52 <= 686.

So the reduction of the strip to (G), (W1'), (W3') is correct as stated. (G) is genuinely weaker than (D), (D1) on S. Given (D), (D1) with log-free constants and C_D <= 1/2, one gets eta Var_S(g) <= 32(eta C_D + C_D1 A_core)/(9 M_0).

### 3.2 Theorem B.1, (W2) by reflection

Correct, and a nice argument. psi_+(x) = r0(2Q - x) solves the full-line equation with source a1(2Q - .) >= 0, and for u > Q >= 0 one has |2Q - u| <= |u|, so (I - T_Q) psi_+ >= g_+. Positivity of (I - T_Q)^(-1) gives h <= psi_+ + psi_-, and the ratio bound 2 exp(-2b(Q - |x|)) per term follows from cosh(bs) >= e^(bs)/2. Numerically (Nystrom, own code, 1 + Delta in {0.5, 0.1, 1e-2, 1e-3}, Q in {0.3, 1, 2 or 3}) h >= 0, h <= psi_+ + psi_- with positive slack in every case, and r_Q/r0 >= 1 - 4e^(-2b d) at every node.

### 3.3 Lemma B.2 and Proposition B.3

Lemma B.2. |a1'| = 2 sinh(2|x|) a1/(cosh 2x + cos eta) <= 2 tanh(2|x|) a1 for cos eta >= 0, and (I - T_Q) eps_Q' = -a a1' (the boundary terms vanish since eps_Q(+-Q) = 0). Correct. Numerically max |eps_Q'|/(2a tanh(2Q) r_Q) <= 0.75.

Proposition B.3. (i) and (ii) follow from Lemma B.2 and C136 (5.3). I also checked the structural identities (I - T_Q) r' = a1' - r(Q) J_0, (I - T_Q) Z' = -Z(Q) J_0 and v_Q(Q) = a(Z'/Z - r'/r)(Q) (from eps_Q = a(r(Q) Z/Z(Q) - r)). Correct. Spot values of the two assumptions (check_att.py): v_Q is monotone on [0, Q] in all 12 sampled (eps, Q), eta Z_Q(Q)^2/tanh(2Q) is at least 0.965 (1 + Delta = 1e-3, Q = 0.3) and tends to about pi/2 at large Q, and e_bulk''/eta^2 <= 1.24. This is consistent with C145's minimum 0.808 at smaller eta and Q/eta of 20 to 70. The large-Q value pi/2 matches the half-line Wiener-Hopf edge value Z(edge)^2 = 1/(1 - That(0)) = B. So (W3') is plausible, but it remains an assumption.

### 3.4 C138 Theorem B and its arithmetic

I checked Lemma 1 (midpoint Peano kernel, Var_cell(f')/(8L)), the core (product rule, sigma^2 >= 9r^2/16, int |r'|/r^2 <= 2/rho, the split of int dx/r), the tail (int_y^inf r0 <= e^(-by)/(gamma b) <= 2 r0(y)/b, including the case rho/kappa > r0(0), where 1/4 <= 1/pi), the straddling cells and the mismatch. All correct. Arithmetic (check_constants.py): the terms of (4.3) are 7.480, 2.094, 96.941, 16.755, 1.604, sum 124.874 <= 125, and 5 * 125 + pi^2/4 = 627.47 <= 628. With the proved kappa(d) and optimized d, the sum is 128.988 and C = 647.41. With (G), (W1'), (W3') it is 686. All reproduced.

### 3.5 Tagging

Part D.1 of C145 labels as [Proved] two obstruction computations ("the computation behind the obstruction" for route (a), and the continuity argument). These rest on the asymptotic form calR'(u) about -1/(pi eta u) and on informal charge counting. They are correct as heuristics about why those particular methods fail, but they are not theorems. Retag as [Heuristic]. The "shell holds of order log(1/eta)/eta roots" count in A.1 is near-half-filling bookkeeping and also not load-bearing.

## 4. What I checked numerically and what I only read

Checked numerically with my own code: Bethe versus ED on both sides at L <= 12; Hessian positivity along the repulsive continuation; flux minimality by twisted ED at L <= 10; the ED doubling defect at L <= 6; every Part I inequality on Bethe roots up to L = 256; every Part I, Part II, C145 and C138 constant (closed forms and grids); calR(0) and its bound; Lemmas II.4 and II.5 by Nystrom; Theorem II.6 on Bethe roots including the Fermi point; Theorem B.1, Lemma B.2, and spot values of (W1'), (W3').

Only read, not recomputed: the C129/C136 field, DM, flux and M >= M_0 reductions (refereed earlier); C139 Lemmas 3.1 to 3.5 and C140 Lemmas 3.1 to 3.3 (refereed earlier; I rechecked only the parts that Part I changes or relies on in a new way); C145's shell numerics (check_shell*.json), which I did not rerun at L up to 16384.

## 5. Overall verdicts

The claim "L(E_(2L) - 2E_L) <= 4078 M on 0 <= J_z/J <= 1, every field, DM angle and L >= 2". Correct after minor repairs. Every new step of C142 Part I holds as written, and the constant 4078 follows from the analytic majorants (4077.47, or 4058.96 with the exact gamma-dependence). The repairs are the gamma_* wording in Theorem I.4, the n = 1 clause (i), the Q = infinity sentence in Lemma I.5, and the location of the 43.42 maximum. None of them touches a constant. The proof relies on the earlier-refereed C139/C140 lemmas and the C129/C136 reductions.

The conditional strip theorem (C145: under (G), (W1'), (W3'), |F_L(n) - L e_bulk| <= 136.61/L and L(E_(2L) - 2E_L) <= 686 M on -1 < J_z/J <= 0, hence 4078 M on the whole U(1) class). Correct as a conditional statement. The only correction is C_G1 <= 13.20 in place of 13.19, which leaves 685.52 <= 686. (W2) is now proved (Theorem B.1). (G), (W1') and (W3') remain unproved hypotheses, supported by numerics with wide margin, so the strip -1 < Delta < -1 + delta is still not covered unconditionally. Also confirmed: C142 Theorem II.6, the log-weak (D) and (D1), is correct and holds pointwise up to the Fermi point.

## Files

ref_core.py (ED and both Bethe solvers), check_ed.py (check_ed.json), check_rep.py (check_rep.json), check_att.py (check_att.json, log_att.txt), check_flux.py (check_flux.json), check_constants.py (check_constants.json). Each runs locally in under a minute.
