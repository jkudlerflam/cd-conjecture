# Referee report on C136 (all fillings on Delta in [-1+delta, 0])

2026-11-07. Independent referee pass on PROOF.md Sections 1 to 7 and 9, with RESULT.md, CHECKS.md, C134's PROOF.md and REFEREE.md, and C129's PROOF.md Sections 1 to 6 read for the baseline. Section 8 (the crossover obstruction) makes no positive claim and was only skimmed. Check script `referee_checks.py`, output `referee_checks.json` (the moment check was rerun separately with `--moment-only` after fixing an overflow in its first version; the numbers below are from the rerun).

## Summary

I re-derived every identity and inequality that enters K_all(delta) and the assembled constant. I did not find a mathematical error, a hidden log L, a constant that depends on the Fermi boundary Q, or a step that silently assumes half filling. The issues below are expository. None of them changes a hypothesis, a constant, or a conclusion.

## What I checked

### Sector identification at general n (PROOF Section 1)

The Yang-Yang recession slope on the top-k ray is k(pi-gamma)(L-2n+2k), which is positive for every n <= L/2, not only at half filling. The Hessian argument does not depend on n. At gamma = pi/2 the quantum numbers I_j = j-(n+1)/2 are integers for odd n and half-integers for even n. That is the Jordan-Wigner boundary condition exp(ikL) = (-1)^(n-1) of C129 Section 7 for K_L, and the occupied block is the symmetric Fermi sea, so the Bethe vector is the free sector ground state. The largest pair sum of momenta is 2pi(n-2)/L < pi, so the contact denominators are nonzero. Perron-Frobenius in sector n (C129 Section 6) holds for every 0 < n < L. The flux minimizer Lpi mod 2pi is proved in C129 for every sector by the modulus argument. Nothing here uses n = L/2.

### Finite-interval resolvent and density (Section 2)

T_Q has norm at most q on both L1 and L-infinity for every Q, so the Neumann bound B is uniform in Q, including Q -> 0 and Q -> infinity. Termwise comparison gives a1 <= r_Q <= r_infinity and 1 <= Z_Q <= B. The level-set proof that r_Q and Z_Q are even and nonincreasing on [0,Q] is correct. I re-derived (2.3): the boundary integrals are (T_Q Z_Q)(+-Q) = Z_Q(Q)-1, giving m' = 2 r_Q(Q) Z_Q(Q). For (2.4), the inside deficit integrates to at most Bq times the outside mass, so 1/2-m <= (1+Bq) times the outside mass, and 1+Bq = B. The tail integral is at most (2/pi)exp(-bQ). Correct, with the stated direction.

### Cumulative-root bound uniform in density (Section 3)

2c-nu = (alpha/pi)(1-2nu) is correct. C_nu = P_nu(z_nu) holds because s_Q = r_Q inside and z_nu(Q) = m/2 = nu/2. The identity z_L-z_nu = T*(C_mu-C_nu) has no constant of integration by oddness. The fixed-point argument ||C_mu-C_nu|| <= 1/(2L)+q||C_mu-C_nu|| gives B/(2L) with no inverse endpoint density, so the bound is uniform down to nu -> 0 and up to nu = 1/2. This is the key new step and it is right.

### Finite-boundary quadrature identity (Section 4)

(I-T_Q)(mu-r dx) = R_L on [-Q_L,Q_L] holds because mu is supported strictly inside. I re-derived (4.3) from e-h = (I-T_Q)epsilon_Q and self-adjointness: F_L = L[E(Q)+h(nu-m)] + L integral epsilon_Q R_L exactly. The chemical-potential subtraction h = -a r_Q(Q)/Z_Q(Q) makes epsilon_Q(+-Q) = 0, which is what removes the first-order endpoint error. (4.2) follows from Stieltjes integration by parts with the sawtooth vanishing at both endpoints, and Var(Z) <= 2(B-1). Numerically the identity (4.3) closes to 7e-14 on every case I ran.

On the question of a square-root singularity at the Fermi point: none enters. In the rapidity variable epsilon_Q vanishes linearly at Q, with epsilon_Q'(Q) > 0, and the quadrature is done in the counting coordinate t = z_L(x), where f' = epsilon_Q'/sigma_L. The square-root behaviour of the energy as a function of momentum near the band bottom never appears. In the dilute limit epsilon_Q'(Q) -> 0 linearly in Q, which only helps.

### Bulk curvature (Section 5)

dE/dQ = 2 r_Q(Q)[(I-T_Q)^(-1)e](Q) and e_bulk' = h are correct. (5.1) uses epsilon_Q(+-Q) = 0 to kill the boundary terms. Positivity of epsilon_Q' on (0,Q] follows from the odd reduction with kernel T(x-y)-T(x+y) >= 0. The bound |a1'| <= 2a1 holds for Delta <= 0. partial_Q epsilon_Q = -h' Z_Q and the boundary-value differentiation give e_bulk'' = epsilon_Q'(Q)/[2 r_Q(Q) Z_Q(Q)^2], which lies in [0,a]. The checks h(0) = -(1+Delta) and e_bulk = -sin(pi nu)/pi at Delta = 0 are correct. So the density mismatch |nu-m| <= B/L costs at most aB^2/(2L), with a sign. Numerically e_bulk''/a never exceeded 0.49999 (Delta = -0.5) or 0.168 (Delta = -0.9).

The restricted constant (5.7) also checks: m(Q_L) <= 1/2-eta/2 for L >= 2B0/eta, hence Q_L <= R0, and the sup bounds S1, U1, U2, the bound on ||T'||_infinity by 2/(pi s0^3), and |a1''| <= 8a1 are all correct.

### Exponential resolvent comparison (Section 6.1 to 6.2)

The moment formula (6.2) is correct. I checked it by direct quadrature at Delta = -0.999, -0.99, -0.9, -0.5, -0.1, with agreement to 2e-15. The L-infinity/L1 convolution bound with weight exp(beta x) gives R_infinity <= C0 exp(-beta|x|). The decay is real but not sharp. The pole of 1/(1+Khat) nearest the real axis is at w = i pi/gamma, so the true decay rate is b = pi/gamma in [1,2). Numerically R_infinity(x)exp(b x) is constant for x >= 2 (0.478 at Delta = -0.5, 1.379 at Delta = -0.9, 14.3 at Delta = -0.999). Since only b+beta >= 3/2 is used, beta = 1/2 is enough. The pointwise ratio R_infinity/[C0 exp(-|x|/2)] stayed below 0.25.

The derivative bounds (6.4) follow because differentiating the first kernel of each Neumann term gives at most H_j times the same positive term. The source-to-interval kernel (I-T_Q)^(-1)T(.-y) is bounded by R_infinity(.-y) termwise. The bound on d_Q and its derivatives, the ratio bounds with P_j = (8/3)C0 H_j, and the core radius c are correct.

The global lower bound (6.6) is uniform in Q. For Q < c+1 it uses a1 on a bounded interval. For Q >= c+1 the boundary value r_Q(Q) is fed from a unit interval inside the core, and the boundary layer uses monotonicity of r_Q and cosh(bx)/cosh(bQ) >= exp(-b(Q-x)). The near-free branch |Delta| <= 1/16 uses the pointwise constant (see Issue 2). Numerically min r_Q/r_infinity was 0.77 (Delta = -0.5) and 0.40 (Delta = -0.9), far above kappa.

### Velocity variation and clipping (Sections 6.3 to 6.4)

The identity r_Q'/r_Q-ell = (-d_Q'+ell d_Q)/r_Q and its derivative bound with W are correct, as are |r_Q'| <= H1 r_Q and the formula Z_Q' = Z_Q(Q)[R_Q(x,-Q)-R_Q(x,Q)]. The bound |h|Z_Q(Q) = a r_Q(Q) <= a r_infinity(Q) combines with r_infinity(Q)/r_infinity(x) <= 2exp(-b(Q-|x|)) to give the stated exponential localization. All integrals of exp(-(3/2)(Q-|x|)) over [-Q,Q] are at most 4/3, independent of Q. This gives Var(v_Q) <= Vv uniformly in Q. Sampled Var(v_Q) never exceeded 8.2 at Delta = -0.5 and 3.2 at Delta = -0.9.

For (6.9), the restricted densities are pointwise ordered in Q, so the uncentered cumulative functions differ monotonically by at most |nu-m|. That gives ||G|| <= B/(2L)+B/L. I re-derived the tail masses 4B0/L and 5B0/L, including H(+infinity) = nu-m. Clipping interacts with the Fermi boundary only through A_eff = min(A,Q). If Q <= A, then r_Q >= kappa r_infinity(A) = 2D/L on all of [-Q,Q], so sigma_L >= r_Q/2 there and no jump or tail occurs. If A < Q, the clipped piece lies inside the occupied interval and its cost uses epsilon_Q(Q) = 0 and |epsilon_Q'| <= 2a r_Q. The integral of 1/(L r_Q) over the core is at most 1/(bD) through kappa, so no log L appears in either case. Vclip and (6.12) follow as written.

### Constants (Section 6.5)

Each majorization holds. One small point: sin(alpha0/2) = sqrt(delta/2) exactly, which is what makes 1-Mbeta >= sqrt(delta/2) work. kappa is the third entry, and 1/kappa < 5*10^5 delta^(-19/6). The dominant term is 128 pi D^2/kappa^2 <= 2.88*10^14 pi delta^(-25/3). On 400 log-spaced delta in [1e-6, 1], max K_all(delta) delta^(25/3)/10^15 = 0.0015, which matches the authors' 0.001512.

### Assembly (Section 9)

The field drops out by C129 Step 1, with n_* chosen at the physical flux, and E_L >= F_L(n_*) + field. The doubled ring uses (2L,2n_*) at flux 2Phi, and twist costs JL'(1-cos(pi/L'))/2 with L' = 2L, which gives pi^2 M_0/4 after multiplying by L. M(B) >= M_0 = J/2 on |Delta| <= 1. So L(E_2L-2E_L) <= 5K_all M_0+(pi^2/4)M <= [5K_all+pi^2/4]M <= 6*10^15 delta^(-25/3) M. Spin flip covers n_* > L/2. J = 0 and saturated sectors are handled separately. I also checked the optional high-field placement (8.4) to (8.5). The decomposition K_L|_n = Delta L/4 + tn + (1-t)T_hop + t(T_hop-A) and the bound sin(pi n/L)/[n sin(pi/L)] <= cos(pi/8) < 15/16 for L >= 8 and n/L >= 1/4 are correct.

### Numerical test of the sector theorem

With bethe_scan.solve and an independent composite Gauss-Legendre Nystrom solution of the restricted equation, I computed Lambda = L(F_L(n) - L e_bulk(n/L)) at Delta = -0.5 and -0.9, L in {16,17,32,33,64,65,128,129,256,257,400,401}, and n in {1,2,3, floor(sqrt L), L/8, L/4, 3L/8, L/2-floor(sqrt L), L/2-3, L/2-2, L/2-1, floor(L/2)}. Lambda is negative, varies smoothly in nu, and is bounded uniformly by 0.341 (Delta = -0.5) and 0.1333 (Delta = -0.9). It tends to the half-filled Casimir value -pi v/6 (0.34009 and 0.13324) as nu -> 1/2 at both parities, and to zero in the dilute limit. There is no growth with L near the Fermi-boundary crossover to Q_L = infinity. The other ratios to their proved bounds were at most L|nu-m|/(B-1) = 0.079, L||G||/(3B/2) = 0.24, and doubling defect L[F_2L(2n)-2F_L(n)] = 0.51 J. None of this certifies anything, and the proved constants carry slack of more than ten orders of magnitude.

## Issues

1. **PROOF Section 6.4, last paragraph, exact half filling.** Severity: minor (exposition). "The same clipped identities hold by integrability or by taking the preceding bounds to the full-line limit" is vaguer than needed. Repair: say that at Q = infinity every step of 6.4 holds verbatim with r_Q = r_infinity, Z = B (constant, so Z' = 0), h = 0, epsilon = -a r_infinity, A_eff = A, R_L the full-line midpoint residual with z_L ranging over exactly (-nu/2, nu/2), and nu-m = 0. This is C134 Sections 4 to 7 with the new (larger) constants, and needs no limit.

2. **PROOF Section 6.2, the |Delta| <= 1/16 branch.** Severity: minor. The "pointwise version of C0" is not defined. Repair: define C0(Delta) = 4T(0)/(1-M_beta(Delta)) with T(0) = |Delta|/(pi sin gamma) and M_beta(Delta) = sqrt(2)sin(gamma/2-pi/4), note that (6.3) to (6.5) hold with C0(Delta) in place of C0, and state that the arithmetic 4/[(16)(3)(9/10)(9/10)] uses pi > 3 and 1-M_beta >= 9/10.

3. **PROOF (6.4), source kernel.** Severity: minor. "The same bounds hold for the kernel carrying a source outside [-Q,Q] to a point inside it" should carry its one-line reason. Repair: write the kernel as (I-T_Q)^(-1)T(.-y) = T(.-y)+R_Q T(.-y), bounded termwise by R_infinity(.-y), with x-derivatives falling on the first factor.

4. **PROOF Section 6.1, decay rate.** Severity: remark, not an error. The bound exp(-|x|/2) is valid but the true rate is b = pi/gamma >= 1, from the pole of 1/(1+Khat) at w = i pi/gamma. C0 is also of order 1/delta where the true prefactor scales like 1/alpha, of order delta^(-1/2). A sharper beta would improve the exponent 25/3, if anyone cares.

5. **Equation numbering.** Severity: minor (housekeeping). Section 8 numbers its equations (7.1) to (7.3), and Section 9 mixes (8.1) to (8.5) with (9.1) and (9.2), so equation numbers are offset by one from their sections there. RESULT.md points to "Sections 6 and 9" for the constants, which is right. Repair: renumber Sections 8 and 9 consistently. Also, the merged JSON in this folder records the moment check from the `--moment-only` rerun.

6. **PROOF Sections 5 and 6.4, small-size branch.** Severity: minor. The statement that e_bulk(nu) lies in the bond spectral interval only needs rational nu = n/L. At nu = 1/2 it uses C134, and below 1/2 it uses the large-L branch along multiples of the denominator. Saying so explicitly would avoid the appearance of an appeal to continuity at irrational nu.

7. **RESULT.md, PROOF.md, STATE.md status lines.** Severity: minor (housekeeping). They say the new proof has not had an independent referee pass. Update them after this report.

## Verdict

**Correct as stated.** For every 0 < delta <= 1, Delta in [-1+delta, 0], L >= 2 and 0 <= n <= L/2, the bound |F_L(n) - J L e_bulk(Delta, n/L)| <= 10^15 delta^(-25/3) J/L follows from the argument as written. So does L(E_2L - 2E_L) <= 6*10^15 delta^(-25/3) M on the stated U(1) class, via C129's field, DM and twist reductions. The Neumann inverse on the finite interval is bounded by B uniformly in the Fermi boundary. The resolvent comparison with the sech density holds uniformly in Q, with a valid (non-sharp) exponential decay. The chemical-potential subtraction removes the first-order endpoint error, the curvature bound makes the mismatch quadratic, and the uniform variation bound together with clipping inside the occupied interval leaves no log L. The listed items are expository. The most substantive is Issue 1, the one-line treatment of exact half filling. The theorem is not uniform as delta -> 0, which the authors state themselves.
