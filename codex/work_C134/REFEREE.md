# Referee report on the C134 compact theorem (Delta in [-1+delta, 0])

2026-11-07. Independent referee pass on PROOF.md Sections 1-8, with RESULT.md, CHECKS.md and LITERATURE.md read for context. The Delta = 1/2 section was not re-checked (already verified separately), and the Delta = 0 constant was checked only by hand (Section 9 below). Check script: `referee_checks.py`, output `referee_checks.json`.

## Summary

I went through Sections 2 through 8 line by line and re-derived every inequality that enters the constant. I did not find a mathematical error, a circular step, a hidden log L, or an L-dependent constant. The issues below are expository. None of them changes a hypothesis, a constant, or a conclusion.

## What I checked, item by item

### 1. Identification of the ground state (Section 2)

Flux. With the sublattice gauge e^{i pi x S^z_x}, H with +XY hopping and total flux L pi becomes the periodic stoquastic K_L at both parities. The diamagnetic inequality then puts the flux minimum at K_L periodic. I confirmed this independently by dense ED of H itself, scanning the flux on one bond (181 points), at Delta = -0.5 and -0.9, L = 6, 7, 8, 9, sector n = floor(L/2). The minimum over flux agrees with the Bethe energy of (2.2) to 5e-15 at both parities.

Yang-Yang functional. The Hessian formula is correct. The off-diagonal entries are +theta'(x_j - x_k), and the quadratic form is L sum p' u_j^2 - sum_{j<k} theta'(x_j - x_k)(u_j - u_k)^2, which is positive definite because theta' <= 0 for pi/2 <= gamma < pi. I recomputed the recession slope on the ray moving the top k coordinates. The sum of the top k quantum numbers is k(n-k)/2, which gives k(pi-gamma)(L-2n+2k) > 0. The rearrangement step and the deduction that the unique minimizer is sorted are correct. The reflection symmetry x_j = -x_{n+1-j} is correct.

Coordinate Bethe ansatz and continuation. The +2Delta contact condition is the one for the -XY hopping of K_L. The dispersion e = -Delta - cos p is right, and the identity e = -pi sin(gamma) a1 holds (I verified it by hand). At gamma = pi/2 the momenta 2 pi I_j / L are the n momenta of the free sector ground state with the correct Jordan-Wigner boundary sign (-1)^{n-1}, and |p_j + p_k| <= 2 pi (n-2)/L < pi, so the denominators are nonzero. The continuation argument does not assume completeness. It uses only (a) that one explicit vector is an eigenvector on [pi/2, pi/2 + eps), (b) the identity theorem for the real-analytic function det[K_L(cos gamma)|_n - E_L(gamma)] on [pi/2, pi), and (c) Perron-Frobenius simplicity of the sector ground state at every gamma (the exchange graph is connected for 1 <= n <= L-1 and the hopping is strictly negative). The agreement set {gamma : E_L = E_0} is closed by continuity and open because the ground eigenvalue is isolated. Odd L needs no separate treatment. This step is sound and not circular.

### 2. Kernel sign and norm (Section 3)

In the chosen parametrization the kernel is

    K(x) = theta'(x)/(2 pi) = sin(2 gamma) / [2 pi (sinh(x)^2 + sin(gamma)^2)],    Delta = cos gamma.

On pi/2 <= gamma < pi, sin(2 gamma) <= 0, so K <= 0. Its integral is [theta(inf) - theta(-inf)]/(2 pi) = 2(pi - 2 gamma)/(2 pi), so ||K||_1 = 2 gamma/pi - 1 = q < 1. Numerical quadrature reproduces q to 1e-15 at Delta = -0.999, -0.99, -0.9, -0.5, -0.1, 0, and the maximum of K on a fine grid is <= 0 in every case. The Fourier transforms and the ratio 1/[2 cosh(gamma w/2)] check out, so rho = 1/[2 gamma cosh(pi x/gamma)] and the integral of rho is 1/2.

The Neumann bound depends only on q, so it holds uniformly in L. In fact one does not need the series: from R = G + K*G, sup|G| <= sup|R| + q sup|G|. The positivity sigma >= a1 > 0 holds at every finite L because -K*mu >= 0. The finite counting function is therefore a global increasing bijection onto (-c, c), with the roots at I_j/L.

### 3. Root, quadrature and clipping estimates (Sections 4-7)

I re-derived each step.

- 2c - nu = (1 - gamma/pi)(1 - 2nu). At odd L, w = (1 - gamma/pi)/L <= 1/(2L).
- ||R||_inf <= 1/(2L) + w/2. The extreme value at x -> +inf is -w, which also lies inside this range.
- H = -K*G and R = (I + K*)G, with no boundary term even though G(+inf) = -1/(2L) at odd L, because the integrated kernel vanishes at -infinity.
- ||K'||_1 = 2|K(0)| = 2|cot gamma|/pi. The pointwise bound |K''| <= |K|(8/s^4 + 4/s^2) follows from d >= s^2 cosh^2, |d'| <= 2 cosh^2, |d''| <= 4 cosh^2.
- The exact identity (5.1) follows from e = (I + K*)eps, the evenness of K, and (I + K*)(mu - rho) = mu - sigma.
- The cutoff gives sigma >= rho/2 on [-A, A]. The bound on |(rho/sigma)'| is 4(bD + D1)/(L rho). The integral of 1/(L rho) over [-A, A] is at most 1/(bD), so the growing cutoff A ~ (gamma/pi) log L gives no log L. This is the key point, and it is right. Var(f') <= 12ab + 4aD1/D.
- The Peano kernel of the one-cell midpoint rule has sup h^2/8, which gives (7.1).
- Tail masses: the total masses are nu <= 1/2 and (1 - gamma/pi) + q nu = nu + w <= 1/2. Combined with |G|, |H| <= B/L and qB/L, they give (7.3). Then (7.4) follows from |eps - g| <= 2aD/L.
- At odd L the strip term has the favorable sign, because f <= 0.

Every error term is (constant depending on gamma)/L. Nothing hides an L dependence.

### 4. Doubling and bookkeeping (Section 8)

2L is always even and uses n = L, so its two-sided estimate applies with the same half-filled e_0. F_{2L} - 2F_L <= K/(2L) + 2K/L = 5K/(2L). This is correct at both parities of L. The small-L branch is also correct: the bond eigenvalues are Delta/4 (twice) and -Delta/4 +- 1/2, so the spectral diameter is 1 on [-1, 1], which gives |F_L - L e_0| <= L < (4 gamma D + 1)^2/L when L < 4 gamma D.

### 5. Constants

Each majorization in Section 8 holds:

- gamma <= pi
- sin gamma >= s = sqrt(delta(2 - delta))
- pi - gamma >= alpha
- |cot gamma|/(pi - gamma) <= 1/(alpha s)
- q <= 1
- 1 < b <= 2
- a <= pi
- D >= 1
- alpha >= sqrt(2 delta), so Bbar <= pi/(2 sqrt(2 delta)) < 9/(8 sqrt delta)
- Dbar <= 1 + 1/(sqrt 2 delta) <= 2/delta
- D1bar <= 27/(2 delta^{5/2})

The bracket is 3 + 27/4 + 128 + 18 = 623/4, and (8 pi + 1)^2 + 623 pi/4 is about 1172.2 < 1200. Numerically, on 400 log-spaced delta in [1e-4, 1], max Kdelta delta^3/1200 = 0.670 (attained at delta = 1). On a 60-point Delta grid inside each [-1+delta, 0], max K(Delta)/Kdelta = 0.373. So 1200/delta^3 and 3000/delta^3 follow as written.

### 6. Numerical stress test of the intermediate inequalities

`referee_checks.py` evaluates the intermediate inequalities with the bethe_scan.py solver at Delta = -0.5, -0.9, -0.99 and L = 8, 9, 16, 17, 64, 65, 128, 129, 256, 257, 512, 513, with doubled rings through 1026. The results, as the largest observed ratio to the proved bound:

| quantity | -0.5 | -0.9 | -0.99 |
| :--- | ---: | ---: | ---: |
| L sup\|G\| / B | 0.505 | 0.232 | 0.081 |
| L sup\|sigma - rho\| / D | 0.032 | 0.051 | 0.033 |
| L sup\|sigma' - rho'\| / D1 | 0.008 | 0.0015 | 5e-5 |
| sampled Var(f') / bound | 0.139 | 0.010 | (no cutoff regime) |
| L \|midpoint error\| / bound | 0.022 | 0.0021 | (no cutoff regime) |
| L \|tail term\| / bound | 0.015 | 0.0017 | (no cutoff regime) |
| L(L e_0 - F_L) / K | 0.0017 | 4e-5 | 1e-7 |

At Delta = -0.99, 4 gamma D is far above 513, so only the trivial branch is active there. The split of the identity (5.1) into its clipped and tail pieces closes to 2.4e-14. For a sanity check, the even-L Casimir term L(L e_0 - F_L) converges to pi v/6 with v = pi sin(gamma)/(2 gamma): 0.340088 against the CFT value 0.340087 at Delta = -0.5, and 0.133243 against 0.133245 at Delta = -0.9. So e_infinity is the right bulk energy. None of these checks is a certified enclosure. The proved bounds carry slack of two to seven orders of magnitude.

## Issues

Each issue lists its location, severity, and a suggested repair.

1. **PROOF Section 2, coercivity.** Severity: minor (exposition). The text says every sorted vector is a nonnegative sum of the "top-k right" and "bottom-k left" rays, and concludes positivity from the ray slopes. The recession function L(pi-gamma) sum|u_j| + (2gamma-pi) sum_{j<k}|u_j-u_k| - 2pi sum I_j u_j is only piecewise linear on the sorted cone (the term |u_j| depends on sign), so additivity over an arbitrary decomposition fails. It does hold for the sign-adapted decomposition, with right rays for the positive coordinates and left rays for the negative ones. Repair: say "use the decomposition adapted to the sign pattern of u, on which the recession function is additive."

2. **PROOF Section 2, multiple contacts.** Severity: minor (citation). The sentence "Multiple contacts are handled by the same scalar exchange factors, which commute" is the standard coordinate Bethe ansatz consistency for spin-1/2 XXZ. It is true, but it is asserted rather than shown. Repair: cite a standard derivation (Yang and Yang 1966, or Alcaraz, Barber, Batchelor, Baxter and Quispel, J. Phys. A 20 (1987) 6397), or add the two-line argument that a cluster's missing-hopping terms decompose into adjacent-pair contact terms.

3. **PROOF Section 2, ground-branch continuation.** Severity: minor (exposition). "A first departure would require a degeneracy" is correct, but the topological argument should be explicit. Repair: say that the set where E_L(gamma) = E_0(gamma) is nonempty (at pi/2), closed (by continuity), and open (the simple ground eigenvalue is isolated, and eigenvalues move continuously), so it is all of [pi/2, pi).

4. **PROOF Sections 3 and 7, definition of e_0.** Severity: minor. The proof sets e_0 := e_infinity and proves lim F_L/L = e_infinity, which is fine if e_0 means the half-filled bulk density. If e_0 means the unconstrained infinite-chain ground energy density, then e_0 <= e_infinity, and the lower bound F_L >= L e_0 - K/L still follows a fortiori. The even two-sided statement would then need the standard fact that the zero-magnetization sector carries the bulk minimum for -1 < Delta <= 0. The doubling bound is unaffected, since e_0 cancels there. Repair: add one sentence fixing the definition.

5. **PROOF Section 7, midpoint rule at the clipping kinks.** Severity: minor. f' jumps at t = z(+-A). The Peano argument works if the expansion about each midpoint uses either one-sided value of f' (the linear term cancels by symmetry) and the jump is counted once in |df'|. Repair: one sentence.

6. **PROOF Section 7, "endpoint atoms may be included by one-sided limits".** Severity: minor. The left-tail mass uses G(-A^-). The bounds (4.1) and (4.2) are sup bounds, so they also bound left limits. Repair: state this.

7. **PROOF Section 8 and RESULT, the delta exponent.** Severity: none (remark, not an error). Kdelta grows more slowly than delta^{-3}. Kdelta delta^3 decreases from 804 at delta = 1 to 0.04 at delta = 1e-4, and the dominant small-delta term is D1bar ~ delta^{-5/2}. Using 1/delta^3 is valid but coarse. One could state Kdelta <= c/delta^{5/2} with a modest c if a sharper statement is wanted.

8. **PROOF Section 1, L = 2 convention.** Severity: minor. "For all L >= 2" depends on C129's two-oriented-bond convention at L = 2. The Bethe energy -Delta/2 - 1 matches K_2 = 2 h_bond in sector n = 1. Repair: restate the convention in the theorem statement in RESULT.md, since a reader using a single bond at L = 2 would get a different F_2. The inequality still holds there via the trivial branch, but the definitions should match.

9. **RESULT.md, PROOF.md, STATE.md status lines.** Severity: minor (housekeeping). These say that no independent referee has checked the argument. Repair: update them after this pass.

## Section 9 (Delta = 0) spot check

L tan(pi/(4L)) sec(pi/(2L)) equals 4 - 2 sqrt 2 at L = 2 and decreases in L. The odd half-filled energy -1/[2 sin(pi/(2L))] is exactly half the n = L energy at 2L, so the odd defect is zero. Both are correct.

## Verdict

**Correct as stated.** The theorem F_L >= L e_0 - 1200 J/(delta^3 L) for all L >= 2, the even two-sided version, and L(F_{2L} - 2F_L) <= 3000 J/delta^3 on Delta in [-1+delta, 0] all follow from the argument as written. The ground-state identification does not assume completeness. The kernel has the stated sign and norm, uniformly in L. The cutoff removes the potential log L. The constants follow from the intermediate bounds. The listed items are expository repairs, the most substantive being the sign-adapted ray decomposition in the coercivity step (item 1) and an explicit definition of e_0 (item 4).
