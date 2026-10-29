# No universal kappa: the two-site distinguishability D_2 is the trace norm of the bond-gradient of A, and it is unbounded at d = 2

Author: Claude (Anthropic, Fable 5.1), 2026-10-29. Closed book apart from the files named in Section 7. Numerics are numpy and scipy only, in this folder, all run locally (nothing here needed the cluster). Every claim is tagged proved, numerical, heuristic or assumed. Nothing was delegated to Codex. Natural logarithms unless a base is written.

## Claim

Setting as in `../fable_bound_proof/README.md`: a one-site translation-invariant nearest-neighbour chain on (C^d)^L with bond h, infrared a finite set of unitary CFT sectors (c_a, v_a), Casimir coefficient A = sum_a c_a v_a, gauge-optimized bond half-spread M, two-site marginals rho_2^L (ring) and rho_2^inf (line), and

    D_2(h) = liminf_L L^2 || rho_2^inf - rho_2^L ||_1 .

The question was whether a universal kappa with D_2 <= kappa log_2 d exists, or more weakly D_2 <= f(d). The answer is no, already at d = 2 and c = 1.

1. [Proved, under the infrared assumption of Section 0] Gradient identity. For every two-site Hermitian O,

       lim_L L^2 Tr[ (rho_2^inf - rho_2^L) O ] = (pi/6) d/ds A(h + s O) |_{s=0} ,

   where the derivative is taken at fixed charge sectors (fixed zero modes). Hence L^2 (rho_2^inf - rho_2^L) -> (pi/6) grad_h A and

       D_2 = (pi/6) || grad_h A ||_1 .

   Proposition 1 of the previous folder (pi A/(6L^2) = Tr[h~ (rho_2^inf - rho_2^L)]) is the Euler identity for the degree-one homogeneous function A: Tr[h grad_h A] = A. In CFT language, (1/(2 pi)) dv_a/ds along O is the coefficient of (T_a + Tbar_a) in the infrared expansion of O.

2. [Proved, given the Bethe-ansatz velocity] Closed form on the XXZ line h = S^xS^x + S^yS^y + Delta S^zS^z, Delta = cos(gamma), A = v = pi sin(gamma)/(2 gamma), M = 1/2:

       D_2(Delta) = (pi/3) [ |v'(Delta)| + max( |v'(Delta)|, |v - Delta v'| ) ] ,

   which reproduces the two previously known values, D_2(0) = (pi+2)/3 and D_2(1) = pi^2/6, and is minimal at the Heisenberg point. Near the ferromagnetic point, v ~ sqrt((1+Delta)/2) and

       D_2(Delta) ~ (pi sqrt 2 / 6) (1 + Delta)^{-1/2} -> infinity  as Delta -> -1^+ .

   So sup over d = 2 chains of D_2 is infinite, M stays at 1/2, and A -> 0 along the family. The Hölder bound A <= (6/pi) M D_2 becomes infinitely loose: at Delta = -0.99 it reads 0.074 <= 7.7.

3. [Numerical] Exact diagonalization of rings to L = 22 at Delta = 0.5, 0, -0.5, -0.8, -0.9, -0.95 confirms both the gradient identity (the S^zS^z and the xy directions separately, to three or four digits after 1/L^2 extrapolation) and the closed form, including the growth toward the ferromagnetic point (D_2 extrapolates to 3.09 at Delta = -0.9 against 3.097, and to 3.98 at -0.95 against 4.04).

4. [Proved] No h-independent norm can replace the trace norm. If N is any norm on the constraint space W (traceless, equal one-site marginals) with |Tr h~ X| <= M(h) N(X) for all bonds h and all X in W, then N dominates the trace norm on W, so N(grad A) is unbounded on the XXZ family as well. The weighting must depend on h, and the only h-dependent pairing available, Tr[h~ grad A] = A, is the statement to be proved.

5. [Heuristic, perturbation theory, and numerical for decoupled flavours] The unboundedness has nothing to do with many slow sectors. In the campaign's slow-mode constructions a sector whose logical hopping arises at second order through a port with coupling J_a contributes to grad A in proportion to dv_a/dJ_a = 2 v_a/J_a, which is of order sqrt(v_a), so slow sectors are suppressed in D_2 (less than in A, but suppressed). For n decoupled flavours (d = 2^n) D_2 is velocity-blind and grows like roughly 1.5 sqrt(n), not n. The divergence found here is instead a single c = 1 sector whose velocity vanishes as the square root of the distance to a z = 2 gapless point, which makes A continuous but not Lipschitz in the bond.

What this implies for the velocity-weighted conjecture A <= pi M log_2 d: it is a statement about a degree-one homogeneous, gauge-invariant function A(h) that is dominated by the norm pi M(h) log_2 d but whose gradient is unbounded near gapless non-conformal points. No argument that controls A through a norm of its gradient, which is what every two-site Hölder route does, can prove it. What survives is the integrated form, Section 5.

## Status

Items 1 and 4 are proved at the level of first-order perturbation theory (Hellmann-Feynman) plus finite-dimensional duality, under the infrared assumption stated in Section 0, which is the same one the previous folder used plus uniformity in the bond parameter on compact sets. Item 2 is proved given the des Cloizeaux-Gaudin velocity of the XXZ chain and the Yang-Yang ground-state energy; a mean-value variant (Proposition 3) that needs only the Casimir asymptotics, not the interchange of the derivative with the L limit, is also proved. Item 3 is exact diagonalization. Item 5 is heuristic for the port constructions and exact for decoupled flavours. The conjecture itself is untouched by this note except for the negative conclusion about the route.

## 0. Setting and the infrared assumption

[Assumed, same as before] For each bond h in the class, E_0(L) = e_0 L - pi A/(6L) + o(1/L) on vacuum-compatible rings, and local expectation values converge to their line values with corrections of order L^-2 governed by the stress tensors. [Assumed, new] For a smooth one-parameter family h(s) staying in the class, the o(L^-2) remainders in e_0 - E_0(L)/L are uniform in s on compact intervals, together with their first s-derivatives. For the XXZ line this uniformity is what Section 3 checks numerically; the mean-value Proposition 3 avoids the derivative part of it.

[Definition] For a direction O (two-site Hermitian), A_O'(h) = d/ds A(h + sO) at s = 0 at fixed charge sectors: the ring sequence L and the conserved-charge sector of the ring ground state are held fixed while s varies. When O couples to a conserved charge (a chemical potential), the ring ground state at fixed L does not follow a smooth filling change, and the first-order response of the ring energy is exactly the line response (Section 2, remark), so the fixed-sector derivative is the one Hellmann-Feynman produces.

## 1. The gradient identity

[Proved, Proposition 1] Hellmann-Feynman on the ring and on the line: for H(s) = sum_i (h + sO)_{i,i+1},

    d/ds E_0(L; s)/L = Tr[ rho_2^L(s) O ] ,      d/ds e_0(s) = Tr[ rho_2^inf(s) O ] ,

the first exact at every L where the ring ground state is non-degenerate (it is, in the vacuum sector), the second standard for the line (e_0 is concave in s, differentiable wherever the two-site marginal is unique, which is the case for a unique translation-invariant ground state). Therefore

    Tr[ (rho_2^inf - rho_2^L) O ] = d/ds [ e_0 - E_0(L)/L ] = d/ds [ pi A(s)/(6L^2) + r_L(s) ]

with r_L = o(L^-2). Under the uniformity assumption, L^2 d r_L/ds -> 0, so L^2 Tr[(rho_2^inf - rho_2^L) O] -> (pi/6) A_O'(h). Since this holds for every O and the two-site operator space is finite dimensional, L^2 (rho_2^inf - rho_2^L) converges to the operator (pi/6) grad_h A defined by Tr[(grad_h A) O] = A_O'(h), and D_2 = (pi/6) ||grad_h A||_1 with the liminf a limit.

[Proved, Corollary (Euler)] A(lambda h) = lambda A(h) for lambda > 0 (velocities scale, central charges do not), so Tr[h grad_h A] = A. With the gauge directions o (x) 1 - 1 (x) o and the identity lying in the kernel of grad_h A (A is gauge invariant and insensitive to constants), this is Proposition 1 of the previous folder, pi A/(6L^2) = Tr[h~ (rho_2^inf - rho_2^L)] + o(L^-2).

[Proved, Corollary] For any O, D_2 >= (pi/6) |A_O'(h)| / ||O||. This is the lower bound used in Section 2.

[Heuristic, dictionary] In the infrared the lattice operator O expands as <O>_inf + sum_a [t_a(O) T_a + tbar_a(O) Tbar_a] + primaries + higher descendants. On the ring only the identity-module operators have one-point functions, <T_a>_L = <Tbar_a>_L = -pi^2 c_a/(6L^2), so L^2(<O>_inf - <O>_L) -> (pi^2/6) sum_a c_a [t_a(O) + tbar_a(O)]. Adding s sum_i O_i to the Hamiltonian adds s [t_a(O) T_a + tbar_a(O) Tbar_a] to the infrared Hamiltonian density of sector a, whose density is (v_a/2pi)(T_a + Tbar_a), so dv_a/ds = pi [t_a(O) + tbar_a(O)] for the average velocity. The two expressions agree: the gradient identity is the statement that the ring one-point function of the stress tensor is the same object as the velocity response. This is point (ii) of the task: in terms of operator content,

    D_2 = (pi^2/6) || sum_a c_a (Theta_a^T + Theta_a^Tbar) ||_1 ,

where Theta_a^T is the two-site operator dual to the projection onto T_a, Tr[Theta_a^T O] = t_a(O). The two-site algebra having dimension d^2 bounds the number of independent Theta's, not their trace norms, and Section 2 shows a single Theta with unbounded trace norm at d = 2.

## 2. The XXZ line and the ferromagnetic point

[Proved] By U(1), spin flip, reflection and translation, the two-site marginal of the XXZ chain (ring or line, even L, Delta > -1, where the S^z_tot = 0 ground state is unique by Perron-Frobenius after a sublattice rotation) is

    rho_2 = 1/4 + 4 g_xx (S^xS^x + S^yS^y) + 4 g_zz S^zS^z ,      g_xx = <S^xS^x>, g_zz = <S^zS^z>.

Applying Proposition 1 with O = S^zS^z (A_O' = dA/dDelta = v') and with O = S^xS^x + S^yS^y (A_O' = A - Delta v' by Euler),

    L^2 (g_zz^inf - g_zz^L) -> a = (pi/6) v',      L^2 (2g_xx^inf - 2g_xx^L) -> b = (pi/6)(v - Delta v').

The eigenvalues of L^2 (rho_2^inf - rho_2^L) are then a, a (on |00>, |11>), b - a (triplet-0) and -b - a (singlet), so

    D_2(Delta) = 2|a| + |b - a| + |b + a| = 2|a| + 2 max(|a|, |b|).

[Proved, given the Bethe ansatz] With v = pi sin(gamma)/(2 gamma), Delta = cos(gamma),

    v'(Delta) = (pi/2) (sin gamma - gamma cos gamma) / (gamma^2 sin gamma) .

Checks: at Delta = 0, v' = 2/pi, a = 1/3, b = pi/6, D_2 = 2/3 + pi/3 = (pi+2)/3, the free-fermion value of the previous folder. At Delta -> 1, v' -> pi/6, a = pi^2/36, b = pi^2/18, D_2 = pi^2/6, the Heisenberg value of Proposition 4 there. On the whole line |b| >= |a| (so D_2 = 2|a| + 2|b|), the minimum of D_2 is at the Heisenberg point, and D_2 is monotonically increasing toward Delta = -1. As gamma -> pi, sin gamma -> sqrt(2(1+Delta)), v -> sqrt((1+Delta)/2), and v' -> 1/(2 sqrt(2(1+Delta))), so

    D_2 -> (2pi/3) v' = (pi sqrt 2/6) (1+Delta)^{-1/2} .

`xxz_closed_form.py` (output `xxz_closed_form.out`) tabulates this; D_2 sqrt(1+Delta) reaches 0.7412 at 1+Delta = 10^-6 against pi sqrt2/6 = 0.7405. The bond spectrum is {Delta/4, Delta/4, 1/2 - Delta/4, -1/2 - Delta/4}, half-spread 1/2 for all Delta in (-1, 1), so M <= 1/2 and the Hölder bound (6/pi) M D_2 is 2.3, 3.0, 7.7, 23 at Delta = -0.8, -0.9, -0.99, -0.999 while A = 0.38, 0.25, 0.074, 0.023.

[Proved, Proposition 3, mean-value form without interchanging derivative and limit] Let f_L(Delta) = L^2 (e_0 - E_0(L)/L), which by the Casimir assumption converges to (pi/6) v(Delta) uniformly on [Delta_1, Delta_2] inside (-1, 1). Then f_L'(Delta) = L^2 (g_zz^inf - g_zz^L) <= (1/4) L^2 ||rho_2^inf - rho_2^L||_1, and by the mean-value theorem for every L there is a Delta_L in the interval with

    L^2 || rho_2^inf - rho_2^L ||_1 (Delta_L) >= 4 |f_L(Delta_2) - f_L(Delta_1)| / (Delta_2 - Delta_1) -> (2pi/3) (v(Delta_2) - v(Delta_1))/(Delta_2 - Delta_1).

Taking Delta_1 = -1 + eps, Delta_2 = -1 + 2 eps gives a right side of order (pi sqrt2/6)(1+Delta)^{-1/2}. So the scaled distinguishability is unbounded over the XXZ family uniformly in L, using only the Casimir asymptotics; the pointwise statement about D_2(Delta) at each fixed Delta is what needs the derivative interchange, and that is what Section 3 checks.

[Remark on zero modes, proved for free fermions] The divergence requires the velocity to vanish at fixed charge sector. For the XX chain in a field at fixed filling k_F, rho_2 depends on the field only through g = sin(k_F)/pi, and the L^-2 difference is (pi sin k_F/6) times a bounded operator, so D_2 -> 0 in the dilute limit at fixed filling even though dv/dB = -B/sqrt(1-B^2) diverges at saturation; that derivative changes the filling, which the ring at fixed L does not follow (the first-order ring response to a field is exactly -(L/2 - N), the line value). The ferromagnetic XXZ point is different because it is reached at fixed S^z_tot = 0, so the vanishing stiffness shows up in the velocity response rather than in the filling.

[Heuristic, generality] At a generic z = 2 gapless point reached at fixed charge (a ferromagnetic or Lifshitz transition where the spin stiffness vanishes linearly in the control parameter), v^2 is proportional to the stiffness, so v ~ sqrt(delta) and D_2 ~ delta^{-1/2}. The spin-1 bilinear-biquadratic chain approaching its SU(3) ferromagnetic point from the SU(3)_1 side (d = 3, c = 2) should behave the same way. I did not compute it.

## 3. Numerical checks

[Numerical] `xxz_ed_check.py` (output `xxz_ed_check.out`): Lanczos in the S^z_tot = 0 sector on rings L = 8 to 22, with rho_2^inf exact from the Yang-Yang energy e_0(Delta) = Delta/4 - (sin gamma/2) int dx sinh((pi-gamma)x)/(sinh(pi x) cosh(gamma x)) and its Delta-derivative (checks: e_0(0) = -1/pi, e_0(1/2) = -3/8, e_0(1) = 1/4 - ln 2). `xxz_ed_extrapolate.py` fits x(L) = x_inf + a/L^2 + b/L^4 over the four largest sizes:

| Delta | L^2 dg_zz, ED | (pi/6) v' | L^2 d(2g_xx), ED | (pi/6)(v - Delta v') | D_2, ED | D_2 closed | A, ED | v |
|---|---|---|---|---|---|---|---|---|
| 0.5 | 0.2971 | 0.2966 | 0.5317 | 0.5319 | 1.6574 | 1.6569 | 1.2991 | 1.2990 |
| 0.0 | 0.3334 | 0.3333 | 0.5236 | 0.5236 | 1.7139 | 1.7139 | 0.9999 | 1.0000 |
| -0.5 | 0.4142 | 0.4142 | 0.5472 | 0.5472 | 1.9228 | 1.9229 | 0.6495 | 0.6495 |
| -0.8 | 0.5708 | 0.5708 | 0.6543 | 0.6542 | 2.4502 | 2.4499 | 0.3773 | 0.3773 |
| -0.9 | 0.7436 | 0.7448 | 0.8025 | 0.8035 | 3.092 | 3.097 | 0.2545 | 0.2545 |
| -0.95 | 0.973 | 0.989 | 1.016 | 1.031 | 3.98 | 4.04 | 0.1743 | 0.1737 |

Convergence is from above at every Delta, and slower near the ferromagnetic point as expected (the corrections scale with v'/v). The previous folder's unextrapolated value 1.669 at Delta = 0.5, L = 18, was heading to 1.657, consistent with its guess of a value between the Heisenberg and XX numbers; the closed form now fixes it.

[Numerical] `multiflavour_product.py` (output `multiflavour_product.out`): n decoupled XX chains stacked into one site, d = 2^n, bond h = sum_a lambda_a h_a. The ground state is independent of the lambda_a, so D_2 is independent of the velocities. From the exact two-mode Gaussian marginals at L = 2000:

| n | d | D_2 | D_2/n | (pi^2/6) log_2 d |
|---|---|---|---|---|
| 1 | 2 | 1.7139 | 1.714 | 1.645 |
| 2 | 4 | 2.2953 | 1.148 | 3.290 |
| 3 | 8 | 2.5595 | 0.853 | 4.935 |
| 4 | 16 | 3.1993 | 0.800 | 6.580 |
| 5 | 32 | 3.4431 | 0.689 | 8.225 |

The flavour contributions add like independent random signs (trace norm roughly 1.5 sqrt(n)), well inside kappa log_2 d for any kappa >= (pi+2)/3. Decoupled slow modes therefore never threatened a kappa bound; the threat is item 2.

## 4. Point (i): slow modes versus the actual mechanism

[Proved for decoupled flavours] A sector realized by a direct lattice hopping lambda_a h_a has dv_a/dlambda_a = v_a/lambda_a, of order one however small lambda_a is: D_2 does not see the velocity of such a sector. This is the velocity-blindness the task anticipated, and the table above shows it costs at most sqrt(n) at d = 2^n.

[Heuristic, second-order perturbation theory] In the port constructions of `../routeA_final/` and `../small_d_note/`, the logical hopping of flavour a is t_a^eff ~ J_a^2/Delta_a with J_a the port coupling and Delta_a the virtual gap, and v_a is proportional to t_a^eff. The bond-gradient of A along the port direction O_a (the operator multiplying J_a in h) is dA/dJ_a = c_a dv_a/dJ_a = 2 c_a v_a/J_a ~ c_a sqrt(v_a Delta_a). Along the perturbation itself (O = the whole second-order coupling) Hellmann-Feynman gives a response proportional to v_a, which is the Euler identity restricted to that sector. So a slow sector contributes to D_2 at order sqrt(v_a), between its weight v_a in A and the order-one weight it would have if realized at first order. With the geometric velocity schedules of the constructions, sum_a sqrt(v_a) converges and D_2 stays bounded as n grows, while sum_a c_a diverges. A schedule with v_a ~ a^-2 would make D_2 grow like log n at fixed d, but that is a weaker and more model-dependent divergence than item 2, and it was not pursued. Conclusion for (i): many slow modes do not make D_2 large; a single mode whose velocity is non-Lipschitz in the bond does.

## 5. Point (iii): what replaces the trace norm

[Proved, Proposition 4] Let N be a norm on W such that |Tr h~ X| <= M(h) N(X) for every bond h (with any gauge representative h~) and every X in W. By Proposition 3 of the previous folder, for every X in W there is a traceless Hermitian h_X with M(h_X) = 1 and Tr[h_X X] = ||X||_1 (the dual pair is attained in finite dimension). Hence ||X||_1 <= N(X) for every X in W, and N(grad_h A) >= ||grad_h A||_1 = (6/pi) D_2 is unbounded on the XXZ family. A weighted norm that "pairs only with bond-like operators" cannot be defined independently of h: at d = 2 every element of W is the gradient direction of some bond. If the inequality is only demanded for the pairs (h, grad_h A), the weakest admissible weighting is N_h(X) = |Tr h~ X|/M(h), and then N_h(grad A) <= pi log_2 d is the conjecture itself.

[Proved, identity] The integrated form survives. Along any smooth path h(t) in the class,

    A(h(1)) - A(h(0)) = int_0^1 Tr[ h'(t) grad A(h(t)) ] dt <= (6/pi) int_0^1 ||h'(t)|| (pi/6) D_2(h(t)) dt = int_0^1 ||h'(t)|| D_2(h(t)) dt .

Because the singularity of D_2 at a z = 2 point is (1+Delta)^{-1/2}, it is integrable: on the XXZ line int_{-1}^{1} D_2 dDelta = 4.010, and with ||S^zS^z|| = 1/4 and the 6/pi the path bound reads pi/2 = 1.571 <= 1.915 (`xxz_ed_extrapolate.out`). So A is of bounded variation along bond-space paths even where it is not Lipschitz. A d-only bound by this route would need, for every h, a path from a bond with A = 0 (a z = 2 point, or a gapped or trivially decoupled bond) along which int ||h'|| D_2 dt <= C M(h) log d, with D_2 allowed to blow up integrably at the endpoint. I have no argument that such paths exist in general; the statement is recorded as the correct weakening of "D_2 <= kappa log_2 d", not as progress on it.

[Heuristic, channel remark] The trace-norm route works in the Casimir channel, where A = sum c_a v_a. The one lattice inequality that carries a log d is the entropy bound s(T) <= log d, which lives in the thermal channel where the CFT quantity is sum_a c_a/v_a, with velocities inverted; combined with the end of the CFT regime at T ~ v_min it gives only crossover-type bounds (the c <= O(xi_cross log d) of C69). The conjecture is the statement that the modular S transformation of this trivial bound survives on the lattice at order 1/L, and Lemma 1 of the previous folder shows the naive lattice version of S loses 2 beta M. Nothing here changes that assessment.

## 6. Remaining gap

For the question asked, none: a universal kappa does not exist, and D_2 <= f(d) fails for every f at d = 2, given the des Cloizeaux-Gaudin velocity and the standard infrared assumption with uniformity in Delta (checked numerically to L = 22). For the conjecture, the gap is unchanged in content but sharper in form: A(h) is a degree-one homogeneous, gauge-invariant function dominated by the norm pi M(h) log_2 d, with square-root cusps (unbounded gradient) at gapless non-conformal points, so any proof must bound A without passing through a norm of grad A. The two things not done here: a second example of the divergence in an interacting d = 3 chain (the bilinear-biquadratic chain near its SU(3) ferromagnet is the natural one), and a test of the sqrt(v_a) suppression of Section 4 on an explicit port toy, both of which are only confirmations of a conclusion already fixed by the d = 2 result.

## 7. Sources used

`../../specs/_bound_background.md`, `../fable_bound_proof/README.md` and its scripts (conventions, Propositions 1 to 4, the XX value), `../../STATUS.md` (slow-mode constructions, the paragraph on local boundedness near gapless non-conformal points). The XXZ velocity and ground-state energy are the des Cloizeaux-Gaudin and Yang-Yang results, quoted from memory and checked against three exactly known values of e_0.
