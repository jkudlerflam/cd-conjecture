# cd_free_route: derivations

Conventions of `specs/_bound_background.md`. A bond h is a Hermitian operator on C^d (x) C^d, H_L = sum_x h_{x,x+1} on the ring Z_L, E_L its ground energy, e_0 = lim E_L/L the energy density on Z, and M(h) the gauge-optimized half-spread, M(h) = min over real a and Hermitian o of the half-spread (lambda_max - lambda_min)/2 of h + a 1 + o (x) 1 - 1 (x) o.

## 1. The Casimir variant is a statement about quantum marginals (proved)

Let V be the real vector space of Hermitian operators on C^d (x) C^d and G the gauge subspace spanned by 1 and by o (x) 1 - 1 (x) o. M is a seminorm on V that vanishes exactly on G, hence a norm on V/G.

Define two compact convex sets of two-site density matrices.

- Omega_L, for L >= 2: the nearest-neighbour marginals rho_{x,x+1} of translation-invariant states on the ring (C^d)^{(x) L}.
- Omega_inf: the nearest-neighbour marginals of translation-invariant states on the chain over Z.

Every element of either set has equal one-site marginals, so Tr[g rho] - Tr[g sigma] = 0 for all g in G and rho, sigma in either set. Hence differences rho - sigma lie in the annihilator of G, which is the dual space of V/G, with dual norm ||Delta||_* = sup{Tr(h Delta) : M(h) <= 1}. Since |Tr(h Delta)| <= M ||Delta||_1 after choosing the optimal gauge representative of h centred at zero, ||Delta||_* <= ||Delta||_1.

Lemma 1. E_L(h)/L = min over rho in Omega_L of Tr(h rho) and e_0(h) = min over rho in Omega_inf of Tr(h rho).

Proof. Averaging any ground state over the ring translations gives a translation-invariant ground state, whose energy is L Tr(h rho_{x,x+1}). The infinite-volume statement is the standard variational characterization of the mean energy, with the minimum attained by weak-* compactness of translation-invariant states. []

Proposition 1 (equivalence). For a constant K >= 0 and a length L the following are equivalent.

(i) For every bond h on C^d (x) C^d, E_L(h) >= e_0(h) L - K M(h)/L.

(ii) Omega_L is contained in Omega_inf + (K/L^2) B_*, with B_* the unit ball of ||.||_*. In words, every nearest-neighbour marginal of a translation-invariant L-ring state is within dual distance K/L^2 of the nearest-neighbour marginal of a translation-invariant state on Z.

Proof. By Lemma 1, (i) says that the support function of Omega_L in the direction -h is bounded by that of Omega_inf plus (K/L^2) M(h) = (K/L^2) times the support function of B_* in the direction -h, for all h. Support functions add under Minkowski sums, and for compact convex sets an inequality of support functions in every direction is equivalent to inclusion (Hahn–Banach). Both sets lie in the affine slice where the gauge pairing is constant, so directions in G are trivial on both sides. []

Consequences.

1. Variant (a) of the target with a d-dependent K(d) is equivalent to (ii) for all L with K = K(d). The statement is kinematic: no Hamiltonian appears.
2. Classical states. For diagonal (commuting, classical) two-site marginals, every element of Omega_L is a marginal of the stationary Markov chain with transition matrix p(s'|s) = p(s, s')/p(s), so Omega_L^cl is contained in Omega_inf^cl and K = 0, consistent with C = 0 for classical chains.
3. Quantum obstruction to K = 0. Omega_2 contains the singlet projector, whose distance to Omega_inf is of order one by monogamy, so K >= L^2 times a positive constant at L = 2, harmless for finite K but showing that the inclusion is genuinely quantum.
4. A proof of (ii) is a periodization construction: from a translation-invariant ring state build a translation-invariant state on Z whose nearest-neighbour marginal moves by at most K/L^2 in the dual norm. Its energy for any h then exceeds E_L/L by at most K M/L^2, which is variant (a). The construction may be universal (independent of h), which is what the rest of this note attempts.

## 2. Quasi-free chains: momentum box smoothing (proved)

Take a translation-invariant quasi-free state on the L ring of a chain with m fermion modes per site, written in Nambu form, with symbol Ghat(k), 0 <= Ghat(k) <= 1, at the allowed momenta k_j = (2 pi/L)(j + phi). The offset phi in [0, 1) covers every Jordan–Wigner parity sector, so the construction applies to the spin ring ground state of any quadratic spin chain. Define a symbol on the circle by Ghat_Z(k) = Ghat(k_j) for k in [k_j - pi/L, k_j + pi/L). It obeys 0 <= Ghat_Z <= 1 pointwise, so it is the symbol of a translation-invariant quasi-free state on Z. Its real-space correlations are

    G_Z(r) = (1/L) sum_j Ghat(k_j) e^{i k_j r} sinc(pi r/L) = G_L(r) sinc(pi r/L),   sinc(u) = sin(u)/u.

So G_Z(0) = G_L(0) and G_Z(+-1) = G_L(+-1) sinc(pi/L). For a quadratic nearest-neighbour bond write the ring energy per site as E_L/L = e_on + e_hop with e_on the on-site (r = 0) part and e_hop the r = +-1 part, both evaluated in the ring ground state. The Z state has energy per site e_on + sinc(pi/L) e_hop >= e_0. Hence, with 1 - sinc(u) <= u^2/6,

    E_L >= e_0 L - L [1 - sinc(pi/L)] |e_hop| >= e_0 L - (pi^2/(6L)) |e_hop|.        (2.1)

Consistency check against XX (h = Sx Sx + Sy Sy, e_on = 0, e_0 = -1/pi). (2.1) reads E_L >= -csc(pi/L), and ED gives equality on every even L from 2 to 12 (run_checks.py, check 1). Since |e_hop| tends to |e_0| = 1/pi, variant (a) gives A <= (6/pi)(pi^2/6)(1/pi) = 1, which is the exact XX value A = c v = 1. The box smoothing is therefore sharp for XX. In general (2.1) gives A <= pi |e_hop^inf| for every quasi-free chain, with or without pairing, which is of the same strength as the known A <= 2M.

Interpretation. The Z state is the direct integral over twists theta in [-pi/L, pi/L) of the ring state shifted in momentum by theta. For interacting chains there is no such decomposition of the many-body space, which is why a different, state-level construction is needed (section 3).

## 3. A universal periodization channel (proved)

Let rho be a translation-invariant state on the L ring, l an integer with 2 <= l <= L - 1, rho_l its marginal on l consecutive sites, and rho_{l-1} its marginal on l - 1 consecutive sites. Label a window of l sites as A (first site) and B (the other l - 1 sites), and let C be one further site. The Petz map

    R(X) = rho_l^{1/2} (rho_{l-1}^{-1/2} X rho_{l-1}^{-1/2} (x) 1_C) rho_l^{1/2}

is completely positive and trace preserving on operators supported in supp(rho_{l-1}), and R(rho_{l-1}) = rho_l.

Construction. On the half-line {1, 2, ...} set sigma^(l) = rho_l on sites 1..l, and for n >= l set sigma^(n+1) = (id_{1..n-l+1} (x) R_{n-l+2..n -> n-l+2..n+1})(sigma^(n)), with R acting on the last l - 1 sites and the new site n + 1.

Lemma 3.1. (i) The marginal of sigma^(n) on its last l sites is rho_l for every n >= l. (ii) The marginal of sigma^(n) on sites (k, k+1) is independent of n once n >= k + l, and equals the (A, B_1) marginal sigma_2* of (id_A (x) R)(rho_l) for every k >= 1. (iii) A weak-* limit point of the averages of translates of sigma^(n) is a translation-invariant state on Z with nearest-neighbour marginal sigma_2*. Hence e_0(h) <= Tr[h sigma_2*] for every bond h.

Proof. (i) By induction: tracing the first site of sigma^(n+1)'s last l + 1 sites commutes with R, and gives R(marginal of sigma^(n) on its last l - 1 sites) = R(rho_{l-1}) = rho_l, using the marginal of rho_l on its last l - 1 sites equals rho_{l-1} by translation invariance. The support condition holds at every step for the same reason. (ii) The step n -> n + 1 acts only on sites n - l + 2, ..., n + 1, so the bond (k, k+1) is not touched once k + 1 < n - l + 2. It is last touched at the step n = k + l - 1, whose input has marginal rho_l on sites k..k+l-1 by (i), and the bond's marginal after that step is the (A, B_1) marginal of (id_A (x) R)(rho_l). (iii) The states on finite windows form a compact set; average sigma^(n) over translations by 0..N - 1 on windows far from the boundary and take a limit point. Every nearest-neighbour marginal is sigma_2*, and the variational characterization of e_0 (Lemma 1) gives the inequality. []

So the construction is a universal (h-independent) map Omega_L -> Omega_inf at the level of the full ring state, and by Proposition 1 the variant (a) bound with constant K follows from the key inequality

    delta_L := Tr[h sigma_2*] - E_L/L <= K M / L^2,      at some l = l(L),            (3.1)

evaluated on the translation-invariant ground state (or ground-space mixture) of H_L.

## 4. Reduction of (3.1) to a recovery error (proved)

Let rho_{l+1} be the marginal of the ring state on l + 1 consecutive sites, labelled A B C, with l + 1 <= L, and X = (id_A (x) R)(rho_AB), rho_AB = rho_l. Then Tr_{BC} X = rho_A, and Tr_A X = R(rho_{l-1}) = rho_l, so X and rho_{l+1} have the same one-site marginals on A and on B_1, and the same trace. Hence for every gauge representative h~ of h the gauge terms and the identity drop out, and choosing h~ with spectrum in [-M, M],

    delta_L = Tr[(h~_{A B_1} (x) 1) (X - rho_{l+1})],     |delta_L| <= M || X - rho_{l+1} ||_1.        (4.1)

Note that delta_L >= e_0 - E_L/L is automatic, so (4.1) is a two-sided handle on the finite-size energy.

## 5. Why an information-theoretic closure is impossible (proved, with one standard input)

Write eps = (1/2)||X - rho_{l+1}||_1. Since X = (id_A (x) R)(rho_AB) with R acting on B only, data processing gives S(A|BC)_X >= S(A|B)_rho. The Alicki–Fannes–Winter continuity bound for conditional entropy with |A| = d gives |S(A|BC)_X - S(A|BC)_rho| <= 2 eps log d + (1 + eps) h_2(eps/(1 + eps)). Combining,

    I(A:C|B)_rho <= 2 eps log d + (1 + eps) h_2(eps/(1 + eps)).                      (5.1)

So the trace-norm recovery error is bounded below by the conditional mutual information of the two end sites of an (l+1)-block of the ring state. Within CFT this CMI is the negative second difference of the block entropy, (c/3)(pi/L)^2 / sin^2(pi l/L) at leading order, which at l = L/2 is (c/3) pi^2/L^2 [heuristic, CFT; numerically 3.36/L^2 on Heisenberg at L = 14 against 3.29]. It sees the total central charge c, not A = sum c_a v_a. Inverting (5.1), the trace-norm error obeys eps >= (c pi^2/3)/(L^2 [2 log d + 2 log(L^2) + O(1)]) at leading order, so any trace-norm bound eps <= f(d)/L^2 would force c <= (3/pi^2) f(d) [2 log d + 4 log L + O(1)] on every chain in its CFT regime. Since c is unbounded at fixed d (STATUS.md, conditional constructions with c = n + 1 whose extra sectors are slow), the constant obtained through (4.1) would grow with c up to a logarithm, and at fixed L it cannot be bounded in terms of d alone. Bounds of that unweighted kind give nothing beyond the trivial A <= 4 e M c from Lieb–Robinson. Hence the trace-norm step (4.1) loses the velocity weighting, and a proof of (3.1) must be energy-weighted: it must use h itself, not only the state. This is a precise obstruction for the information-theoretic half of the route, and it is lesson 3 in a new form (slow sectors flood an unweighted quantity).

## 6. Velocity weighting of the construction itself (proved for decoupled stacks)

For a decoupled stack h = h_1 (x) 1 + 1 (x) h_2 on C^{d_1 d_2}, with nondegenerate ring ground spaces or equal mixtures over them, every marginal factorizes, so rho_l^{1/2}, rho_{l-1}^{-1/2} and the Petz map factorize, sigma_2* = sigma_2*^(1) (x) sigma_2*^(2), and delta(h) = delta(h_1) + delta(h_2). Since the ring ground state of eps h_2 does not depend on eps > 0, delta(eps h_2) = eps delta(h_2). So a slow copy enters the defect weighted by its own bond scale, while it enters the CMI of section 5 with its full central charge (run_checks.py, check 7, to 1e-9). The obstruction of section 5 is therefore an obstruction to the trace-norm proof, not to the construction.

## 7. Numerical tests of (3.1) and consistency checks

All values are L^2 delta_L / M at the best l (l = L/2 for pure ground states), from run_checks.py and the logged runs in logs/.

| Chain | L | L^2 delta_L/M | pi A/(3M) |
| :--- | :--- | :--- | :--- |
| XX | 6, 10, 14, 20 | 2.534, 2.167, 2.125, 2.109 | 2.094 |
| XXZ, Delta = 0.5 | 6, 10, 14 | 3.320, 2.825, 2.766 | 2.721 |
| Heisenberg | 6, 10, 14, 20 | 4.054, 3.441, 3.365, 3.333 | 3.290 |
| XXZ, Delta = 2 (gapped) | 6, 10, 14 | 3.286, 2.535, 2.197 | 0 |
| SU(3) swap | 3, 6, 9, 12 | 9.000, 7.494, 5.387, 4.934 | 4.387 |
| FM Heisenberg, multiplet mixture | 4, 8, 12 | 1.333, 1.143, 1.091 (l = L - 1) | 0 |
| TFIM g = 1.6 | 6, 10, 14 | 0.306, 0.083, 0.019 | 0 |
| C106 two-band, eps = 0.03 | 4, 6, 8 | 3.313, 2.530, 2.246 (l <= 4) | |
| fast/slow ladder, eps = 0.01, J = 0.5 | 4, 6, 8 | 3.063, 1.841, 1.365 (l <= 4) | |

[Numerical] Every sequence is bounded and decreasing on the sizes reached, so no loss grows with L. Massive modes are exponentially suppressed (TFIM), and the long-crossover C106 chain and the coupled fast/slow ladder stay uniformly bounded as eps decreases.

[Numerical, heuristic CFT identity] On the critical chains the limit appears to be L^2 delta_L -> pi A/3, that is e_trial = e_0 + pi A/(6L^2) + o(L^-2): the Petz Z state overshoots e_0 by exactly the Casimir deficit. XX at L = 20 gives 2.1086 against 2.0944 and Heisenberg 3.3326 against 3.2899, both still decreasing. If (3.1) is proved with constant K, variant (a) gives A <= (6/pi) K M, and on CFT chains the best possible K is pi A/(3M), so the route loses at most a factor 2.

Consistency checks. XX: section 2 reproduces E_L = -csc(pi/L) exactly on even L and A = 1; the Petz route gives L^2 delta/M -> 2pi/3, twice the exact Casimir coefficient pi A/(6M) = pi/3, and e_trial >= e_0 = -1/pi throughout. Heisenberg: L^2 (e_0 - E_L/L) = 0.840, 0.835 at L = 10, 12 against pi^2/12 = 0.822, and the Petz defect tends to pi^2/3 = 4 times pi^2/12, again a factor 2 in the Casimir coefficient. On Heisenberg the trace-norm bound (4.1) is satisfied with room, L^2 delta = 1.695 against M L^2 ||X - rho||_1 = 2.369 at L = 12.

## 8. Status of the route

[Proved] Sections 1 to 6.

[Assumed sufficient, open] Lemma P. There are K(d) and l(L) such that for every bond h on C^d (x) C^d and every L >= 2, the translation-invariant ground state (or ground-space mixture) of H_L satisfies

    Tr[h~_{A B_1} ((id_A (x) R_{B->BC})(rho_AB) - rho_ABC)] <= K(d) M(h) / L^2.

Lemma P implies E_L >= e_0 L - K(d) M/L for all L and A <= (6/pi) K(d) M. A d = 2 version would give C(2)-type control in the Casimir form.

[Proved obstruction] Lemma P cannot be obtained from the trace norm of the recovery error, by (5.1) and the unbounded central charge at fixed d. Any proof must be energy-weighted, using that rho is the ground state of the same h. The natural energy-weighted mechanism, heuristic, is that the Petz map acts near the A–B boundary as an imaginary modular flow by the mismatch K_BC - K_B of the block modular Hamiltonians, which in CFT is a weighted energy density with weight of order x^2/L^2 and inverse velocity, so its effect on the bond energy scales as v/L^2 per sector. Making this lattice-rigorous needs locality of block modular Hamiltonians, which is not available for critical chains.
