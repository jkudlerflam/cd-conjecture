# C141 derivations: the fast-sector budget

Run by a Claude Opus agent, 2026-11-07. Standard: theoretical physics (RG, perturbation theory, Lieb–Robinson, sum rules), every step tagged proved, derived (physics standard, assumptions listed), numerical, heuristic or assumed. Conventions are those of `specs/_bound_background.md`: A = sum_a c_a v_a, M the gauge-optimized bond half-spread, V = 4eM the inherited Lieb–Robinson bound on every sector velocity, N_c(u) = sum_{v_a >= u} c_a, and A = int_0^V N_c(u) du (C133, eq. 10).

## 0. Standing assumptions

- (IR) [assumed] Each fixed bond has finitely many gapless unitary CFT sectors (or is gapped), so A(h) < infinity bond by bond. This is the infrared assumption of every earlier run. Section 7.3 notes what happens if it is dropped.
- (LR) [proved, inherited, C109] v_a <= V = 4eM.
- (S-add) [assumed] Thermal entropy densities of distinct sectors add, and the nonconformal part of the spectrum contributes nonnegatively.

## 1. What any budget must weight

### 1.1 The weight lemma

[proved] If sum_a c_a (v_a/M)^q <= B with 0 <= q < 1, then for every u, N_c(u) (u/M)^q <= sum_{v_a >= u} c_a (v_a/M)^q <= B, so N_c(u) <= B (M/u)^q, and C133's integration gives A <= B (4e)^{1-q} M/(1 - q). A bound of the same form with q = 1 is the target itself. With q > 1 it gives N_c(u) <= B (M/u)^q, which is not integrable at u = 0. So a budget has to weight slow sectors more heavily than A does. A sum rule weighting sectors by v^k with k >= 1 can never supply it, whatever its constant.

### 1.2 Per-sector spectral weights of the energy density

[derived within CFT, exact for one sector] Let eps_q = sum_x e^{-iqx} eps_x be the Fourier mode of the bond-energy density on a ring of length L, q = 2 pi n/L. In a sector with central charge c and velocity v, H = (2 pi v/L)(L_0 + Lbar_0 - c/12) and eps_q = (2 pi v/L)(L_n + Lbar_{-n}). The Virasoro algebra gives <0|L_n L_{-n}|0> = c(n^3 - n)/12, and the only intermediate states sit at omega = v|q|. Hence, per unit length and for small q,

    S_eps(q)   = (1/L)<eps_q eps_{-q}>             = c v^2 |q|^3/(24 pi) (1 - n^{-2}),
    m_1(q)     = (1/L) sum_n omega_n |<n|eps_q|0>|^2 = c v^3 q^4/(24 pi)  (1 - n^{-2}),
    chi_eps(q) = (2/L) sum_n |<n|eps_q|0>|^2/omega_n = c v q^2/(12 pi)    (1 - n^{-2}).

So the frequency moment m_j(q) = int omega^j S(q, omega) weighs sector a by c_a v_a^{2+j}. The static susceptibility (j = -1) weighs by v, as A does, and the q^2 coefficient of chi_eps is at least A/(12 pi) once all sectors are conformal at scale 1/q.

[numerical, xx_chi.py] For the XX chain (c = v = A = 1) the three ratios S/(q^3/24 pi), m_1/(q^4/24 pi), chi/(q^2/12 pi) equal 1 - 1/n^2 to five digits (0.75000 at n = 2, 0.93750 at n = 4, 0.99609 at n = 16, L = 16000). This checks the three formulas including the Virasoro n^3 - n.

### 1.3 Dimension matching: why only charge sum rules are crossover-free

[derived] On the lattice, every nested-commutator sum rule controls an odd positive moment. For the energy density, m_1(q) = (1/2L)<[[eps_q, H], eps_{-q}]> <= 8 M^3 q^2 (C100), and in general m_{2j+1}(q) <= K_j M^{2j+3} q^2. The lattice side is O(q^2) because eps_x is a local operator with no conservation constraint at the lattice scale. The CFT side of m_1 is O(q^4), because eps has scaling dimension 2. A sector therefore enters the inequality only at wavevectors q <~ 1/xi_a, where it is conformal. With that restriction,

    sum_{a: xi_a < 1/q} c_a v_a^3 <= 192 pi M^3/q^2,

which is a v^3/xi^2 weighting. It has the wrong power of v (section 1.1, k = 3) and carries the crossover length. This is the precise reason for C100's M^3/T_*^2.

[derived] For a conserved U(1) density rho_q the CFT side is m_1 = K v q^2/(2 pi) (from S_rho = K|q|/(2 pi) at omega = v|q|) and the lattice f-sum is also q^2 times a ground-state expectation. Both sides scale as q^2, so the comparison holds at q -> 0 with no crossover, and with D_a = v_a K_a/pi it is the Kohn budget of section 9.75, sum_alpha D_alpha <= 4M, i.e. sum_a K_a v_a <= 4 pi M. It is and it is the only crossover-free sum rule of this type, because only spin-one currents have scaling dimension equal to their engineering dimension. The criterion is: a sum rule is crossover-free if and only if the probed density is conserved and its IR scaling dimension is 1. Energy density (dimension 2), order parameters and neutral primaries fail it.

## 2. Candidate (a): the ground-state energy

[derived] For free fermions the budget A <= pi|e_hop| is a Kohn budget, not an energy-lowering statement. Each Fermi point pair is a sector with K = 1 exactly, its Drude weight is v/pi, and the optical f-sum gives sum D <= -<kinetic>, with every quantity exact at every scale. C106's avoided crossing does not hurt it, because the f-sum never asks the sector to be conformal.

[derived scaling, heuristic in the O(1) constant] For an interacting sector, the ground-energy density attributable to it is the zero-point energy of modes it reorganizes below its crossover, roughly

    delta e_a ~ -c_a int_0^{1/xi_a} (dk/2 pi) [sqrt(v_a^2 k^2 + (v_a/xi_a)^2) - v_a k] ~ -c_a v_a/xi_a^2,

with energy density in units of velocity/length^2. Any bound that uses the ground energy (e_0 >= -||h||, a variational comparison with a gapped competitor, or a sum over sectors of energy lowered) therefore controls sum_a c_a v_a/xi_a^2, and a fast sector with a long crossover costs essentially no energy. The natural generalization of A <= pi|e_hop| is the Kohn budget sum K_a v_a <= 4 pi M of section 9.75. It fails exactly for small Luttinger weights and for neutral sectors (section 1.3). Candidate (a) yields nothing beyond what section 9.75 already proved.

## 3. Candidate (b): the zero-temperature energy-current f-sum

[derived, negative] At T = 0, m_1(q) for the energy density is a ground-state expectation of a local double commutator, bounded by 8 M^3 q^2. By section 1.3 its CFT content is sum c_a v_a^3 q^4/(24 pi), so the weight is v^3, and the q^2 mismatch is filled only for xi_a < 1/q. Both defects are structural. A Cauchy–Schwarz transfer to lower moments goes the wrong way: m_0^2 <= m_{-1} m_1 bounds the v-weighted m_{-1} (which carries A) from below, never from above. The v-weighted quantity chi_eps is a second-order response to a modulated bond. It is bounded below by A/(12 pi) q^2 and above only through the near-critical massive modes of lesson 2 (C111: a gapped spectator of gap 2 delta contributes about 1/delta). Moments with j <= -2 weigh by v^{<= 0}, which would be the right direction, but no commutator controls them. m_{-2} is velocity-blind (weight c), and the c = n + 1 constructions show that velocity-blind counts are unbounded.

## 4. Candidate (c): entanglement generation

[derived] With small incremental entangling (rate <= 18 M log d per cut) and the boundary-quench formula dS/dt = pi A/(12 ell), valid only once ell exceeds every crossover, A <= (216/pi) ell M log d. Restricted to sectors with xi_a <= ell it gives sum_{xi_a <= ell} c_a v_a <= (216/pi) ell M log d. Even if crossovers were bounded, this is a q = 1 statement (N_c(u) <= C ell M log d/u), the borderline case of section 1.1. Entanglement rates weigh sectors by v, so they give p = 1 at best, and the crossover enters through ell. This confirms the expectation of the spec.

## 5. The compactness reduction and the emergence lemma (framework for Tasks 1 and 3)

### 5.1 Compactness

[proved] Let B_d be the set of bonds modulo scalar shifts, gauge h -> h + o (x) 1 - 1 (x) o and on-site unitaries, normalized to M = 1. It is compact. Its generic dimension is d^4 - 1 - (d^2 - 1) - 1 - (d^2 - 1) = d^4 - 2d^2, which is 8 at d = 2 (the 8-parameter family of STATUS). If a function F: B_d -> [0, infinity] is locally bounded at every point (each h0 has a neighbourhood U on which sup F < infinity), then sup F < infinity by Heine–Borel.

[derived consequence] Under (IR), sup_h A(h)/M < infinity at fixed d if and only if A is locally bounded on B_d. Any family h_N with A(h_N)/M -> infinity has a subsequence converging to some h* in B_d, and A is not locally bounded at h*. So every hidden tower is, at the top level, a perturbative problem around a single accumulation bond h*, with h_N = h* + delta_N and delta_N -> 0. (C131 is the rigorous version of this for each fixed L. Here it is used for the asymptotic A.)

### 5.2 Regular bonds

Call h0 regular if its IR is gapped or a finite CFT reached by a finite crossover length xi0.

[derived; assumes perturbative RG and a c-theorem for weak perturbations of a multi-velocity CFT] For |delta| small, the flow of h0 + delta agrees with that of h0 up to scales ~ xi0, then is CFT_0 perturbed by couplings O(delta). Its IR has c_tot <= c_tot(h0) (no new sectors from a weakly perturbed CFT), and velocities <= V. So A(h0 + delta) <= c_tot(h0) V on a neighbourhood. Regular bonds are points of local boundedness, whatever their crossover lengths, and long crossovers that are continuations of a regular limit bond (C106, where the limit epsilon -> 0 is a regular two-Dirac-cone bond with c = 2) are harmless.

[assumed, flagged] The c-theorem step is standard for Lorentz-invariant flows. For flows with several velocities coupled by relevant cross terms it is a physical expectation, not a theorem. Hybridization through marginal couplings can lift slow velocities above a threshold u, but by eigenvalue interlacing a rank-r coupling moves at most r velocities across any threshold, so it cannot create an unbounded number of fast sectors from a bounded number of couplings.

### 5.3 Singular bonds with a gapped degenerate manifold: the emergence lemma

Call h* singular if its low-energy spectrum contains a nonconformal part: an extensively degenerate ground manifold D, a z > 1 band, or an unbounded hierarchy (section 7.3).

[derived; assumptions: D separated by a gap Delta* from the rest of the spectrum of h*, and |delta| << Delta*] Schrieffer–Wolff gives an effective Hamiltonian on D,

    H_eff = P (delta H') P + O(delta^2/Delta*),

which is quasi-local, with terms of range r carrying norm <~ |delta| (|delta|/Delta*)^{r-1}. The Lieb–Robinson velocity of such a Hamiltonian is bounded by C(r_0) |delta| M (1 + O(delta/Delta*)). Leakage out of D is O(delta/Delta*). Every sector that emerges from D (one not continuously connected to a sector of h*) therefore has

    v_em <= C |delta| M.                                                   (5.1)

At k-th order emergence (first k - 1 orders vanish on D or are classical), v_em <= C_k |delta|^k M^{1-k}.

[derived consequence] Near such an h*, A(h* + delta) <= A_cont + N_em(delta) C |delta| M, where A_cont <= c_cont V comes from the conformal part of h* (section 5.2). Local boundedness at h* is equivalent to N_em(delta) |delta| bounded. Write H_eff = |delta| h_eff, with h_eff in the compact class of normalized quasi-local effective bonds on D. Then N_em(delta) = c_tot(h_eff) and the problem recurses: local boundedness at h* follows from local boundedness of c_tot (not only of A) in the effective class, or from A(h_eff) bounded together with the factor |delta|. Each level of recursion multiplies velocities by a factor O(delta_j/Delta_{j-1}) < 1.

[numerical, fib_emergence.py] At the d = 2 Fibonacci point h* = n (x) n (D = no two adjacent up spins, Delta* = 1) perturbed by delta(s+ s- + h.c.), the emergent sector's velocity, measured from the lowest momentum 2 pi/L excitation in the ground-state particle sector, satisfies v/delta -> v_eff(projected model) linearly in delta. For example, at L = 20: 2.7427, 2.7171, 2.6142 for delta = 0.01, 0.03, 0.1 against v_eff = 2.7544. The thermodynamic value from the exact map to free fermions on L - N sites is v_eff = 2 sin(pi rho')/(1 - rho) = 2.792 at the optimal density rho = 0.3008. (Even L/odd N* rings, L = 16 and 18, pick a different finite-size excitation at momentum 2 pi/L and read about twice this. The comparison of full and projected models is internally consistent there too.)

### 5.4 Singular bonds with a z > 1 band

[derived scaling] If the nonconformal part of h* is a z > 1 band (ferromagnetic magnons, Lifshitz points), a perturbation delta opens a Fermi sea of size k_F ~ delta^{1/z} and gives v ~ M k_F^{z-1}, so v_em -> 0 as delta -> 0. The number of emergent sectors is bounded by the number of z > 1 band minima of h*, which is finite. Local boundedness holds. (The U(1) d = 2 case is C136/C138: v ~ sqrt(1 + Delta) at the ferromagnetic point.) Emergent fast sectors with long crossover cannot come from z > 1 regimes, because v xi^{z-1} ~ M forces xi <~ (M/v)^{1/(z-1)}.

### 5.5 Where a fast sector with long crossover can come from

[derived from 5.2 to 5.4] Combining the three cases: a sector with v ~ M that becomes conformal only at a length xi -> infinity along a family must be either (i) a continuation of a sector of a regular limit bond (harmless, section 5.2), or (ii) emergent from a gapped degenerate manifold, which (5.1) forbids because v_em = O(delta), or (iii) emergent from a z > 1 band, which section 5.4 forbids, or (iv) created by a multi-velocity z = 1 flow that increases the fast central charge, which is a failure of the assumed c-theorem of section 5.2. The fast, late-emerging sector of the question is therefore excluded except through (iv) or through a failure of (IR) (section 7.3). Unbounded A without fast sectors (C133's v_k = M/k tower) remains possible in principle, through unbounded N_em at fixed level. Section 7 treats it.

## 6. Task 1: d = 2

### 6.1 The singular strata at d = 2

[derived, completeness heuristic] Singular bonds in B_2 are of three types.

1. Classical (commuting) bonds with an extensive ground space. Up to a common on-site rotation and gauge, a classical bond is h = J Z Z + b (Z (x) 1 + 1 (x) Z). Its four configuration energies are J + 2b, -J, -J, J - 2b. Extensive degeneracy occurs only for J > 0 and b = +-J, where the ground space is the Fibonacci shift (no two adjacent up spins, or down spins), with entropy log phi. At J < 0, b = 0 the degeneracy is twofold, not extensive. So there is one classical singular stratum, the Fibonacci point and its images. (The fully degenerate h = 0 is excluded by M = 1.)
2. Frustration-free gapless bonds. By Bravyi and Gosset (arXiv:1503.04035, J. Math. Phys. 56, 061902), frustration-free nearest-neighbour qubit chains are gapped unless the projector is onto a maximally entangled state, which up to an on-site unitary is the ferromagnetic Heisenberg (twisted) case with z = 2 magnons. Saturation points in a field (the polarized product state degenerate with the one-magnon band minimum) are of this z = 2 type.
3. Quadratic Lifshitz points (XY at gamma = 0, |h| = 1), also z = 2. A one-band quadratic qubit chain has no zero-energy flat band unless h = 0.

[heuristic] No other nonconformal gapless or extensively degenerate qubit bonds are known, and the campaign's survey (about 103,000 sampled bonds, C6_d2) found none. First-order coexistence between two gapped phases is not singular in this sense, because domain walls cost finite energy at h*.

### 6.2 Local boundedness at the d = 2 strata

[derived; single-species input] At each stratum the emergent degrees of freedom consist of one species: up spins in the Fibonacci space, or magnons over a product state (the one-magnon space of a nearest-neighbour qubit chain has one state per site, so one cosine band). Under U(1)-preserving perturbations, a one-species lattice gas with nearest-neighbour interactions has c <= 1. Bound pairs cannot coexist with singles as a second fluid, because a cluster of n adjacent particles has energy (n - 1)V: any attraction strong enough to bind pairs binds larger clusters more strongly per particle, and the gas phase-separates. In the U(1) class this is confirmed by the Bethe ansatz (C134 to C139). Under U(1)-breaking perturbations the one-species gas has p-wave pairing, which is relevant at the free point and gives a gap or an Ising transition (c = 1/2). The Fibonacci emergent model at first order is the hard-rod gas with constrained hopping and transverse field (PXP-type, Fendley–Sengupta–Sachdev), which has c <= 1 in both cases. In all strata v_em = O(delta^alpha M) with alpha = 1 (Fibonacci, by (5.1)) or 1/2 (z = 2).

[derived; nesting depth] A new singular stratum inside an effective class at level j requires the level-j effective bond to be singular, an independent fine-tuning condition. For example, P (Z Z) P vanishes on the Fibonacci space up to a chemical potential, so the first-order classical effective model has no new degenerate point, and the next one (the sub-shift of no two particles within distance 2) needs a second-order next-nearest-neighbour term tuned against a first-order chemical potential. With 8 parameters in B_2, nesting depth is at most 8 if each level imposes at least one independent condition [heuristic: transversality of the conditions is assumed]. At d = 2 the conformal part of each singular stratum is empty, so emergent levels do not coexist. A one-species level that is conformal has no further degenerate manifold, and one that is degenerate has no conformal part.

### 6.3 Result for d = 2

[derived, conditional] Under (IR), the weak-perturbation c-theorem of section 5.2, the completeness of the singular list of section 6.1 and the single-species analysis of section 6.2, A is locally bounded at every point of the compact B_2, hence

    C(2)_asymptotic = sup_h A(h)/M < infinity.

With the additional non-coexistence statement of section 6.2 and the absence of c > 1 at regular bonds, c_tot <= 1 and A <= V = 4eM. Compactness alone gives finiteness without a value, so the second statement needs more than the first: a strongly coupled regular qubit bond with c = 3/2 (Ising plus Luttinger) is not excluded by compactness, only by its absence from the search and by the single-species heuristic. No such bond is known. The nearest known c > 1 qubit models (interacting Majorana chains with c = 3/2 phases, O'Brien–Fendley) need range-3 terms or two-site unit cells, which is d = 4 after blocking.

[open] The step with the least support is completeness of the singular list at strong coupling. A qubit bond whose IR were a nonconformal gapless state of a new kind (neither degenerate nor z > 1) would require a separate analysis.

## 7. Candidate (d), Task 3: rank, second-order emergence and hidden towers

### 7.1 First order and the rank lemma

[proved, C133] An exact first-order quadratic product-block compression carries at most floor((d^2 - 1)/2) complex sectors, by operator-Schmidt rank of one physical cut and the degree of the Fermi determinant. The rank of an effective coupling across one cut is bounded by d^2 - 1 only at first order.

### 7.2 Second order: a singular-value budget

Setup [assumed model]: cells of n physical sites with a quadratic, number-conserving logical model, inter-cell hopping matrix T generated at second order through intermediate states of the left and right cells with energies E_i, F_j in [a, b], a > 0 the smallest intermediate excitation, b <~ M. Through nu port operators (nu <= d^2 - 1, the operator-Schmidt rank of the cut bond), T = sum_{p=1}^{nu} G_p Ct_p G'_p, with Ct_p(i, j) = x_i y_j/(E_i + F_j), visibilities x, y and contractions G, G'.

[proved, cited] Ct_p has displacement rank one (diag(E) Ct + Ct diag(F) = x y^T), so by Beckermann and Townsend (SIAM J. Matrix Anal. Appl. 38, 1227, 2017; arXiv:1609.09494) its singular values obey s_{k+1} <= Z_k s_1 with Z_k <= 4 exp(-pi^2 k/log(4b/a)). [numerical, cauchy_sv.py] This bound holds with ratio at most 0.38 on 16 test matrices (box and log-uniform spectra, random and end-visibility weights, a/b from 1e-2 to 1e-6).

[derived] Summing, ||Ct_p||_1 <= s_1 [1 + 4/(exp(pi^2/log(4b/a)) - 1)] <= s_1 [1 + (4/pi^2) log(4b/a)]. The quadratic trace lemma of section 9.71 gives, in cell units, A_cell = (1/2) sum_FP |v| <= 2 ||T||_1, and one cell is n sites, so

    A_em <= 2 n nu s_1 [1 + (4/pi^2) log(4b/a)].                            (7.1)

Second-order perturbation theory requires port coupling lambda <= eta a with eta < 1, and s_1 <= lambda^2/(2a) sum x_i^2 <= eta^2 a/2 by the visibility sum rule sum_i x_i^2 <= 1. Therefore

    A_em <= nu eta^2 n a [1 + (4/pi^2) log(4b/a)].                          (7.2)

This is a crossover-free budget at second order. It bounds A, and so N_c(u), directly. Its constant is governed by n a, the cell length times the smallest intermediate gap, and by the logarithm of the cell's dynamic range.

[derived consequences]
- Gapped cell interior (a ~ Delta_int fixed): logical modes that must cross the cell see visibility e^{-n/xi_int}, s_1 is exponentially small and (7.1) is tiny. Logical modes delocalized over the cell have wall amplitude ~ n^{-1/2}, so s_1 ~ lambda^2/(n a), the n cancels, and A_em <= nu eta^2 a (1 + 0.41 log(4b/a)) <= C nu b. That is bounded.
- Critical cell interior with velocity v_int: a <= pi v_int/n, so n a <= pi v_int and (7.2) gives A_em <= nu eta^2 pi v_int [1 + (4/pi^2) log(4 b n/(pi v_int))]. This grows at most like log n, and only if wall-localized logical modes have O(1) visibility to every intermediate level.
- The route-A encoder (C133) sits far inside (7.2): its payload velocities carry D^{-2n-4} and V_n <= (n+2)^{-2}.

So at second order a hidden tower is possible only logarithmically in the cell length, and only with critical cell interiors that are cut into finite pieces by walls while staying visible at the walls. The Cauchy example of C133 (full rank at second order) does not defeat the budget: rank is unbounded but singular values decay geometrically, at a rate set by log(b/a).

[open] Order k >= 3 gives products of k - 1 resolvents. The displacement structure is then of rank growing with k, and an analogous Zolotarev bound with nu^{k-1} ports is expected but was not derived. The quadratic logical model is also an assumption. Interacting or neutral logical payloads need a nonquadratic replacement for the trace lemma.

### 7.3 The hierarchical tower of C133 at the physics level

[derived] Take sector k with v_k ~ M conformal beyond b^k, k = 1, ..., N, along a family h_N. By section 5.1, h_N -> h* along a subsequence. Two branches:

(i) For fixed k, sector k's crossover stays bounded as N -> infinity. Then h* itself carries all sectors with v_k ~ M and crossovers b^k, so c(h*) = infinity and (IR) fails at a single bond. Its ring energy would be E_L = e L - (pi/(6L)) sum_{b^k < L} v_k, so L(E_{2L} - 2E_L) ~ (M/log b) log L, which already violates the finite-size C(d) at one fixed bond. Its ground-state entanglement would grow like (log l)^2. Entropy capacity does not forbid it (C133, eq. 9). The Casini–Huerta entropic c-theorem does forbid it if the flow at every scale is Lorentz invariant with a single velocity. So branch (i) needs the not-yet-conformal parts at every scale to break Lorentz invariance while staying z = 1 (z > 1 parts produce only slow sectors, section 5.4). That is a multi-velocity z = 1 regime that increases the fast central charge, case (iv) of section 5.5.

(ii) Sector k's crossover diverges with N at fixed k. Then sector k emerges at h* from a nonconformal part. By (5.1) and section 5.4, emergence from a gapped degenerate manifold or a z > 1 band gives v_k -> 0, which contradicts v_k ~ M. Branch (ii) also reduces to case (iv).

[derived conclusion] At the physics level the equal-velocity hierarchical tower is equivalent to a z = 1, Lorentz-violating RG flow that increases the central charge carried at velocities >= u. Within the perturbative mechanisms available at a fixed accumulation bond (Schrieffer–Wolff emergence, z > 1 bands, finite-rank hybridization) the fast count cannot grow. The dangerous object that remains is the slower harmonic tower (v_k ~ M/k), realized through many sectors at a single emergence level. At second order with a quadratic payload, (7.2) caps it at logarithmic growth in the cell length.

### 7.4 What the budget candidates say about the towers

- (a) and (b) see a tower sector only through xi_k^{-2} (energy) or q <= 1/xi_k (f-sum). Both are blind to it.
- (c) gives p = 1 and needs ell > b^N. Blind.
- (d) at first order caps the number of sectors at (d^2 - 1)/2. At second order, (7.2) caps A_em by n a log(b/a), which is finite for every fixed cell structure and at most logarithmic in the cell length. This is the only candidate that touches the towers.
- The Kohn budget charges every charged tower sector v_k K_k. An equal-velocity tower with bounded Luttinger weights is impossible (sum K_k v_k <= 4 pi M), so it must be neutral or have K_k -> 0.

## 8. Conclusion

[derived, conditional] C(2) < infinity in its asymptotic form follows at the physics standard from compactness of B_2, the weak-perturbation c-theorem at regular bonds, the Schrieffer–Wolff/Lieb–Robinson emergence lemma (5.1), the z > 1 scaling of section 5.4, and the classification and single-species analysis of the d = 2 singular strata. With non-coexistence it sharpens to c_tot <= 1 and A <= 4eM. The weakest link is completeness of the singular strata at strong coupling.

[derived, conditional] For general d the same reduction shows that C(d) < infinity follows from two statements: (P1) a multi-velocity c-theorem, that a z = 1 flow does not increase the central charge carried above any velocity u beyond a finite-rank hybridization; and (P2) a single-level emergence budget, that at each singular bond the emergent A stays bounded as delta -> 0. (P2) is derived at second order for quadratic payloads in the form (7.2), with a constant that grows at most like log(b/a). It is open at higher orders and for interacting or neutral payloads.

[derived, negative] No sum rule can be the budget. Commutator sum rules weigh sectors by v^{>= 3} with a crossover (section 1.3), energy-based bounds weigh by v/xi^2 (section 2), and entanglement rates weigh by v with a crossover (section 4). The only crossover-free sum rule is the Kohn budget, which weighs by K v and is blind to neutral and weakly charged sectors. A budget with the right weight (v^q, q < 1) has to be structural, as (7.2) is, not spectral.

[bearing on the rigorous program] The compactness reduction says that a rigorous C(d) needs local bounds only near singular bonds. At d = 2 these are the three explicit strata of section 6.1. Two of them (the ferromagnetic one and the quadratic Lifshitz one) are already largely covered by C129, C132, C135 to C138. The Fibonacci stratum is new to the rigorous list. Its first-order model maps exactly to free fermions on L - N sites, so a Lemma-P or window-type argument plus a Schrieffer–Wolff error should close it.

## 9. Consistency checks

| Case | Expected | What the derivations give |
| :-- | :-- | :-- |
| Free fermions | A <= 2M, A <= pi|e_hop| | Kohn budget with K = 1 (section 2); first-order rank lemma; (7.1) with T the bare hopping and n = 1 reduces to A <= 2||T||_1, the section 9.71 trace lemma |
| XXZ | A/M = pi sin gamma/gamma | regular bond, c = 1, A <= V; Kohn: K v = pi/(2(pi - gamma)) * pi sin gamma/(2 gamma) <= 4 pi M |
| Heisenberg | A = pi M, the maximum at d = 2 | regular (marginal crossover only), c_tot = 1 <= 1, A = pi M <= 4eM |
| SU(N) chains | A = 2 pi (N-1) M/N | regular; covered by the WZW stiffness budget (C104); no emergence |
| C106 | v/M fixed, crossover diverging, A/M = 2r | continuation of a regular limit bond (epsilon = 0, c = 2 Dirac point), harmless by section 5.2 |
| c = n + 1 encoder | c unbounded, A bounded | one emergence level with n sectors and v ~ D^{-2n-4}, far inside (7.2) |
| Fibonacci point | v_em = O(delta) | ED: v/delta -> v_eff linearly in delta (section 5.3) |
| XX energy density | per-sector weights | xx_chi.py, ratios 1 - 1/n^2 to five digits |
