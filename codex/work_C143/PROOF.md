# C143 derivations: is the list of singular qubit bonds complete?

Run by a Claude Opus agent, 2026-11-07 (not by Codex). Standard: theoretical physics. Every claim is tagged [proved] (a complete argument, or a cited theorem whose hypotheses were checked), [derived] (physics standard, assumptions named), [numerical] (sizes and scripts given), [heuristic] or [assumed]. Conventions follow `specs/_bound_background.md` and `specs/_d2_background.md`:

    h = sum_{ab} J_ab s^a (x) s^b + sum_a B_a (s^a (x) 1 + 1 (x) s^a)/2,

with Pauli matrices s^a, J a real 3 x 3 matrix and B the uniform field (per site, H_L contains sum_j B.s_j). Gauge terms o.(s (x) 1 - 1 (x) s) and constants cancel on rings. B_2 is the 8-dimensional compact space of such bonds modulo scale and simultaneous on-site rotations. M is the gauge-optimized half-spread of C141.

## 0. Summary

The list of C141 PROOF 6.1 is not complete as stated. Three corrections and one new stratum.

1. Stratum 2 (frustration-free gapless) is larger than "the twisted ferromagnetic Heisenberg case". By Bravyi and Gosset it is the two-parameter set of rank-one projectors onto psi with |mu_1| = |mu_2| != 0 (eigenvalues of T_psi), which includes non-maximally entangled psi, plus a rank-two family. It contains a curve of XYZ bonds in a uniform field (J_x = -J_y, at the Kurmann-Thomas-Muller field), which is a staggered-U(1) class not covered by the Bethe-ansatz runs C129 to C140. All of it has z = 2 (Gosset-Mozgunov). [proved, section 2.5]
2. Stratum 3 (quadratic Lifshitz) is not a set of isolated points. Band-onset (Pokrovsky-Talapov, commensurate-incommensurate) transitions are hypersurfaces of codimension one in B_2, generic even without any symmetry, wherever the gap minimum sits at an incommensurate momentum. In the free locus they are explicit (section 2.4). In the U(1) class they include the lower critical field h_c1 of XXZ at Delta > 1. [proved on the free locus, literature-proved in U(1), derived elsewhere]
3. A new singular point: the XZ bond h = X (x) Z (images under on-site rotations, and its reflection Z (x) X). It is the uniform form of the isotropic quantum compass chain. Its spectrum is a critical Majorana chain (c = 1/2) tensored with an exactly degenerate manifold of dimension 2^(L/2 - 1). It is singular in the sense of Task 1, but of neither structural type C141 assumes: the degenerate manifold is not separated by a gap from the rest of the spectrum, and the conformal part is not empty, so C141's emergence lemma (5.3) and its non-coexistence statement (6.2) do not apply as written. [proved, sections 2.2, 2.4, 4.2]
4. The XZ point is isolated: it is the only non-classical bond with a strictly local integral of motion of range two, and the only flat-band point of the free locus. Strictly local integrals of motion of any range require rank J <= 1. [proved, section 2.2]

Its neighbourhood is controlled. The two-parameter family h = XZ + alpha ZX + beta YY is exactly an anisotropic Ashkin-Teller chain written as two commuting Majorana chains with velocity ratio alpha and energy-energy coupling beta. The emergent sector is a second Majorana chain with velocity O(alpha + beta <a>) at first order, so v_em = O(delta M) and c_tot = 1 at generic directions. On the tuned cone where the first-order hopping cancels, the emergent chain is strongly interacting, and its central charge is not settled analytically. Small-ring numerics there do not show c > 1 but are not converged (section 4.2).

## 1. Task 1: an operational criterion for "singular"

### 1.1 Definitions

For a bond h and ring length L let E_0 <= E_1 <= ... be the spectrum of H_L, and define the low-state count and its smooth version

    N_h(kappa, L) = #{n : E_n - E_0 <= kappa M / L},
    S_L(h; kappa) = sum_n exp(-L (E_n - E_0)/(kappa M)).

Definition. h is regular if, along vacuum-compatible lengths, sup_L N_h(kappa, L) < infinity for every kappa > 0 (equivalently S_L bounded in L for every kappa). h is singular otherwise.

[derived] Behaviour of the criterion in each case.

1. Gapped, q-fold (quasi-)degenerate ground space: N -> q for any fixed kappa. Regular.
2. Finitely many CFT sectors with velocities v_a > 0, reached at a finite crossover length: for L >> xi, E_n - E_0 = (2 pi/L) sum_a v_a (x_n^(a) + ...), so N(kappa, L) tends to the number of CFT levels with sum_a v_a x^(a) <= kappa M/(2 pi). This is finite (each tower has finitely many levels below any bound). S_L tends to a product of CFT characters at a fixed modular parameter tau_a = i v_a/(kappa M). Regular.
3. A z > 1 band with dispersion e(k) ~ M k^z near its minimum: single-particle levels within kappa M/L number about L (kappa/L)^(1/z)/(2 pi), which grows like L^(1 - 1/z). N grows without bound. Singular.
4. Extensively degenerate manifold (exact or with splitting o(1/L)): N >= exp(s L). Singular.
5. A CFT coexisting with an extensive manifold (the XZ point): N >= exp(s L). Singular, of a mixed type.

### 1.2 What C141 actually uses

[derived, from C141 PROOF 5.2 to 6.2] Lemma 5.2 (regular bonds are harmless) uses exactly the negation of "singular" above, in the stronger form that beyond a finite xi_0 every state below kappa M/L is described by finitely many CFT towers. The singular side of the argument uses more than the criterion. It assumes that every singular bond is of one of two structural types:

- (G) a gapped, extensively degenerate manifold D separated by a gap Delta* from everything else (Schrieffer-Wolff in 5.3 needs the gap), or
- (Z) a nonconformal part made of finitely many z > 1 band minima (5.4),

and, at d = 2, that the conformal part of every singular bond is empty (6.2, non-coexistence). So a completeness test has to check two things: that N is bounded away from the listed strata, and that at every singular point the structure is (G) or (Z) with no coexisting conformal sector. The XZ point fails the second test.

### 1.3 Practical pitfalls (found in this run)

[numerical] Incommensurate (spiral or floating) bonds have ground-state level crossings at particular ring lengths. A single-increment indicator log S_L2 - log S_L1 is maximized by them (version 1 of the search found such bonds at 7 of 9 starts, all with an accidental fourfold near-degeneracy at L = 12 only, probe_v1.out). The usable finite-size indicator is the minimum over consecutive increments,

    g(h) = min_{L in {L1, L1+2, ...}} [log S_{L+2} - log S_L],

which requires growth at every step. Calibration at L = 8, 10, 12 and kappa = 3 (calibrate_v2_8_12.out): Fibonacci 0.96, XZ 0.69, Bravyi-Gosset real psi (p = 0.2) 0.47, FM Heisenberg 0.33, XX saturation 0.12, gapped Neel 0.04, XX, AF Heisenberg, critical Ising and the Ising ferromagnet within 0.01 of zero. The z = 2 onset signal is weak at these sizes (0.12), which limits what small-ring searches can resolve.

## 2. Task 2: analytic classification by family

### 2.1 Organization of B_2

The symmetry classes used below are: classical (commuting) bonds; the free locus F (quadratic after a Jordan-Wigner map along some axis); frustration-free bonds; the uniform U(1) class (XXZ, field and DM along one axis); the staggered U(1) class (conserved sum_j (-1)^j s^z_j, see 2.7); XYZ in a uniform field along a principal axis; and the residual strata of `_d2_background.md` (one unitary Z2, antiunitary only, no symmetry). The first five admit exact statements.

### 2.2 Strictly local integrals of motion

Definition. A strictly local integral of motion (IOM) of range r is an operator Q supported on r consecutive sites, acting nontrivially on both end sites, with [H_L, Q] = 0 on a ring with L >= r + 2. (Translation-invariant sums such as the Heisenberg charges are not of this type.)

Theorem 2.2. [proved] (a) If h has a strictly local IOM of some range r >= 1, then rank J <= 1. (b) For r <= 2 the IOMs exist exactly in two cases: the classical bonds J = lambda n n^T with B parallel to n (any r), and the XZ point J = lambda m n^T with m orthogonal to n and B = 0, where Q = (n.s)_j (m.s)_{j+1}.

Proof of (a). Write h_{j-1,j} = sum_{mu,nu} h_{mu nu} s^mu_{j-1} s^nu_j with mu, nu in {0, x, y, z}. If Q is supported on [j, j + r - 1], the only term of [H_L, Q] that acts on site j - 1 is [h_{j-1,j}, Q]. Its component along s^mu_{j-1}, mu != 0, is s^mu_{j-1} (x) [R_mu, Q] with R_mu = sum_nu h_{mu nu} s^nu acting on site j, and these components are linearly independent, so [R_mu, Q] = 0 for mu = x, y, z. If the R_mu generate the full matrix algebra M_2 on site j, Q acts as the identity on site j, contradicting the support assumption. So the traceless parts of R_x, R_y, R_z, which are sum_b J_{ab} s^b, are all proportional to one n.s: every row of J is parallel to n. The same argument at the right end (site j + r) shows every column of J is parallel to one vector m. Hence J = lambda m n^T.

Proof of (b). Let J = lambda m n^T with M = m.s, N = n.s, and u = B/2 + o, w = B/2 - o the one-body coefficients of the left and right sites of the bond (o the gauge). By (a) applied at both ends, Q = alpha N_j + beta M_{j+1} + gamma N_j M_{j+1}. The terms h_{j-1,j} and h_{j+1,j+2} commute with Q except through their one-body parts, and collecting, the condition [H_L, Q] = 0 reduces to

    [B.s_j + B.s_{j+1} + lambda M_j N_{j+1}, alpha N_j + beta M_{j+1} + gamma N_j M_{j+1}] = 0

(the gauge o cancels, as it must). Put c = m.n and v = m x n, so MN = c + i v.s. Using [M_j N_{j+1}, N_j M_{j+1}] = 2 i c (V_j - V_{j+1}), V = v.s, the condition splits into one-site and two-site parts:

    one-site:  lambda gamma c v + alpha (B x n) = 0,   -lambda gamma c v + beta (B x m) = 0,
    two-site:  lambda alpha v n^T - lambda beta m v^T + gamma [(B x n) m^T + n (B x m)^T] = 0.

Case m orthogonal to n (c = 0). In the orthonormal frame (m, n, v), B x n = b_m v - b_v m and B x m = -b_n v + b_v n, so the two-site equation reads lambda alpha v n^T - lambda beta m v^T + gamma (b_m v m^T - b_v m m^T - b_n n v^T + b_v n n^T) = 0. The six matrices are independent, so alpha = beta = 0 and gamma B = 0. A nonzero IOM therefore exists iff B = 0, and it is Q = N_j M_{j+1}. This is the XZ point.

Case m parallel to n (c = +-1, v = 0). The equations give alpha (B x n) = beta (B x n) = 0 and gamma [(B x n) n^T + n (B x n)^T] = 0, so an IOM exists iff B is parallel to n: the classical bonds.

Case 0 < |c| < 1. Expand the two-site equation in the (non-orthogonal) basis m, n, v of R^3; the nine dyads are independent. The coefficients of v n^T and m v^T give alpha = beta = 0. If gamma != 0, the remaining coefficients force B x n parallel to n and B x m parallel to m, hence B x n = B x m = 0 and B = 0 (m, n independent). Then the one-site equation gives lambda gamma c v = 0, impossible. No IOM. QED.

[numerical] On the rank-one family (m = x, n = cos t x + sin t z, any B with |B| < 3), a Nelder-Mead minimization of the smallest singular value of Q -> [H_7, Q] over range-3 operators found only the classical points (ioms_rank1_r3.out). The singular value is zero to machine precision at XZ and at classical bonds (sanity checks in the same file). So the range-two classification appears to hold at range three as well.

Consequence. [proved] Exact extensive degeneracy produced by local symmetries (the mechanism of both known examples) occurs in B_2 only at the classical Fibonacci points and at the XZ point. [heuristic] Exact extensive ground degeneracy without strictly local IOMs would need a frustration-free constraint (classified in 2.5, giving only Fibonacci) or a fine-tuned fragmentation mechanism for which no nearest-neighbour qubit example is known.

### 2.3 Classical bonds

[proved] h is classical (some gauge representative has [h_12, h_23] = 0) iff J = lambda n n^T and B is parallel to n. Indeed the three-body part of [h_12, h_23] is gauge independent and equals sum_{a,d} s^a (x) (2i (row_a(J) x col_d(J)).s) (x) s^d, which vanishes iff every row is parallel to every column, i.e. J = lambda n n^T. The two-body parts then require the one-body coefficients on both sides of the middle site to be parallel to n. Up to rotation h = J ZZ + b (Z (x) 1 + 1 (x) Z) (C141's normalization, b = B/2). The configuration energies J + 2b, -J, -J, J - 2b give an extensive ground space only for J > 0, b = +-J (no two adjacent up, or down, spins; entropy log phi). Every other classical bond has ground degeneracy at most 2 on even rings and a gap. This confirms C141 stratum 1.

### 2.4 The free locus F

A bond is in F (after a rotation) iff J_zz = J_xz = J_zx = J_yz = J_zy = 0 and B is along z. With the Jordan-Wigner map along z,

    a = J_xy - J_yx,  c = J_xx + J_yy,  g = |(J_xx - J_yy) - i (J_xy + J_yx)|,  B = B_z,

the Bogoliubov quasiparticle energy is E(k) = 2 f(k) with

    f(k) = -a sin k + sqrt((c cos k + B)^2 + g^2 sin^2 k).                     (2.1)

[derived, checked] Derivation: s^+ s^- hopping amplitude t = c + i a, pairing amplitude J_xx - J_yy - i (J_xy + J_yx), so the even part of the normal dispersion is 2(c cos k + B), the odd part -2 a sin k, and the pairing 2 g |sin k|. [numerical, free_pt2.out] In the paramagnetic region (c = 1, g = 0.5, B = 1.5) the ED many-body gaps at L = 12, 14 agree with min over the allowed momentum grids of 2f to 1e-4 (0.75083 against 0.75060, 0.74288 against 0.74281), and at the critical a_c the L = 14 gap is 0.00444 against 0.00438.

F is three-dimensional in B_2 (parameters (a, c, g, B) up to scale). Its singular set is explicit.

Proposition 2.4. [proved, given (2.1)]
(i) Gapless region. For |a| <= g, f >= 0 with zeros only on the Ising surfaces |B| = |c| (at k = 0 or pi, linear dispersion with chiral slopes g -+ a, z = 1). For |a| > g, the chain is gapless with a Fermi sea of quasiparticles iff B^2 < a^2 - g^2 + c^2.
(ii) Flat band. f vanishes on an interval iff c = B = 0 and |a| = g. This is a single point of F (up to signs), the XZ point.
(iii) z > 1 points. Away from (ii), f has a degenerate zero exactly on the two Pokrovsky-Talapov surfaces
    PT1: a^2 = g^2 with B^2 < c^2  (touching at cos k0 = -B/c, f ~ c^2 sin k0 (k - k0)^2/(2g)),
    PT2: B^2 = a^2 - g^2 + c^2 with a^2 > g^2  (double root at cos k0 = -c B/(a^2 - g^2 + c^2)),
and at their intersections with the Ising surfaces. The order of the touching is at most two away from k0 in {0, pi}: f vanishes to the same order as the quadratic Q(u) = (a^2 - g^2 + c^2) u^2 + 2 c B u + B^2 - (a^2 - g^2) in u = cos k, which has at most a double root unless it vanishes identically (the XZ point), and u - u0 is linear in k - k0 when sin k0 != 0. On the curves |B| = |c|, |a| = g (where PT1, PT2 and the Ising surface meet) the touching is one-sided cubic, f ~ c^2 k^3/(8g) for a sin k > 0 and f ~ 2 g |k| on the other side [numerical check, free_z3_check.out: f(k)/k^3 = 0.2494, 0.2498, 0.2500 at k = 0.1, 0.05, 0.025 for c = 1, B = -1, a = g = 0.5]. There is no quartic onset anywhere in F.
(iv) Every singular point of F has a single quasiparticle band, and c <= 1 throughout F (the free-fermion theorem of the d = 2 program).

Proof. Gaplessness means f(k) <= 0 somewhere, i.e. (a^2 - g^2) sin^2 k >= (c cos k + B)^2 with a sin k > 0. If a^2 < g^2 the left side is negative unless sin k = 0, which needs c cos k + B = 0, i.e. |B| = |c|; near such k, (c cos k + B) = O(k^2) and f = (-a sgn + g)|k| + O(k^2). If a^2 > g^2, put u = cos k. The condition (a^2 - g^2)(1 - u^2) >= (c u + B)^2 is a quadratic inequality in u with discriminant proportional to (a^2 - g^2)(a^2 - g^2 + c^2 - B^2); it has solutions in (-1, 1) iff B^2 <= a^2 - g^2 + c^2 (boundary included). For (ii), f = 0 on an interval forces the identity (a^2 - g^2)(1 - u^2) = (c u + B)^2 in u, whose coefficients give c^2 = -(a^2 - g^2), c B = 0, B^2 = a^2 - g^2, hence c = B = 0 and a^2 = g^2. For (iii), at a^2 = g^2 and |B| < |c| we have f = g sin k [sqrt(1 + xi^2/(g^2 sin^2 k)) - 1] ~ xi^2/(2 g sin k) with xi = c cos k + B vanishing linearly at k0, which is a quadratic touching. On PT2 the quadratic inequality has a double root, so f has a double zero. QED.

[numerical, free_pt.out] Slice c = 1, g = 0.5: a_c = 0.5 = g for B = 0, 0.3, 0.6 (PT1, k0/pi = 0.5, 0.597, 0.705, f''(k0) = 2.00, 1.91, 1.60), and a_c = 1.2247 for B = 1.5 (PT2, k0/pi = 0.732, f'' = 1.37). All quadratic.

So inside F the singular set is two two-dimensional surfaces plus the XZ point. C141 stratum 3 lists only the U(1) points g = 0, B^2 = a^2 + c^2 (XX saturation, and its DM images). The PT surfaces with g > 0 are not U(1) symmetric. They bound the chiral gapless phase, and their gap minimum sits at an incommensurate k0, where pairing between k0 and k0 is forbidden by momentum and pairing between k0 and -k0 is off resonance because f(-k0) > 0.

### 2.5 Frustration-free bonds (Bravyi-Gosset)

[proved, cited, hypotheses checked] Bravyi and Gosset, J. Math. Phys. 56, 061902 (2015), arXiv:1503.04035. Every frustration-free nearest-neighbour qubit chain is, up to the null-space equivalence of their footnote 2, a sum of projectors of rank 1, 2 or 3.

- Rank 1, H = sum |psi><psi|. With T_psi = [[<psi|01>, <psi|11>], [-<psi|00>, -<psi|10>]] and eigenvalues mu_1, mu_2: the open-chain gap is at most 1/(n - 1) if |mu_1| = |mu_2| != 0, and bounded below otherwise (Theorem 1). Ground space: n + 1 dimensional for entangled psi (Prop. 2), the image of the symmetric subspace under 1 (x) T (x) T^2 (x) ...; for product psi = |1 v_perp> with v not parallel to |0>, n + 1 domain-wall states (Prop. 1); for symmetric product psi = |a a>, the Fibonacci space. On rings: n + 1 states if T^n is proportional to the identity, two otherwise (Prop. 3).
- Rank 2 (Appendix, Theorem 3): ground space of dimension 0, 1 or 2 for n >= 4; gapped except case 4 with |f| = 1 (span of a polarized state and a twisted one-flip state), where the gap is at most 1 - cos(pi/n).
- Rank 3: unique product ground state, gap 1.

[cited, recalled and not re-read in this run] Gosset and Mozgunov, J. Math. Phys. 57, 091901 (2016), arXiv:1512.00088: a translation-invariant frustration-free chain that is gapless in the thermodynamic limit has open-chain gap O(n^-2). So every gapless frustration-free qubit bond has z >= 2 in this sense.

[proved] Dimension count. T_psi is linear and bijective in psi, and psi -> (U (x) U) psi acts as T -> (det U)^-1 U T U^dagger. Modulo scale and conjugation, T has Schur form [[mu_1, x], [0, mu_2]] with x >= 0, three real parameters (mu_2/mu_1 and x/|mu_1|). The gapless condition |mu_2/mu_1| = 1 leaves two: the twist angle phi = arg(mu_2/mu_1) and the non-normality x/|mu|. The twisted FM Heisenberg bonds are x = 0, a one-parameter subfamily. So stratum 2 is two-dimensional in B_2, not one-dimensional as C141's description suggests. Bravyi and Gosset stress the same point ("the naive intuition is wrong"): psi = sqrt(1 - p)|00> + sqrt(p)|11> is gapless for every 0 < p < 1.

[derived] The real family psi_p = sqrt(1 - p)|00> + sqrt(p)|11> is, up to gauge and a constant, the XYZ bond with J = diag(s/2, -s/2, 1/4), s = sqrt(p(1 - p)), in the uniform field B = (0, 0, (1 - 2p)/2). (Pauli expansion: |00><00| and |11><11| give ZZ/4 and (1 - 2p)(Z1 + Z2)/4, the coherence gives s (XX - YY)/2.) This is the Kurmann-Thomas-Muller factorizing field h^2 = 4 (J_x + J_z)(J_y + J_z) for J_x = -J_y, and it conserves sum_j (-1)^j Z_j (section 2.7).

[numerical, ff_check.out] Rings, p = 0.2: exact ground degeneracy L + 1 on even rings and 2 on odd rings (T has eigenvalues +-i s, so T^L is proportional to 1 iff L is even), gap with L^2 gap = 10.9, 11.6, 11.9, 12.1 at L = 8, 10, 12, 14 (z = 2). Singlet (FM Heisenberg): L^2 gap = 18.7 to 19.4 (exact 1 - cos(2 pi/L) normalization). A non-normal example (twist pi/3, A = [[1, 0.7], [0, 1]]) has L + 1 states at L = 6, 12 and 2 otherwise, as Prop. 3 predicts.

Consequence for the list. Stratum 2 should read: the Bravyi-Gosset gapless set (rank one with |mu_1| = |mu_2| != 0, two-dimensional, z = 2; rank two case 4 with |f| = 1, z = 2), together with the gapped extensive point (product psi = |aa>, which is the Fibonacci point of stratum 1). The rank-two case 4 is the frustration-free form of saturation (a polarized state degenerate with a twisted one-flip band minimum), so it overlaps stratum 3 on its U(1) part.

### 2.6 The uniform U(1) class

h = J_perp (XX + YY) + D (XY - YX) + J_z ZZ + B (Z1 + Z2)/2. The DM term is a twist (C129): Delta_eff = J_z/sqrt(J_perp^2 + D^2). [proved in the literature, standard Bethe-ansatz completeness assumed] (Yang and Yang 1966; Griffiths 1964; Takahashi's book.) The phases are: polarized (gapped, z = 2 onset at the saturation field B_s = 2 (J_z + sqrt(J_perp^2 + D^2)) for Delta_eff > -1); Luttinger liquid with v > 0 at every magnetization strictly between the edges; for Delta_eff > 1 a gapped Neel plateau at B < B_c1 with a commensurate-incommensurate (z = 2, square-root magnetization) onset at B_c1 = 2 sqrt(J_perp^2 + D^2) sinh(gamma) sum_n (-1)^n/cosh(n gamma), cosh gamma = Delta_eff (des Cloizeaux-Gaudin gap); the ferromagnet Delta_eff <= -1 (gapped for B != 0, twofold at B = 0 for Delta_eff < -1, the FF multiplet at Delta_eff = -1, B = 0). Singular set: saturation surfaces, the h_c1 surfaces, the FM point, and the Fibonacci limit J_perp = D = 0, B = 2 J_z. [proved] No quartic onset occurs in this class: the one-magnon band over the polarized state is exactly 2 J_z-type constant plus 4 sqrt(J_perp^2 + D^2) cos(k - theta), quadratic at its minimum unless its amplitude vanishes (J_perp = D = 0, the flat band of the Fibonacci limit). The h_c1 surface is missing from C141's list. It is of type (Z): a single species (spinons, the Luttinger parameter tends to 1/4 in the convention where free fermions have 1 at the edge), v ~ M k_F.

[numerical, scan_uniform.jsonl] On the slice D = 0, (J_perp, J_z) = (cos t, sin t), B in [0, 3] (24 x 25 grid, L = 8, 10, 12), g exceeds 0.12 on 7 percent of the grid: a ridge just below B_s (0.1 to 0.2 below, the finite-size offset expected at L <= 12) and the Fibonacci peak (g = 0.955 at t = pi/2, B = 2). The h_c1 onset is not resolved at these sizes (its signal is weaker than the saturation one, since the onset is of spinon pairs on even rings).

### 2.7 The staggered U(1) class (new to the program)

Bonds that commute with s^z (x) 1 - 1 (x) s^z (up to gauge) conserve Q_st = sum_j (-1)^j Z_j. Up to rotation they are h = g (XX - YY) + J_z ZZ + B (Z1 + Z2)/2: a two-dimensional family in B_2. On even rings, a pi rotation about x on odd sites maps it to the FM-side XXZ chain g (XX + YY) - J_z ZZ in a staggered field B. This is not Bethe-integrable, and it is not in the scope of C129 to C140, whose U(1) class is the uniform one.

[proved] Its singular set contains the Bravyi-Gosset curve B = 2 sqrt(J_z^2 - g^2), J_z >= g > 0 (section 2.5, the psi_p family with g/J_z = 2s), which joins the FM Heisenberg point (g = J_z, B = 0) to the Fibonacci point (g = 0, B = 2 J_z). On it the ring ground space is exactly (L + 1)-fold on even rings and z = 2.
[derived, standard bosonization] Elsewhere it is regular. For J_z < g (Delta = -J_z/g > -1 in the rotated frame) the staggered field couples to the longitudinal Neel operator cos(2 phi), of scaling dimension K in the convention K = 1 for free fermions. It is relevant for K < 2 (gap, unique ground state) and irrelevant for K > 2, i.e. close to the FM point (Delta below about -0.71), where the Luttinger liquid survives with z = 1. Either way the bond is regular away from the FM point. For J_z > g the chain interpolates between the twofold ferromagnet (small B) and the pinned Neel state (large B), with the transition at the Bravyi-Gosset curve.
[numerical, scan_staggered.jsonl] On a 24 x 25 grid (t with (g, J_z) = (cos t, sin t), B in [0, 3]) the indicator g exceeds 0.12 on 5.7 percent of the grid, all of it within one grid spacing (0.125 in B) of the predicted curve or of the FM point. Peak values: B = 1.000, 1.375, 1.625, 1.875, 2.000, 2.000 against the predicted 0.948, 1.379, 1.664, 1.854, 1.964, 2.000 at t/pi = 0.286, 0.329, 0.372, 0.414, 0.457, 0.500. No other singular region.

### 2.8 XYZ in a uniform field along a principal axis

This family (J diagonal, B along z, three-dimensional in B_2) contains 2.3, the free locus J_z = 0, the uniform U(1) class J_x = J_y and the staggered one J_x = -J_y. [derived, from Bravyi-Gosset] Its Kurmann-Thomas-Muller surface B^2 = 4 (J_x + J_z)(J_y + J_z) (two product ground states, a disorder surface inside the ordered phase) is frustration free with a rank-two projector of Bravyi-Gosset case 2 (G = span{alpha^n, beta^n}), hence gapped, except on the curve J_x = -J_y where the projector drops to rank one with T_psi of equal-modulus eigenvalues and the chain becomes gapless with L + 1 ground states. [numerical, ktm_check.out, ktm_rank.out] With H = sum j_x XX + j_y YY + j_z ZZ + b sum Z and b_f = 2 sqrt((j_x + j_z)(j_y + j_z)): at (j_x, j_y, j_z) = (0.6, 0.2, 0.3), (0.4, -0.2, 0.3), (0.3, -0.25, 0.3), (0.5, 0.5, 0.3), (0.6, 0, 0.3) the gauge-optimized bond has a twofold minimum (projector rank 2), E_0/L equals the bond minimum to 1e-15, the ring ground space is twofold, and the gap is 0.555 (converged), 0.200 (exact) and, for (0.6, 0.2, 0.3), 0.283, 0.246, 0.223 at L = 12, 14, 16 (slowly converging, consistent with the large correlation length of case 2 when |<alpha|beta>| is close to 1). At (0.3, -0.3, 0.3), the zero-field (FM Heisenberg) end of the J_x = -J_y curve, the ground space is L + 1 dimensional and L^2 gap = 22.5, 22.9, 23.2 at L = 8, 10, 12. This answers the spec's question about "XYZ points where the factorized ground state becomes degenerate": they are exactly the staggered-U(1) Bravyi-Gosset curve of 2.7, z = 2, single species. The other singular points of the family are the band-onset (PT and saturation) surfaces and Ising-PT junctions, which `contrib/d2_multicritical/` found numerically to be "band edge or Lifshitz points with z > 1, or a classical multiphase point" (the Fibonacci point), consistent with this list.

### 2.9 The residual strata (one Z2, antiunitary, none)

No exact solution exists. [derived, standard theory of commensurate-incommensurate transitions: Pokrovsky and Talapov 1979; Japaridze and Nersesyan 1978; Schulz 1980] A gapped phase whose lowest excitation sits at an incommensurate momentum k0 (generic without inversion or time reversal, and also forced in chiral and spiral phases) has, near the gap closing, a dilute gas of quasiparticles whose number is conserved at low energy: creating two at k0 costs momentum 2 k0 != 0 mod 2 pi, and the pair (k0, -k0) costs e(k0) + e(-k0) > 0. The onset is a free-fermion z = 2 fixed point (the contact interaction psi^dagger d psi^dagger psi d psi has dimension 4 > 1 + z = 3 and is irrelevant), so the transition is a codimension-one hypersurface with z = 2 and one emergent species. These PT hypersurfaces exist in every stratum, bound floating and chiral c = 1 phases, and the free-locus surfaces PT1 and PT2 are examples. Inside the residual strata:

- Mechanisms for z > 1 are (a) band onsets (PT, codimension one, single species), (b) frustration-free ferromagnetic structure (2.5, codimension 6), (c) vanishing velocity at multicritical points (codimension >= 2; in F they are the one-sided junctions of 2.4(iii)).
- Mechanisms for exact extensive degeneracy are local IOMs (2.2: only Fibonacci and XZ in all of B_2) and frustration-free constraints (2.5: only Fibonacci).

[heuristic] No other mechanism for a nonconformal gapless point or an extensive degeneracy is known for nearest-neighbour qubit chains. The completeness of (a) to (c) at strong coupling is the remaining assumption, now in a sharper form: every z > 1 point outside F, FF and U(1) is a band onset of a single quasiparticle band or a multicritical junction of such onsets with z = 1 lines.

### 2.10 What is proved, and on which set

[proved] The amended list (section 0) is complete on: the classical bonds (2.3), the free locus F (2.4), the frustration-free bonds (2.5), and for every bond whose singularity is an exact extensive degeneracy produced by strictly local IOMs of range two, or of any range if rank J >= 2 (2.2). [proved in the literature, modulo Bethe-ansatz completeness] on the uniform U(1) class (2.6). [proved on a curve, derived elsewhere] on the staggered U(1) class (2.7).

An "open dense set" statement of the form "bonds outside the list are regular on an open dense set" is true but empty: the PT hypersurfaces show that the singular set has codimension exactly one, so its complement is open and dense whatever the list. The informative statement is the classification above. No new stratum outside the amended list was found analytically.

## 3. Task 3: numerical search

See section 3 of RESULT.md for the current numbers and NUMERICS notes below.

[numerical] Search design. Version 1 (single increment) and version 2 (random Nelder-Mead starts) are documented as failures in PROGRESS.md: the first is fooled by level crossings at one length, the second stalls on the g = 0 plateau of gapped phases (8 of 8 local v2 starts that finished ended at g between -0.001 and 0.035). Version 3 (screen_refine.py) screens random bonds uniformly on the 11-sphere of (J, B), keeps the best, and refines with Nelder-Mead.

[numerical] Screens. Locally 900 random bonds and on Typhon (job 17579868) 2880 random bonds, at L = 8, 10, 12 with kappa = 3. Quantiles of g (50, 90, 99 percent, max): local -0.041, 0.063, 0.205, 0.279; Typhon -0.036, 0.061, 0.210, 0.347. No random bond reaches the calibration values of the extensive points (Fibonacci 0.96, XZ 0.69) or of the frustration-free set (0.33 to 0.47). This is expected, since all singular strata have codimension at least one, and it shows that no open region of B_2 behaves singularly at these sizes.

[numerical] Refinement (local, L = 8, 10, 12, the 12 best screened bonds, maxfev 800, local_runs/v3_refined.jsonl). Final g between 0.26 and 0.55. None is classical (comm3 >= 0.73), none has a U(1) or staggered-U(1) symmetry (relative defects >= 0.94), none is frustration free (E_0/L above the gauge-optimized bond minimum by 0.1 to 0.6 in units of the bond norm), none is near the XZ orbit, and none has an exact ground degeneracy (deg = 1 at L = 8, 10, 12). All twelve carry a large DM vector (|D| between 0.70 and 1.10 for unit coefficient norm), i.e. they are chiral, which is where incommensurate band onsets are generic. The two strongest were analysed along the cut x(s) = (J, (1 + s) B):
- Point 3 (g = 0.55). The gapped side (s >= 0.05) has a single-quasiparticle band whose minimum sits at the incommensurate momentum k/pi = 1.29 to 1.33 (equivalently -0.67 to -0.71) and closes as s decreases: 0.107, 0.055, 0.025, 0.010 at s = 0.2, 0.15, 0.1, 0.05, identical at L = 12 and 14 (pt_cut_p3.out), roughly linear in s - s_c with s_c near 0.03, as for a Pokrovsky-Talapov onset (gap linear in the detuning, z nu = 1). On the other side (s <= 0) the low spectrum is a dense incommensurate liquid with many states below M/L (probe_v3.out). The classical product-state energy equals E_0/L to within 0.002 per site and has a unique uniform minimizer (classical_landscape_v3.out, uniform_minima_p3.out), so this is a weakly quantum chain next to a PT surface of the amended list (stratum 3).
- Point 4 (g = 0.36). Along the cut a k = pi level falls to 2e-4 (a twofold ordered ground state appears, regular) while a low mode at small momentum softens (pt_cut_p4_p0.out).
- Continuing the maximization (follow.py, kappa = 1, 2500 evaluations). Point 3 moves by 0.12 and stays a chiral generic bond with g = 0.57 and n(<M/L) = 16, 29, 48 at L = 8, 10, 12 (growth ratio falling, 1.81 then 1.66). Point 4 moves by 0.28 toward a nearly rank-one J with nearly orthogonal vectors (singular values 0.905, 0.072, 0.065; |a.b| = 0.083; |B|/s_1 = 0.46), with n(<M/L) = 8, 16, 32, a doubling every two sites, i.e. the 2^(L/2) structure of the XZ manifold. Minimizing the width of the lowest 2^(L/2 - 1) states at L = 10 from there (flatten.py, flatten_p4.out) drives the bond onto the XZ orbit: width 6.0e-2 -> 3.5e-4 -> 3.2e-4, J singular values -> 0.9975, 0.040, 0.034, |a.b| -> 5e-4, |B|/s_1 -> 0.048. So the strongest near-extensive signal found by the unbiased search is the basin of the XZ point, and no other zero-width point appeared.

[numerical] The Typhon refinement at L = 10, 12, 14 (48 starts) was still running when this file was written. Its results are in /scratch/lustre/jkudlerflam/c143/v3_refined.jsonl. See RESULT.md for the state at the end of the run.

[numerical, family scans] The exactly characterized families confirm the indicator: section 2.6 (uniform U(1)) and 2.7 (staggered U(1)).

Assessment. [numerical] The search did not find a singular point outside the amended list. It did rediscover the XZ basin from a random start. Its resolving power is limited: z = 2 onsets give g of order 0.1 to 0.5 at L <= 12, comparable to finite-size noise in incommensurate liquids, so the search can locate onset surfaces only roughly and cannot by itself certify z. Neither isolated point was hit directly by the 3780 random bonds, as expected for points of codimension 8, but the XZ basin was reached by refinement (point 4).

## 4. Task 4: local analysis at the new strata and checks of C141's claims

### 4.1 PT hypersurfaces and h_c1

[derived] C141 5.4 applies unchanged. Near a PT point h*, a perturbation delta opens a Fermi sea of size k_F ~ (delta/M)^(1/2) of one quasiparticle species, with c = 1 (a floating or chiral Luttinger liquid) on one side and a gap on the other, and v = M k_F ~ M (delta/M)^(1/2) -> 0. One species, so N_em = 1, and no fast sector with long crossover can come from it (v xi ~ M forces xi ~ (M/v)). The same holds at h_c1 (spinons) and at the one-sided junctions (one chirality).

[heuristic] Coexistence of a PT onset with an independent conformal sector would be a c > 1 point (a second species on top of a c >= 1/2 phase). In F this does not happen (one Bogoliubov band). Outside F it is the same question as the c > 1 searches of the d = 2 program (REPORT 9.37, `contrib/d2_multicritical/`), which found none.

### 4.2 The XZ point

Exact structure. [proved, xz_algebra.out] Put a_j = X_j Z_{j+1} and b_j = Z_j X_{j+1}. Then a_j^2 = b_j^2 = 1, neighbouring a's anticommute and distant ones commute (likewise for b), [a_i, b_j] = 0 for all i, j, and

    Y_j Y_{j+1} = a_j b_j.

With a Jordan-Wigner map along y (Y_j = -i g_{2j-1} g_{2j}), a_j = i g_{2j-1} g_{2j+1} and b_j = -i g_{2j} g_{2j+2}: the odd Majoranas form one chain and the even Majoranas another. At XZ, H = sum_j a_j is a uniform Majorana chain in the odd Majoranas (critical, c = 1/2) and the even Majoranas are free zero modes.

[numerical, xz_check.out] Rings L = 4 to 12: ground degeneracy 2^(L/2 - 1) exactly (2, 4, 8, 16, 32); L times the first gap = 3.18, 3.17, 3.16 at L = 8, 10, 12 (z = 1); the local IOMs Q_j = Z_j X_{j+1} commute with H to machine precision. E_0/L = -0.6533 at L = 8 (tending to -2/pi). Under an on-site Hadamard on alternate sites the uniform XZ chain becomes the isotropic quantum compass chain (sum X_{2i} X_{2i+1} + Z_{2i+1} Z_{2i+2}), for which macroscopic degeneracy and a gapless spectrum at the isotropic point are known (Brzezicki, Dziarmaga and Oles, PRB 75, 134415 (2007); Brzezicki and Oles, arXiv:0805.3904, which calls it the XZ model).

Isolation. [proved] By Theorem 2.2 and Proposition 2.4(ii), XZ is the only non-classical bond with a range-two IOM and the only flat-band point of F. [numerical, xz_dirs.out] Of the 11 tangent directions, only X(x)Y and Y(x)Z keep the 2^(L/2-1) degeneracy exactly at finite epsilon; both are infinitesimal on-site rotations of XZ. All others split it: Y(x)Y and Z(x)X at first order, the remaining seven at second order (widths ~ epsilon^2 for the quadratic directions XX, ZZ, B_y; apparent exponent 1.6 to 2 for YX, ZY, B_x, B_z, where the splitting competes with the chain gap 3.2/L). So XZ is an isolated point of B_2 up to images.

The anisotropic Ashkin-Teller neighbourhood. [proved] On the two-parameter family h = XZ + alpha ZX + beta YY,

    H = sum_j [ a_j + alpha b_j + beta a_j b_j ],

two commuting Majorana chains with hopping ratio alpha and an energy-energy coupling beta. For alpha = 1 this is the self-dual quantum Ashkin-Teller chain in Majorana form (Kohmoto, den Nijs and Kadanoff 1981). Equivalently, a 45-degree rotation about y turns XZ + ZX into X'X' - Z'Z', and a pi rotation about x' on odd sites (even rings) gives the XXZ chain about the y axis with Delta = -beta: c = 1 for |beta| < 1, with the FM Heisenberg point at beta = 1. For alpha != 1 it is the same model with unequal velocities. XZ is the point alpha = beta = 0 where the second chain has no dynamics.

First-order emergent model. [proved at first order] At XZ on a finite ring the A chain has a unique ground state in the relevant sector and the B algebra acts on the degenerate manifold. Projecting, P (a_j b_j) P = <a_j>_A b_j, so the first-order Hamiltonian on the manifold is

    H_eff^(1) = (alpha + beta <a>) sum_j b_j,      <a> = E_0(L)/L -> -2/pi.

[numerical check] The first-order splitting widths at L = 8 are 10.45 epsilon (ZX) and 6.83 epsilon (YY); their ratio 0.653 equals |E_0/L| = 0.6533 at L = 8.

So for generic directions the emergent sector is a uniform Majorana chain with velocity v_B proportional to |alpha + beta <a>| = O(delta M): c_em = 1/2, total c = 1, A = (v_A + v_B)/2 + O(delta). This is C141's (5.1) with the Schrieffer-Wolff gap replaced by the exact commutation of the two algebras. [derived] The energy-energy coupling eps_A eps_B has dimension 2 and is marginal; the Ising OPE eps x eps = 1 + T contains no eps, so it generates no beta^2 term in its own flow, and the two Majorana sectors survive with c = 1 for small beta, as for the isotropic Ashkin-Teller line.

[numerical, xz_at_sym.out, symmetry-resolved ground states, L = 8 to 16] c estimates from S(L/2) at L and L + 4: (alpha, beta) = (1, 0): 0.975, 0.983, 0.988; (1, 0.5): 0.980, 0.986, 0.990; (0.3, 0): 0.975, 0.983, 0.988; (0, 0.3): 0.996, 1.025, 1.011; (0, -0.3): 0.996, 1.006, 1.011. So c = 1 near XZ off the tuned cone, as predicted.

The tuned cone. On alpha + beta <a> = 0 (alpha = (2/pi) beta to leading order) the first-order hopping vanishes. The B chain then acquires dynamics only at second order, through the fluctuations of a_j: H_eff^(2) contains next-nearest hopping (from b_j b_{j+1} = -i g g') and four-Majorana terms b_j b_l (|j - l| >= 2) with couplings set by the connected eps_A correlator, which decays as a power. This is an interacting Majorana chain of the Rahmani-Zhu-Franz-Affleck type (PRL 115, 166401 (2015); PRB 92, 235123 (2015)) with long-range tails and with a ratio of four-fermion to two-fermion couplings of order one. Such chains have, besides Ising, a tricritical Ising point (c = 7/10) and an incommensurate c = 3/2 phase. [open] Whether the tuned cone reaches one of them, which would give c_tot = 6/5 or 2 at a nearest-neighbour qubit bond, is not settled. [numerical, xz_at_sym.out, symmetry-resolved, unique ground states except where noted] On the tuned cone the half-ring entropies are large and oscillate with the parity of L/2:
(alpha, beta) = (0.064, 0.1): S = 1.26, 1.84, 1.79, 1.94, 1.94 at L = 8, 10, 12, 14, 16;
(0.19, 0.3): 1.30, 1.79, 1.64, 1.92, 1.77;
(-0.19, -0.3): 1.30, 1.42, 1.64, 1.68, 1.77;
(0.38, 0.6): 1.12, 1.30, 1.28, 1.40, 1.40 (twofold at L = 10, 14);
while far out on the cone, (0.64, 1.0) and (1.27, 2.0), S saturates at 0.79 and 0.71 (gapped; xz_at3_tuned.out), and on the axis beta = -1 it falls (c_est 0.79 to 0.62), beta = +-2 is gapped. The excess over the c = 1 values (about 1.2 at these sizes) is close to log 2 at the smallest couplings, which is what a B chain still inside its nearly degenerate regime (crossover length beyond L = 16, since its hopping is O(beta^2)) would give. The c estimates from L and L + 4 decrease with L (2.5, 1.1 to 1.5 at (0.19, 0.3); 2.5, 2.3, 1.4 at (-0.19, -0.3)). These data are not converged. They neither establish nor exclude c_B > 1/2 on the cone.

Consequences for C141. [derived] (i) v_em = O(delta M) still holds at XZ (the emergent Hamiltonian on the B algebra is O(delta) at first order and O(delta^2) on the tuned cone), so no fast emergent sector arises and local boundedness of A at XZ reduces, as in C141 5.3, to boundedness of c for the emergent Majorana chains. (ii) The non-coexistence statement and the bound c_tot <= 1 at d = 2 are not established at XZ: c_tot = 1/2 + c_B, with c_B = 1/2 generically and unknown on the tuned cone.

### 4.3 Checks of C141's single-species claims at the other strata

- Cluster argument (no second fluid of bound pairs) [proved in U(1) by Bethe ansatz, derived otherwise]. In the uniform U(1) class the ground state at any magnetization is a single sea of real rapidities for Delta_eff > -1 (Yang-Yang); string (bound-state) solutions never enter the ground state, and for Delta_eff < -1 the attraction collapses the system into the ferromagnet. This is the content of the (n - 1) V argument: in a one-species nearest-neighbour lattice gas, an attraction that binds pairs binds n-clusters more strongly per particle ((n - 1)V/n increases with n), so the gas phase separates instead of forming a pair fluid. The argument is correct for nearest-neighbour attraction; it does not exclude pair fluids stabilized by longer-range repulsion, which at first order near a qubit singular point is absent.
- p-wave pairing under U(1) breaking [derived]. At a z = 2 onset of spinless fermions near k = 0 or pi, the pairing psi d_x psi has scaling dimension 2 < 3 = 1 + z and is relevant; the chemical potential (dimension 1) is more relevant. The flow ends in a gapped phase or an Ising transition with v proportional to the pairing amplitude (free check: XY in a field has v_Ising = 2 gamma). At an incommensurate onset (PT) pairing is momentum-forbidden and the transition stays z = 2. Both agree with Proposition 2.4.
- PXP at the Fibonacci point [proved for the exactly solvable cases, literature]. Perturbing n (x) n by U(1)-preserving hopping gives free fermions on L - N sites (C141), c = 1. Perturbing by a transverse field gives the PXP / hard-boson model of Fendley, Sengupta and Sachdev (PRB 69, 075106 (2004)), whose phases are disordered, Z2-ordered and the Ising (c = 1/2) and tricritical Ising (c = 7/10) transitions; with second-order next-nearest repulsion, Z3 (Potts c = 4/5, chiral) and floating (c = 1) phases appear. All have c <= 1 and v_em = O(delta) (fib_emergence.py of C141). The Fibonacci stratum is of type (G) with an empty conformal part, so C141's analysis stands there.
- Frustration-free stratum [derived]. The ground manifold is the image of a symmetric multiplet under a product map, the low excitations are one-flip ("magnon") states, one band. Under perturbation: a field-like term selects a product state (gapped), an easy-plane term gives a Luttinger liquid with v ~ M (delta/M)^(1/2) (C136, C138 in the U(1) case). Single species.

## 5. Conditional derivation of C(2) < infinity

[derived, conditional] Assumptions, each tagged:

- (IR) [assumed] Each fixed bond has finitely many unitary CFT sectors or is gapped (C141 0).
- (LR) [proved, C109] v_a <= 4eM.
- (WCT) [assumed] Weak-perturbation c-theorem at regular bonds (C141 5.2): a weakly perturbed multi-velocity CFT does not increase c_tot.
- (COMP) [proved on classical, F, FF, local-IOM degeneracies; literature-proved on uniform U(1); proved on a curve and derived elsewhere on staggered U(1); heuristic in the residual strata] Every singular bond of B_2 lies in the amended list: Fibonacci points; the Bravyi-Gosset gapless set; band-onset (PT, saturation, h_c1) hypersurfaces and their multicritical junctions; the XZ point.
- (EM-G) [derived, C141 5.3] Emergence at gapped degenerate manifolds (Fibonacci): v_em <= C delta M.
- (EM-Z) [derived, C141 5.4 and section 4.1] Emergence at z > 1 onsets: one species per band minimum, v_em ~ M (delta/M)^(1/z).
- (EM-XZ) [derived at first order, open on the tuned cone] At XZ the emergent B chain has v_em = O(delta M) (first order) or O(delta^2 M) (tuned cone), and its central charge c_B is bounded, uniformly in the direction of approach.
- (NEST) [heuristic] Finite nesting depth (C141 6.2), now including the XZ recursion: the emergent model at XZ is an interacting Majorana chain, whose own singular points would be nested strata.

Derivation. B_2 is compact (C141 5.1). Under (IR), (WCT) and (LR), A is locally bounded at every regular bond. Under (COMP) every remaining bond is on the amended list. At Fibonacci points (EM-G) and at FF and onset points (EM-Z), A(h* + delta) <= A_cont + N_em C |delta|^alpha M with N_em finite and A_cont = 0 (empty conformal part), so A is locally bounded. At XZ, A(h* + delta) <= (1/2) V + c_B v_B with v_B = O(delta M), bounded under (EM-XZ). With (NEST), Heine-Borel gives sup_h A(h)/M < infinity, i.e. C(2)_asymptotic < infinity. The finite-size C(2) follows from it under C131's uniform remainder lemmas (C131 PROOF 3), which are not re-proved here.

The sharper statement of C141, c_tot <= 1 and A <= 4eM at d = 2, needs in addition (NC) [open at XZ]: c_B <= 1/2 on the tuned cone of the XZ point, and absence of c > 1 at regular bonds [numerical, REPORT 9.37].

## 6. Consistency checks

| Check | Expected | Found |
| :-- | :-- | :-- |
| ED vs dense diagonalization (XX, L = 8) | equal | 2e-14 (ed.py) |
| Heisenberg E_0/4, L = 8 | -3.6510934 | -3.651093408937 |
| Critical TFIM E_0, L = 10 | free-fermion sum | agree to 7e-15 |
| Free dispersion (2.1) vs ED, gapped region | equal | 2e-4 at L = 14 |
| XZ degeneracy | 2^(L/2-1) by Majorana count | 2, 4, 8, 16, 32 at L = 4 to 12 |
| XZ: first-order YY/ZX splitting ratio | |<a>| = |E_0/L| | 0.653 vs 0.6533 (L = 8) |
| Bravyi-Gosset psi_p, even rings | L + 1 ground states, z = 2 | L + 1, L^2 gap -> 12.1 |
| Bravyi-Gosset curve in staggered U(1) | B = 2 sqrt(J_z^2 - g^2) | scan peaks within one grid step |
| Fibonacci limit of the BG curve | g -> 0, B -> 2 J_z | peak g = 0.955 at (t = pi/2, B = 2) |
| Ashkin-Teller line alpha = 1 | c = 1 | 0.98 to 0.99 |
| Two decoupled Majorana chains (alpha = 0.3, beta = 0) | c = 1 | 0.975 to 0.988 |
