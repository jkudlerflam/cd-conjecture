# C154 derivations: qubit bonds without symmetry

This run was done by a Claude Opus agent (2026-11-07), not by Codex.

Standard: theoretical physics. Tags: [proved] (complete argument, or an exact identity checked to machine precision), [derived] (physics-standard derivation with named inputs), [numerical] (sizes, scripts and outputs named), [heuristic], [assumed]. The tag applies to the paragraph or displayed claim it opens.

Conventions. H_L = sum_j [ sigma_j^T J sigma_{j+1} + B . sigma_j ] on a ring of L sites, Pauli matrices, J real 3 x 3, B real. B_2 is (J, B) modulo a common rotation and positive scale, dimension 8. "Symmetry-free" means no on-site symmetry of any kind (uniform or site-dependent, unitary or antiunitary); translation is always exact. Strata T, K, N are those of C68 (T: only Kramers time reversal, B = 0; K: only one square-positive antiunitary reflection; N: no uniform internal symmetry). Prior numerical studies of the same questions are contrib/d2_floating_ising, d2_multicritical, d2_two_fluid and d2_classification, and the analyses C63, C67, C68; this run builds on them and does not repeat them.

## 1. Site-dependent symmetries and the complete list of symmetric strata

C68 classified uniform internal symmetries and stated that "spatial and site-dependent symmetries must be checked separately". This section closes the site-dependent case.

### 1.1 Setup

[proved] Let U = (tensor_j u_j) or U = (tensor_j u_j) K be a product of on-site unitaries, possibly followed by complex conjugation, with U H_L U^{-1} = H_L on a ring of L >= 3 sites. It acts on Pauli vectors by U sigma_j U^{-1} = Q_j sigma_j with Q_j in SO(3) for unitary U and Q_j in O(3) \ SO(3) for antiunitary U (complex conjugation acts as diag(1, -1, 1)). Since the Pauli expansion of H_L on L >= 3 sites is unique (the coefficient of sigma^a_j sigma^b_{j+1} is J_ab and that of sigma^a_j is B_a), invariance is equivalent to

    Q_j^T J Q_{j+1} = J   and   Q_j^T B = B   for all j (mod L).      (1.1)

### 1.2 Polar reduction

Write J = S O with O in SO(3) and S symmetric. (Take the polar decomposition J = (J J^T)^{1/2} O'; if det O' = -1 use S = -(J J^T)^{1/2}, O = -O'.) Let s_1 >= s_2 >= s_3 >= 0 be the singular values of J (the |eigenvalues| of S), with eigenvectors e_1, e_2, e_3 of S, and let G_S = {diag(+-1, +-1, +-1) in the frame (e_1, e_2, e_3)}, the group of orthogonal matrices commuting with S when the s_i are distinct.

Lemma 1. [proved] If s_1 > s_2 > s_3 (rank 2 or 3), every solution of (1.1) satisfies Q_j in G_S and Q_{j+1} = O^{-1} Q_j O for all j.

Proof. Put W_j = O Q_{j+1} O^{-1}, of the same determinant as Q_{j+1}. Then (1.1) reads S W_j = Q_j S. Multiplying by its transpose on either side gives S^2 = Q_j S^2 Q_j^T and W_j^T S^2 W_j = S^2, so Q_j and W_j commute with S^2 = J J^T. Its eigenvalues s_i^2 are distinct, so Q_j and W_j are diagonal (signs) in the frame e_i, hence commute with S, and S W_j = Q_j S gives W_j = Q_j on the range of S. If s_3 > 0 this is everything. If s_3 = 0 the two diagonal sign matrices agree on e_1, e_2, and the third sign is fixed by the common determinant (all Q_j have the same determinant, since U is either unitary or antiunitary as a whole). So W_j = Q_j, i.e. Q_{j+1} = O^{-1} Q_j O. QED.

### 1.3 Classification for distinct singular values

Theorem 1. [proved] Let s_1 > s_2 > s_3. The on-site symmetries of H_L (L >= 3, any site dependence) are Q_j = O^{-j} Q_0 O^{j} with Q_0 in the group A = {Q in G_S : O^{-n} Q O^{n} in G_S for all n}, subject to Q_j^T B = B for all j (and O^{-L} Q_0 O^{L} = Q_0 on the ring). Explicitly:

(i) Q_0 = 1 always; Q_0 = -1 (Kramers time reversal Theta = prod_j (iY_j) K, Theta^2 = -1 on a site) iff B = 0. These are uniform.

(ii) Q_0 = R_pi(e_i) (pi rotation about a principal axis) or its antiunitary partner -R_pi(e_i) (reflection through the plane normal to e_i) lies in A iff the orbit O^{-n} e_i stays among {+-e_1, +-e_2, +-e_3}. Two cases occur.
  (a) O e_i = +-e_i. Then O^{-1} R_pi(e_i) O = R_pi(e_i) and the symmetry is uniform. The unitary one needs B parallel to e_i: this is the parity class of C67 and C152. The antiunitary one needs B perpendicular to e_i: this is C68's stratum K.
  (b) O e_i = +-e_k with k != i. If O^{-1} maps two orthogonal principal axes to principal axes it maps the third one too, so O belongs to the octahedral rotation group O_S of the principal frame (signed permutations with determinant one). Conversely, for O in O_S every element of G_S is in A.

(iii) For O in O_S the possibilities are: O = 1 or a pi rotation about a principal axis (J symmetric in the frame, uniform Klein group: XYZ-type, plus T if B = 0); O a rotation by pi/2 about e_k or by pi about a face diagonal (e_i +- e_j)/sqrt 2 (then R_pi(e_k) is uniform, so the bond is in the parity class, and R_pi(e_i), R_pi(e_j) alternate with period 2 as an extra staggered symmetry when B is parallel to e_k or zero); O a rotation by 2 pi/3 about a body diagonal (+-e_1 +- e_2 +- e_3)/sqrt 3 (the three axes are permuted cyclically and the whole Klein group is twisted with period 3; no nontrivial uniform unitary survives; the field condition forces B = 0, since R_pi(e_1), R_pi(e_2), R_pi(e_3) all occur along the chain and only B = 0 is fixed by all of them, and likewise for the reflections).

Proof. Lemma 1 gives the form of Q_j. The orbit condition in (ii) is the statement that every Q_j lies in G_S, whose nontrivial elements other than -1 are pi rotations about principal axes and their negatives. The rest is the enumeration of the 24 elements of O_S by conjugacy class. QED.

Corollary 1 (the cyclic stratum). [proved] The only site-dependent symmetry not accompanied by a uniform unitary one (for distinct singular values) is the period-3 twisted Klein group, with B = 0. Up to a common rotation and sign flips of principal axes these bonds are exactly

    h = alpha XZ + beta YX + gamma ZY,     J = D P,  D = diag(alpha, beta, gamma),      (1.2)

with P the cyclic permutation P e_z = e_x, P e_x = e_y, P e_y = e_z (a 2 pi/3 rotation about (1,1,1)). This is a 2-dimensional stratum of B_2. It contains the XZ point (1,0,0) and the YX and ZY directions at XZ studied in C144.

Degenerate singular values. [derived, sketched] If two singular values coincide the commutant of J J^T contains O(2), continuous symmetries can occur, and these are exactly the uniform and staggered U(1) charges classified in contrib/d2_gap_reduction/additive_charge.md (XXZ + DM + field, and a(XX - YY) + uZZ + bZ). The discrete part follows the same orbit argument with G_S replaced by the commutant; for generic O (not a rotation about the distinguished axis) only +-1 survive, so a degenerate pair of singular values alone does not create a symmetry. Rank-one J (s_2 = s_3 = 0) has the extensive local symmetries of C143 (classical Ising, XZ), already in the singular list.

### 1.4 The strata of B_2

[proved for the counts, using Theorem 1, C67, C68 and additive_charge.md]

| stratum | symmetry | LSM-type anomaly | dim in B_2 | status of c <= 1 |
| :-- | :-- | :-- | :-- | :-- |
| N | none | none | 8 | open (this run, C155) |
| T | Kramers Theta only, B = 0 | yes, Theta^2 = -1 per cell | 5 | open |
| K | one antiunitary reflection | none (square +1) | 5 | open |
| parity | uniform pi rotation | none | 4 | C152 |
| cyclic (new) | twisted Klein (period 3) + Theta | yes (Klein projective on a 3-site cell, and Theta) | 2 | this run, section 3 |
| U(1) | XXZ + DM + field | at half filling only | 3 | established |
| staggered U(1) | a(XX - YY) + uZZ + bZ | no | 2 | C153 |
| XYZ, B = 0 | uniform Klein + Theta | yes | 2 | Baxter, c = 1 lines |
| free locus | quadratic after Jordan-Wigner | | 3 | proved c <= 1 (C143) |
| singular strata | Fibonacci, FF, XZ, onsets | | various | C141 to C144, this run section 3 for YX/ZY |

Spatial reflections composed with on-site rotations (J^T = R^T J R, R^T B = B) give further strata of dimension at most 5 (for R = 1 they are the bonds with no DM vector). They matter in section 5 because a reflection can protect a chirality Ising order.

### 1.5 What survives without symmetry

[proved, standard] On N the only exact symmetry is translation. Translation carries no 't Hooft anomaly of its own (its infrared image is an internal symmetry generated by one element, and H^3(Z, U(1)) = 0), so no LSM-type statement constrains gaplessness, and any (CFT, automorphism) pair is kinematically allowed as an infrared fixed point. What survives is a selection rule: the infrared image g_T of translation must leave every allowed perturbation invariant, i.e. lattice momentum is conserved. Two consequences are used below. First, when g_T acts through a torus (incommensurate momenta), every operator that is charged under the closure of {g_T^n} is forbidden, which is how floating phases become open regions (C68 section 4). Second, when g_T acts as a Z_p (period-p order) it forbids the order field itself. There is no filling constraint (Oshikawa's argument needs a U(1)).

[proved, literature] On T and on the cyclic stratum, Ogata and Tasaki (Commun. Math. Phys. 372, 951 (2019)) forbid a unique gapped ground state for a half-odd-integer-spin chain with time reversal, or with the Z2 x Z2 pi rotations, and translation. For the cyclic stratum use translation by three sites, which commutes with the twisted Klein group (O = P has order 3) and has spin 3/2 per cell.

## 2. The spiral frame: every bond is a U(1)-class bond plus momentum-carrying harmonics

### 2.1 Exact identity

[proved] Let J = S O as in 1.2, and let O be a rotation by theta in [0, pi] about the unit axis m. The site-dependent rotation sigma_j = O^{-j} tau_j (implemented by on-site unitaries; on a ring it leaves a twist O^L at one bond) maps H_L to

    H_L = sum_j [ tau_j^T (O^j S O^{-j}) tau_{j+1} + (O^j B) . tau_j ].      (2.1)

Decompose S and B by angular momentum about m: S = S^(0) + S^(+-1) + S^(+-2), B = B^(0) + B^(+-1), with S^(0) = s_m m m^T + s_perp (1 - m m^T), s_m = m.S.m, s_perp = (tr S - s_m)/2, B^(0) = (m.B) m. A component of angular momentum q about m picks up the factor e^{i q theta j} in (2.1). So, in the tau frame,

    H_L = H_par + sum_j sum_{q = +-1, +-2} e^{i q theta j} V_q(j, j+1),      (2.2)

where H_par is translation invariant and U(1)_m symmetric: XXZ about m with transverse coupling s_perp, longitudinal coupling s_m, and field m.B. Back in the sigma frame the parent bond is J_par = S^(0) O (XXZ about m plus a DM vector along m, with J_perp = s_perp cos theta, D = s_perp sin theta) and B_par = (m.B) m. This is the uniform U(1) class.

[proved] The straight path J(lam) = (S^(0) + lam (S - S^(0))) O, B(lam) = B^(0) + lam (B - B^(0)), lam in [0, 1], stays in B_2, starts in the U(1) class and ends at the given bond. (spiral.py checks J(1) = J to 6e-16 and the exact U(1)_m symmetry of J(0), B(0) to 2e-16.)

### 2.2 Momentum of the harmonics in the parent liquid

[derived; inputs: bosonization of the XXZ chain in a field, as in C68 section 4] Suppose the parent is in its Luttinger-liquid regime, with Luttinger parameter K, Fermi momentum k_F (in the tau frame, set by the magnetization along m) and bosons (phi, vartheta). A charge-q harmonic of (2.2) couples to vertex operators e^{i(q vartheta + 2 n phi)}, whose tau-frame momenta are q pi + 2 n k_F (mod 2 pi) in the standard expansion, while the harmonic itself carries the spatial phase e^{i q theta j}. The coupling is translation invariant, hence able to grow under RG, only on the resonance set

    q (theta + pi) + 2 n k_F = 0  (mod 2 pi),   q in {+-1, +-2}, n in Z, dimension below 2.      (2.3)

Each resonance is one equation on the parameters (k_F moves with the field and the couplings), so the resonant set is a union of hypersurfaces. Off resonance every charged harmonic is irrelevant at weak amplitude, second order produces only charge-zero terms (which renormalize the parent) and charge +-2, +-4 terms at momenta 2 theta, 4 theta (still off resonance), and the infrared is the parent's c = 1 boson. Translation then acts on (phi, vartheta) by two independent incommensurate shifts, so the emergent symmetry is U(1) x U(1) and every vertex operator is forbidden: this is C68's floating mechanism, now available around every gapless parent rather than one anchor.

### 2.3 Consequence: a canonical path from every bond to the U(1) class

[derived] For every bond whose parent is a Luttinger liquid off resonance there is lam_0 > 0 such that every point of the path with lam < lam_0 is a c = 1 floating liquid connected to the parent through its own critical region, with an emergent U(1) x U(1) at the fixed point. Near a resonance the harmonic locks a commensurate gapped phase (period set by (2.3)), bounded by commensurate-incommensurate transitions (Pokrovsky-Talapov, z = 2, or Kosterlitz-Thouless, c = 1). Every critical point reached in this regime is the parent's c = 1 fixed point perturbed by a weak relevant coupling, so c <= 1 by the c-theorem. If the parent is gapped (Neel along m, or polarized) a weak harmonic leaves it gapped. Parents at singular points (ferromagnet, saturation, Delta = 1 end) are the singular strata of C141 and C143.

[proved, limitation] lam_0 is not uniform, and generic bonds sit at lam = 1 with harmonics that are not small: for C68's symmetry-free floating point 1282 the harmonic fraction |S - S^(0)|/|S| is 0.61 and |B - B^(0)|/|B| is 0.57 (spiral.py). The parent of 1282 is XXZ with Delta = -0.414 and h/J_perp = 0.22, a Luttinger liquid. So section 2.3 is a perturbative statement around the U(1) class, and the path is only a candidate route at lam = 1. Section 6 tests it numerically on 1282 and on the cyclic stratum.

## 3. The cyclic stratum and the YX and ZY directions at XZ

### 3.1 Exact map to a period-3 XYZ chain

Proposition 2. [proved] For h = alpha XZ + beta YX + gamma ZY (J = D P as in (1.2)) the site-dependent rotation sigma_j = P^{-j} tau_j maps the ring of length L = 0 mod 3 to the XYZ ring

    H = sum_j sum_{a = x,y,z} J^a_j tau^a_j tau^a_{j+1},   (J^x_j, J^y_j, J^z_j) = (alpha, beta, gamma), (gamma, alpha, beta), (beta, gamma, alpha) for j = 0, 1, 2 mod 3.      (3.1)

Proof. Bond j becomes tau_j^T P^j D P P^{-j-1} tau_{j+1} = tau_j^T (P^j D P^{-j}) tau_{j+1}, and P^j D P^{-j} is diagonal with cyclically permuted entries. For L not divisible by 3 the tau ring carries the twist P at one bond. QED. Check: cyc_check.py compares the full spectra of the two rings at L = 6, 9 for four parameter points (agreement 3e-14), and verifies that the twisted Klein generators U_n = prod_j (P^{-j} e_n) . sigma_j and Kramers Theta commute with H exactly.

So the YX direction at XZ, (1, eps, 0), is the period-3 chain with bonds XX + eps YY, YY + eps ZZ, eps XX + ZZ, and the ZY direction (1, 0, eps) is its mirror image (bonds XX + eps ZZ, eps XX + YY, eps YY + ZZ). At eps = 0 (3.1) is the three-colour Kitaev chain XX, YY, ZZ, which is the XZ point in this gauge. The period-3 commensuration that C144 saw on rings is the twist P on rings with L not divisible by 3.

### 3.2 Equivalences, the fundamental domain, and special points

[proved] A uniform pi rotation about a principal axis flips the signs of two of (alpha, beta, gamma); a uniform rotation by P permutes them cyclically; spatial reflection composed with a rotation exchanging x and z exchanges beta and gamma (up to the sign flips). Hence the stratum is the 2-sphere of (alpha, beta, gamma) modulo signed permutations with an even number of sign changes, and a fundamental domain is the spherical triangle alpha >= beta >= |gamma|. Its vertices are XZ (1, 0, 0), the antiferromagnetic Heisenberg point (1, 1, 1) (bonds tau.tau, SU(2) exact in the tau frame) and the ferromagnetic Heisenberg point (1, 1, -1). Its edges are beta = gamma (XZ to AF), beta = -gamma (XZ to FM) and alpha = beta (AF to FM through XZ + YX = (1, 1, 0)). The YX line (1, eps, 0) runs from XZ through the interior to the midpoint (1, 1, 0) of the AF-FM edge.

### 3.3 Symmetry, anomaly, and emergent SU(2) in the Heisenberg basin

[proved] The internal image of the symmetries is generated by the Klein group {1, R_pi(x), R_pi(y), R_pi(z)} acting on tau, the cyclic permutation P of the tau axes (from one-site sigma translation, which is tau translation composed with P), and Theta. Klein and P generate A_4, the rotation group of the tetrahedron, acting on (tau^x, tau^y, tau^z) by its real irreducible 3-dimensional representation. Section 1.5 gives the LSM constraint: no unique gapped ground state.

[derived; inputs: the SU(2)_1 description of the Heisenberg chain and its marginal operator (Affleck, Gepner, Schulz and Ziman, J. Phys. A 22, 511 (1989))] At (1,1,1) the infrared is SU(2)_1 with fields g (dimension 1/2), currents J^a_L, J^a_R (dimension 1), and lattice translation acting as the centre g -> -g. In the cyclic stratum sigma translation acts as the centre times the rotation P, and the Klein group acts by pi rotations. Enumerate the perturbations of dimension at most 2 with spin zero:
- tr g (dimerization, dimension 1/2, rotation singlet): odd under the centre, so forbidden by translation.
- tr(g tau^a) (Neel, dimension 1/2, a vector): an invariant combination needs an A_4-invariant vector; the 3-dimensional representation has none (Schur), so forbidden.
- J^a_L J^b_R (dimension 2, marginal), coefficient g_ab: 3 x 3 = 1 + 3 + 5 under SO(3), and 3 and 5 contain no A_4 singlet (5 restricts to 1' + 1'' + 3, with no trivial summand), so g_ab = g delta_ab is forced.
- J^a_L or J^a_R alone (spin +-1) and the traceless part of J^a_L J^b_L (spin 2): no A_4 singlet.
Hence, at every point of the cyclic stratum whose flow enters the SU(2)_1 basin, the only allowed low-dimension coupling is the isotropic marginal g J_L . J_R. For g > 0 it is marginally irrelevant: the fixed point is SU(2)_1 with full emergent SU(2)_L x SU(2)_R and c = 1. For g < 0 it is marginally relevant: a gap opens exponentially slowly with spontaneous dimerization (two ground states), c = 0. There is no third option inside the basin, and in particular no c > 1.

[derived, finite-size diagnostics of the sign of g] With log corrections, and in the normalization where g > 0 is the marginally irrelevant sign (the Heisenberg sign), x_triplet = 1/2 - g/4 + ..., x_singlet = 1/2 + 3g/4 + ... (Affleck et al. 1989), so r_s = x_singlet/x_triplet > 1 decreasing towards 1 like 1/ln L means g > 0, and r_s crossing below 1 means g < 0 (the level-spectroscopy criterion of Okamoto and Nomura, Phys. Lett. A 169, 433 (1992)). The current level stays at x = 1 (r_c near 2).

### 3.4 Numerics on the YX line

[numerical] Momentum-resolved ring ED (cyc_spec.py, ed_sec.py with C155's RingTable; L = 12 and 18 locally, L = 24 on Typhon in sectors m = 0, 1, 4, 12; analysis cyc_analyze.py). Every point from eps = 0.05 to Heisenberg has the SU(2)_1 tower in the same quantum numbers: a unique ground state at momentum 0 (L = 12, 24) or pi (L = 18), an exactly threefold triplet at relative momenta pi + {0, +-2 pi/3} (the three components are permuted by sigma translation, so their degeneracy is exact), a singlet at relative momentum pi, and an exact current sextet at relative momenta +-2 pi/L + {0, +-2 pi/3}. The ZY spectra equal the YX spectra (mirror). Table (x = L (E - E0); v = x_c/2 pi; c_cas from E0/L at consecutive sizes divided by v at the larger size):

| (a,b,g) | r_s at L = 12, 18, 24 | r_c at L = 24 | v (L = 24) | c_cas (12,18) | c_cas (18,24) |
| :-- | :-- | :-- | :-- | :-- | :-- |
| (1, 0.05, 0) | 1.185, 0.918, 0.787 | 1.79 | 0.286 | 3.01 | 1.79 |
| (1, 0.1, 0) | 1.264, 1.090, 1.020 | 1.94 | 0.562 | 1.51 | 1.22 |
| (1, 0.15, 0) | -, -, 1.160 | 2.02 | 0.817 | | |
| (1, 0.2, 0) | 1.386, 1.288, 1.251 | 2.07 | 1.052 | 1.17 | 1.08 |
| (1, 0.3, 0) | 1.471, 1.396, 1.362 | 2.14 | 1.469 | 1.10 | 1.045 |
| (1, 0.5, 0) | 1.567, 1.499, 1.462 | 2.20 | 2.135 | 1.047 | 1.021 |
| (1, 1, 0) | 1.626, 1.555, 1.513 | 2.23 | 3.211 | 1.025 | 1.010 |
| (1, 1, 1) Heisenberg | 1.714, 1.626, 1.574 | 2.27 | 6.335 (exact 2 pi) | 1.016 | 1.004 |

Here r_s = x_singlet/x_triplet and r_c = x_current/x_triplet. Calibration: at Heisenberg the same estimators give c_cas = 1.004 at (18, 24), v within 0.8 percent of 2 pi, and r_s decreasing slowly from above, the known marginally irrelevant behaviour.

Reading the table. [numerical, with the criterion of 3.3]
- For eps >= 0.2 (and on to Heisenberg) r_s > 1 and decreases slowly with L as at Heisenberg, r_c is near 2, and c_cas approaches 1 from above (1.08, 1.045, 1.021, 1.010 at (18, 24)): the SU(2)_1 phase with marginally irrelevant coupling, c = 1.
- For eps = 0.05 the singlet crosses below the triplet between L = 12 and 18 and keeps falling (r_s = 0.918, 0.787), while x_triplet grows (0.744, 0.915, 1.005): the isotropic marginal coupling is negative, which is the spontaneously dimerized phase (two ground states at momenta 0 and pi, exponentially small gap, c = 0 in the infrared). At eps = 0.1, r_s = 1.020 at L = 24 and still falling, so eps = 0.1 is at or just below the transition. The transition eps_c (a point with g = 0, SU(2)_1, c = 1, the analogue of J2c in the J1-J2 chain) is near 0.1 from these sizes; L = 24 runs at eps = 0.07, 0.085, 0.125, 0.175, 0.25 and DMRG dimerization at L up to 192 are in Typhon jobs 17594697 and 17594764 (section 3.6).
- The excess of c_cas above 1 at small eps is the fast A Majorana (velocity 4, c = 1/2) that is gapped only at the scale set by eps (C144 section 1.4): it contributes to the Casimir energy until L exceeds its correlation length. c_cas falls with L at every eps (3.01 -> 1.79, 1.51 -> 1.22).
- The emergent velocity is linear in eps at small eps: v/eps = 5.72, 5.62, 5.44, 5.26 at eps = 0.05, 0.1, 0.15, 0.2 (L = 24), consistent with C141's emergence bound v_em <= C delta M.

Additional L = 12, 18 points (cyc_spec_run2.out), r_s at L = 12, 18: eps = 0.07: 1.218, 0.995; 0.085: 1.241, 1.046; 0.125: 1.299, 1.154; 0.15: 1.331, 1.206 (1.160 at L = 24); 0.175: 1.360, 1.251; 0.25: 1.433, 1.349. The crossing r_s = 1 moves from eps = 0.07 (L = 18) to about 0.095 (L = 24, interpolating 0.05 and 0.1). Extrapolating the decrements of r_s geometrically (ratio near 0.4 per size step at eps = 0.1 and 0.15) gives r_s(infinity) near 0.97 at eps = 0.1 and 1.13 at eps = 0.15, so eps_c is between 0.1 and 0.15 [numerical; the extrapolation is heuristic, since marginal flows are logarithmic, and the A-sector crossover shifts small sizes]. Safe statement: dimerized for eps <~ 0.07 to 0.1, SU(2)_1 for eps >~ 0.15 to 0.2.

### 3.5 The rest of the cyclic stratum

[numerical] cyc_grid.py on 20 points of the fundamental triangle alpha = 1 >= beta >= |gamma| (beta = 0.25, 0.55, 0.85, 1; gamma/beta = -0.9, -0.5, 0, 0.5, 0.9), rings L = 12, 18 (cyc_grid_local.jsonl, cyc_grid_local1.out, cyc_grid_local2.out; the denser 53-point grid is queued in Typhon job 17594765). Two regions appear.
- SU(2)_1 region: every point with gamma/beta >= -0.5 (16 points). Same tower as on the YX line, triplet split below 6e-12 (exact A_4 degeneracy), r_s between 1.23 and 1.95 and decreasing with L, r_c between 2.01 and 2.71, c_S (half-ring entropy slope between L = 12 and 18) between 0.98 and 1.07, c_cas between 0.84 and 1.21 (the extremes are the slowest points, v = 0.48 at (1, 1, -0.5) and v = 0.55 at (1, 0.25, -0.125), where 12/18 Casimir fits are least converged). Marginally irrelevant everywhere tested, so c = 1.
- Ordered region near the XZ-FM edge (gamma/beta = -0.9, and (1, 0.55, -0.495)): the L = 18 ground state is degenerate (gap1 = 0), the ground-state momentum changes, and the "current" level falls far below the triplet (r_c = 0.39): a symmetry-broken, ferromagnet-like phase continuing from the FM Heisenberg vertex (allowed by LSM since it breaks Theta). Its entropy estimators (c_S = 0.53 to 0.66, and a meaningless c_S = 2.58 with c_cas = 48.7 at (1, 0.25, -0.225), where the ground space is a near-degenerate multiplet as at the FM point in C155's calibration) do not describe a CFT. The boundary between the two regions lies between gamma/beta = -0.5 and -0.9 and was not resolved.
- No point shows a tower beyond SU(2)_1, an extra low singlet with r_s well below 1 away from XZ, or c estimates growing with L. In particular there is no sign of a curve of SU(2)_3-type (c = 9/5) or other c > 1 criticality, the only c > 1 option that section 4.3 leaves on this stratum. The grid is coarse (spacing 0.3 in beta), so a narrow feature between grid points is not excluded by it.

### 3.6 Verdict on the YX and ZY directions and on the cyclic stratum

[derived plus numerical] Along YX (and ZY, its mirror) the A Majorana of the XZ point is gapped and the B degrees of freedom form the SU(2)_1 liquid of the period-3 XYZ chain, connected without a transition (same tower, same quantum numbers) to the antiferromagnetic Heisenberg point. Very close to XZ (eps below eps_c, with eps_c between 0.1 and 0.15 from L <= 24) the isotropic marginal coupling changes sign and the chain dimerizes with an exponentially small gap; at eps_c the theory is SU(2)_1 with g = 0. So c = 1 for eps > eps_c, c = 1 at eps_c, c = 0 for eps < eps_c, and c_tot <= 1 everywhere, with v proportional to eps (v = 5.7 eps at small eps). This settles C144's "slow sector with c about 1 at L = 64": it is SU(2)_1 (or its dimerized neighbour), not a c > 1 theory. On the whole cyclic stratum the critical points have c = 1 (SU(2)_1 region, the dimer transition line, and the XZ vertex with c = 1/2 plus degeneracy) or are singular (FM vertex, z = 2), and the ordered and dimerized regions are gapped.

Pending checks (Typhon). Open-chain DMRG at L = 48 has finished for eps = 0.05, 0.1, 0.2, (1,1,0) and Heisenberg (typhon_out/dm_L48.jsonl): parity-resolved fits (1.32, 0.73), (1.42, 0.80), (1.27, 0.76), (1.08, 0.70), (1.04, 0.69), averages 1.03, 1.11, 1.02, 0.89, 0.86, with the strong even/odd alternation typical of SU(2)_1 or dimerized open chains, and centre dimerization 0.052, 0.055, 0.061, 0.117, 0.190 in units of the bond energies (which are 0.65, 0.68, 0.73, 1.19, 1.76 per site). At a single size these do not separate a slow dimerization from SU(2)_1; the L = 96, 144, 192 runs (jobs 17594697 and 17594764) will, through the size dependence of the centre dimerization.

## 4. Codimension and RG count for c > 1 candidates

C68 section 5 counts relevant operators for Ising, tricritical Ising, Potts, Ising x boson, three Majoranas, SU(2)_k and parafermions. This section adds the one ingredient that changes the counts, the action of translation, and draws the conclusion for B_2.

### 4.1 The count

[derived] Let a bond h* flow to a unitary CFT C, and let g_T in Aut(C) be the infrared image of one-site translation (on N there is nothing else). Bonds near h* flow near C with couplings lambda_i(h) to the scalar operators O_i of C; lattice momentum conservation forces lambda_i = 0 unless g_T O_i = O_i. If the matching map h -> (lambda_i) has full rank at h*, the locus of bonds flowing to C (or to its conformal manifold) has codimension

    k(C, g_T) = #{ relevant (Delta < 2), spin-zero, g_T-invariant, non-redundant O_i },      (4.1)

and its dimension inside it is 8 - k minus nothing (exactly marginal invariant operators move along the conformal manifold inside the locus). Lorentz-violating invariant operators of spin +-1 and dimension 1 (uniform currents) shift densities and pitches, not tunings. If the closure of {g_T^n} is a torus T^r (incommensurate action), every operator charged under T^r is forbidden, not only those charged under g_T. A realization needs k <= 8, and the rank condition, which no symmetry argument can certify.

### 4.2 Table (minimum over admissible translation actions, in N)

| CFT | c | translation action | invariant relevant scalars | k |
| :-- | :-- | :-- | :-- | :-- |
| Ising | 1/2 | Z2 (sigma at momentum pi, period-2 order) | epsilon | 1 |
| Ising | 1/2 | trivial (uniform order, liquid-gas endpoint) | epsilon, sigma | 2 |
| compact boson | 1 | torus T^2 (floating) | none | 0 |
| compact boson | 1 | commensurate | at least one locking cosine, unless all are irrelevant | 0 or 1 |
| Ising x boson | 3/2 | Z2 x T^2 (sigma at pi, boson floating) | epsilon (sigma V, epsilon V, sigma J carry momentum) | 1 |
| Ising x boson | 3/2 | trivial x T^2 | epsilon, sigma, and sigma J (dimension 9/8, spin 1) | 2 to 3 |
| SU(2)_2 | 3/2 | maximal torus of SO(4) | Phi^{zz} (the Ising x boson point at one radius) | 1 |
| SU(2)_2 | 3/2 | trivial | tr g, 9 Phi^{ab} | 10 (> 8: impossible in N) |
| two bosons | 2 | T^4 generic | none | 0 |
| Ising^2 x boson | 2 | (Z2 on both) x T^2 | epsilon_1, epsilon_2, sigma_1 sigma_2 (dimension 1/4) | 3 |
| tricritical Ising x Ising | 6/5 | Z2 x Z2 | epsilon, epsilon' (TCI), epsilon (Ising), plus sigma sigma products | >= 3 |
| Z5 parafermion | 8/7 | Z5 | epsilon_1 (4/7), epsilon_2 (12/7), plus chiral terms | >= 2 |

[proved, consequence] Kinematics allows c > 1 in N at codimension 0 (two or more floating bosons) and codimension 1 (Ising x floating boson with a translation-odd order field). The count therefore cannot exclude c > 1; it singles out these two as the only candidates that could occupy open sets or hypersurfaces, and it shows that every other candidate needs at least two or three tunings, i.e. sits on loci of dimension at most 6 that, if they exist, should be found at or near strata where extra symmetry removes relevant operators. This agrees with C68 section 6 ("symmetry-permitted, not excluded").

### 4.3 Where a c > 1 locus would have to sit

[derived] (i) A codimension-0 or 1 locus is large (dimension 7 or 8). A codimension-0 c = 2 region would be bounded by codimension-1 surfaces where one species gaps (Pokrovsky-Talapov or KT boundaries), and those surfaces separate it from c = 1 floating regions; this is pattern P3 of section 5. A codimension-1 Ising x boson surface separates a floating phase from a "floating plus period-2 order" phase; pattern P1. (ii) Loci with k >= 3 have dimension at most 5 in N, but on a symmetric stratum of codimension d_G where the symmetry forbids some of the relevant operators the count is d_G + k_G, which can be smaller. The relevant strata are T and K (d_G = 3), parity (4), cyclic (6). For SU(2)_k with odd k (which carries the Kramers/Klein anomaly) the cyclic stratum would need k_G = 1 (the A_4-singlet thermal operator of dimension 4/(k+2)), i.e. a curve in the 2-dimensional cyclic stratum: section 3.5 scans that stratum and finds no such curve. (iii) In the parity class, the TCI x Ising and Ising^2 x boson candidates are C152's.

## 5. Phase-diagram topology: c = 1/2 and c = 1 loci and how they meet

The coordinator's angle: a c > 1 fixed point h* in B_2 has relevant deformations flowing to lower-c fixed points, so its neighbourhood contains lower-c critical loci glued in a pattern fixed by the CFT. Section 5.1 describes the c = 1/2 and c = 1 loci of N, 5.2 the patterns the candidates of 4.2 require, 5.3 and 5.4 two arguments that the patterns do not form, and 5.5 and 5.6 the two cases topology cannot see.

### 5.1 The lower-c loci in N

[derived, building on C68 and section 2]
- c = 1, codimension 0: floating liquids. Their emergent symmetry is the torus generated by translation (both vertex towers carry incommensurate momenta, (2.3)). Near the U(1) class they are the parents' liquids dressed by off-resonant harmonics (2.3), so they are attached to the U(1) stratum through the lam path. Whether C68's symmetry-free floating point 1282 is attached in this way is being tested numerically (section 6).
- c = 1, codimension 1: Kosterlitz-Thouless boundaries of floating regions at resonances (2.3), where a locking cosine becomes marginal; and commensurate c = 1 lines inherited from symmetric strata (XYZ and cyclic SU(2)_1, Gaussian lines of the parity class), which lie on the boundary of floating regions because a symmetry-breaking field makes the liquid incommensurate rather than gapping it, or gaps it through a locking cosine.
- z = 2, codimension 1: commensurate-incommensurate (Pokrovsky-Talapov) surfaces, C143 stratum 3. Not conformal.
- c = 1/2, codimension 1: Ising surfaces of period-2 order (sigma at momentum pi, so its source is forbidden by translation even in N). The ordered side breaks translation and is gapped. In the spiral picture these come from Neel order of the parent (Delta > 1) or from period-2 locking at the resonance q(theta + pi) + 2 n k_F = pi.
- c = 1/2, codimension 2: Ising critical endpoints of first-order surfaces (uniform order field, both sigma and epsilon tuned), with the emergent Z2 of the Ising fixed point not protected by any lattice symmetry. Near the parity stratum (dim 4) they are the continuation of the parity-class Ising surfaces (dim 3 inside the stratum, where sigma is parity-odd and its source vanishes identically); off the stratum the sigma source is a second tuning, giving a set of dimension 6 that contains the parity-class Ising surface. The same holds near K for antiunitary-odd order fields.

### 5.2 The patterns required by the c > 1 candidates

[derived from 4.2]
- P1 (Ising x floating boson, sigma at pi, k = 1). Its two relevant directions +-epsilon lead to a floating liquid with the period-2 order disordered, and to a floating liquid with spontaneous period-2 order (two vacua, each carrying a c = 1 boson). So P1 needs an open c = 1 region whose ground states break translation spontaneously ("floating + SSB") adjacent to an ordinary floating region.
- P2 (Ising x boson, uniform sigma, k = 2). Its neighbourhood contains a first-order surface inside a c = 1 region (two coexisting floating liquids with different order-field expectation) ending on the c = 3/2 line. The order field must be independent of the boson, since a first-order jump of the boson's own density is a compressibility instability, which in one dimension ends at a z = 2 (ferromagnetic or phase-separation) point, not at an Ising x boson point.
- P3 (two floating bosons, k = 0). An open c = 2 region, bounded by codimension-1 surfaces where one of the two species gaps, beyond which lies a c = 1 floating region (a second Pokrovsky-Talapov onset inside a liquid).
- Every other candidate in 4.2 has k >= 2 and its neighbourhood contains P1, P2 or P3 type sub-loci (for example Ising^2 x boson flows to Ising x boson when one Ising mass is turned on), so it also requires one of these three structures nearby.

So in all cases a c > 1 point needs, somewhere in its neighbourhood, a phase or a coexistence in which a c = 1 liquid lives together with an independent second sector (an ordered Z2 or a second liquid).

### 5.3 One soft species at the instability of a uniform polarized phase

Lemma 3. [proved as algebra; derived as physics at the level of linear spin waves] Let the uniform product state with all spins along n be a classical stationary point of a bond (J, B) and linearize (Holstein-Primakoff, spin 1/2). The quadratic Hamiltonian is sum_k [A_k a_k^dag a_k + (Pi_k a_k a_{-k} + h.c.)/2] with A_k = A_0 + a_c cos k + a_s sin k and Pi_k = pi_0 cos k, and the magnon band is

    omega(k) = a_s sin k + sqrt( f(cos k) ),     f(x) = (A_0 + a_c x)^2 - |pi_0|^2 x^2,      (5.1)

where a_s is proportional to the DM component along n, a_c to the transverse symmetric exchange J_uu + J_vv, pi_0 to the transverse anisotropy (J_uu - J_vv, J_uv + J_vu), and A_0 to the field and J_nn. Suppose the state is stable (omega >= 0 for all k) and marginal (min omega = 0). Then the set of soft momenta is one of: a single k_1; a pair {k_1, -k_1} with a_s sin k_1 = 0, which requires a_s = 0 for k_1 not in {0, pi}; the pair {0, pi}; or an interval (flat band, singular).

Proof. omega(k) = 0 iff g(cos k) = 0 and a_s sin k <= 0, with g(x) = f(x) - a_s^2 (1 - x^2), a polynomial of degree at most 2. Stability for both signs of sin k (each half of the Brillouin zone covers x in [-1, 1]) is g(x) >= 0 on [-1, 1] when a_s != 0, and f >= 0 when a_s = 0. A nonzero polynomial of degree at most 2 that is nonnegative on [-1, 1] vanishes there either at one interior point (a double root) or only at endpoints, so the zero set in x is {x_0} or a subset of {-1, 1}. Each interior x_0 gives one k (the one with a_s sin k < 0) if a_s != 0 and the pair +-k_0 if a_s = 0; endpoints give k = 0 or pi. QED.

[derived] Physical content: at the boundary of a uniform polarized phase the soft sector is (i) one complex incommensurate mode (translation acts on it as an incommensurate phase, so an emergent U(1) appears and the transition is Pokrovsky-Talapov or KT into a single c = 1 floating liquid), (ii) one real mode at 0 or pi (an Ising transition, uniform or period-2), (iii) two real modes at 0 and pi (two Majoranas, c <= 1, the Ashkin-Teller/Gaussian case of C68), or (iv) with a_s = 0, the pair +-k_0, which pairing Pi binds into a single complex Bogoliubov mode, again one boson. Two independent incommensurate modes, which P3 needs, never soften together at such an onset, because a nearest-neighbour bond gives a single cosine harmonic in A_k and Pi_k and so a quadratic g. This is the quadratic, onset-level form of "one species per site". Its limits: it says nothing about instabilities of non-uniform (spiral) references, of states far from any product state, or about bound-state (two-fluid) channels, which contrib/d2_two_fluid studied numerically (two Fermi seas do appear in the two-body problem on a small open set, but a dense canted state preempted the second sea in every case tested).

### 5.4 Translation-breaking order cannot coexist with a single floating boson

[derived, conditional on a one-boson infrared] Suppose the infrared of a region is a single compact boson with translation acting through the torus (2.3). Every local lattice operator O_pi of momentum pi (a Neel or dimer order field) has a continuum expansion in vertex operators e^{i(m vartheta + 2 n phi)} with momenta m(theta + pi) + 2 n k_F; at generic (theta, k_F) none of them has momentum pi, so O_pi has no component in the gapless sector at all, and its expectation value in a translation-broken vacuum would have to come from a gapped sector. A spontaneous period-2 order with a gapless floating boson ("floating + SSB", needed by P1) therefore requires a second, gapped, ordered sector whose own Ising transition is the epsilon direction of P1: an independent species. Within the single-species description (all controlled regions of this report: the U(1) class, the free locus, the lam-neighbourhood of section 2.3, the cyclic stratum) the order field is a vertex of the same boson, which is locked (gapped) whenever the order is long ranged. So P1 cannot form there, and generic junctions of Ising surfaces with floating regions are Lifshitz-type (the Ising surface ends on a Pokrovsky-Talapov surface, z > 1). This is what contrib/d2_multicritical found numerically in XYZ in a field, its DM extension and tilted-field XXZ (every junction non-conformal or c = 1), and what contrib/d2_floating_ising found on six two-parameter cuts in Z2 and N (no neutral level softens anywhere inside a floating liquid; FES c between 0.93 and 1.09 at all 18 followed points).

### 5.5 Case (a): an isolated gap closing with c > 1

[derived for what is covered; open otherwise] A c > 1 fixed point all of whose relevant deformations lead into the same gapped phase (for instance several Majoranas or Dirac fermions whose masses can all be rotated into one another without a symmetry) is an isolated gap closing with no lower-c loci nearby; topology cannot see it. What excludes it:
- on the free locus, the rank bound (C143 section 2.4, C68 section 3): at most two soft Majoranas, c <= 1, at any codimension [proved];
- in a weak-interaction neighbourhood of the free locus, the number of soft Majoranas at a gap closing is fixed by the quadratic part, and four-fermion terms are irrelevant or marginal there, so c <= 1 [derived];
- in the U(1) class (Bethe ansatz, one real-root sea) and in its lam-neighbourhood (section 2.3) [literature and derived];
- at the boundary of every uniform polarized phase (Lemma 3: at most two real or one complex soft mode) [derived at LSW level];
- in the cyclic stratum (section 3) [derived plus numerical];
- near the singular strata (C141 to C144; YX/ZY completed here).
Outside these sets only a direct numerical location of gap closings can decide. C155's adversarial estimator is designed for exactly this (it maximizes min(c_S, c_cas), calibrated to read 1.9 on two decoupled XX chains and 1.39 on cluster-Ising), and the earlier survey of about 103,000 bonds found no c > 1. This run does not exclude case (a) in general.

### 5.6 Case (b): a stable c = 2 floating phase

[derived for the mechanism, open in general] Kinematically a c = 2 phase with two incommensurate modes is stable in N if every vertex operator carries nonzero momentum, which holds when the two Fermi momenta and 2 pi are rationally independent (section 4.1). So the question is dynamical: can a nearest-neighbour qubit chain have two independent incommensurate soft modes?
- No on the free locus (single-harmonic band, at most one Fermi sea) [proved, C143], in the U(1) class [literature], at zero-field XYZ [Baxter], in the cyclic stratum [section 3], in the lam-neighbourhood of the U(1) class [derived, section 2.3: harmonics lock or are irrelevant, they never create a second sea].
- Not at the boundary of a uniform polarized phase: Lemma 3 allows only one complex incommensurate soft mode there. A c = 2 region therefore cannot be entered directly from a polarized phase; it would have to open by a second onset inside a c = 1 floating liquid (P3).
- The candidate second species are a second quasiparticle minimum (excluded at onset by Lemma 3) or bound pairs. The pair channel is pinned to twice the single-particle minimum by U(1) (contrib/d2_two_fluid section 2, exact), and charge-one terms that unpin it lead, in every numerically tested case, to a dense canted state that undercuts the vacuum before the second sea fills.
So the answer is: two independent incommensurate modes are excluded in every controlled region and at the onset of every uniform polarized phase, and every mechanism that could produce them away from those regions has failed in the tests so far; a general exclusion is not available. This is the content of assumption (A6) below.

## 6. Task 2: does every symmetry-free gapless point reduce to a symmetric stratum?

[proved, correction of the premise] The list "U(1), parity, staggered charge, singular" is incomplete. By Theorem 1 the symmetric strata of B_2 are those of the table in 1.4, which adds T and K (both 5-dimensional, C68) and the new cyclic stratum. T and K are not solved classes, so reduction to them does not by itself fix c.

[derived] Emergent symmetry. At every conformal point the infrared image of translation is an internal symmetry of the fixed point: an incommensurate torus U(1) x U(1) for floating liquids, a Z2 for period-2 Ising surfaces, SU(2) x SU(2) on the cyclic stratum's critical region (section 3.3, forced by A_4 and Schur). Ising critical endpoints have the Ising Z2 of the fixed point, which no lattice symmetry protects; it is the second tuning.

[derived] Continuous connection. (i) Floating regions and their KT and C-IC boundaries are connected to the U(1) class along the spiral path of section 2 as long as the harmonics are perturbative. [numerical, incomplete] The path from 1282 (lam = 1) to its parent (lam = 0, XXZ with Delta = -0.414 and h/J_perp = 0.22) was started locally at L = 32 and 64 (path_dmrg.py, path1282.jsonl). The window-averaged entropy Sbar at L = 32 decreases smoothly with lam (0.928, 0.909, 0.866, 0.836, 0.823 at lam = 0, 0.25, 0.5, 0.75, 1), with no kink. Only lam = 0 and 0.25 finished at L = 64 before the local machine became unusable (load average above 100 from other runs): c_mid = 6 [Sbar(64) - Sbar(32)]/ln 2 = 0.99 at lam = 0 (the parent, a c = 1 calibration) and 0.51 at lam = 0.25. The low value at lam = 0.25 is unexplained at these sizes (an incommensurate liquid with a long crossover, or a locked pocket near a resonance (2.3), are both possible); C68's independent finite-entanglement value at lam = 1 is c = 1.008. The full path at L = 32, 64, 128 is in Typhon job 17594765. So the connection of 1282 to the U(1) class is not yet demonstrated. (ii) Period-2 Ising surfaces descend from parent Neel order or period-2 locking, i.e. from the U(1) class through the same path. (iii) Ising critical endpoints near the parity stratum contain the parity-class Ising surfaces (5.1). (iv) The cyclic stratum's critical region is connected to the Heisenberg point (section 3).

[heuristic, not derivable in general] The strong statement "every gapless point of a symmetry-free bond is continuously connected, along its critical locus, to a point of a symmetric stratum" cannot be derived with these tools: the lam path is canonical but leaves the perturbative regime at lam of order one for generic bonds (harmonic fraction 0.61 at 1282), and a critical locus may in principle close up inside N without touching any stratum. Not refuted either: no counterexample is known, every followed case (1282, the cyclic stratum, the d2 cuts) is connected. The statement is therefore kept as assumption (H) of section 7.

## 7. Conditional derivation of c <= 1 on all of B_2

Closure lemma. [derived; inputs: Zamolodchikov's c-theorem and (A2)] Let Sigma be a connected component of the set of bonds flowing to one conformal manifold. Then c is constant on Sigma (marginal deformations preserve c), and at every conformal frontier point h' of Sigma, c(h') >= c(Sigma), because bonds of Sigma near h' are weak relevant perturbations of h' (by (A2)). At a non-conformal frontier point (a singular stratum), c(Sigma) is bounded by the emergent central charge of that point's neighbourhood. So upper bounds on c propagate from frontiers inward.

Assumptions.
- (A1) [assumed, as in C143] Every bond is gapped (possibly with a finite degenerate ground space), conformal with finitely many unitary sectors, or singular (z > 1 or an extensive degeneracy).
- (A2) [assumed, C141's weak c-theorem] Near a conformal bond h0, every bond in some neighbourhood with a conformal infrared has c <= c(h0).
- (A3) [established] c <= 1 at every conformal point of the U(1) class, the free locus and zero-field XYZ.
- (A4) [pending in C152 and C153] c <= 1 on the parity class and on the staggered-charge class.
- (A5) [derived and numerical, C141 to C144 and section 3] Near every singular point the emergent central charge is at most 1 in every direction, now including the YX and ZY directions at XZ (c = 1, SU(2)_1, or 0).
- (A5') [derived plus numerical, section 3] c <= 1 on the cyclic stratum.
- (H) [assumed, the remaining input] Every connected critical component Sigma has a point of its closure in the controlled set (A3) to (A5'), or in the lam-neighbourhood of the U(1) class of section 2.3.

Conclusion. [derived, conditional on (A1) to (A5') and (H)] Every conformal bond of B_2 has c <= 1. Proof: by (H) the closure of its component Sigma meets the controlled set at some point h'. If h' is in Sigma, c(Sigma) = c(h') <= 1. If h' is a conformal frontier point, c(Sigma) <= c(h') <= 1 by the closure lemma. If h' is singular, c(Sigma) <= 1 by (A5). QED.

What (H) contains. (H) fails exactly when there is a critical component of B_2 that never comes near any controlled region: an isolated c > 1 gap closing (case (a), section 5.5), a c = 2 floating region or a "floating + SSB" region whose closure avoids the controlled set (case (b) and pattern P1), or more generally a c > 1 island. By 4.2 and 5.2 every such island needs two independent soft sectors coexisting at a nearest-neighbour qubit bond, which is excluded in every controlled region (5.3 to 5.6) and was not found by any numerical search to date. Equivalently, (H) can be replaced by the one-species statement
- (A6) [assumed] at a conformal point of a nearest-neighbour qubit chain outside the controlled set, the infrared has a single boson or at most two Majoranas,
which gives c <= 1 directly. (H) is weaker than (A6) and is what the topology argument needs.

## 8. Consistency checks

- Theorem 1 against known strata: for O = 1 it returns the uniform Klein group of XYZ (plus Theta at B = 0); for O a rotation about a principal axis it returns C67's parity class (B along the axis) and C68's K stratum (B perpendicular); for O a 2 pi/3 body-diagonal rotation it returns the cyclic stratum, whose symmetries cyc_check.py verifies exactly. The counts reproduce the dimensions 4 (parity) and 5 (T, K) of C67 and C68.
- Proposition 2 limits: (alpha, beta, gamma) = (1,1,1) gives XZ + YX + ZY, which is also XXZ about (1,1,1) with J_perp = -1/2, J_parallel = 1 plus a DM vector of length sqrt 3/2 along (1,1,1); gauging away the DM gives J' = 1 and Delta = 1, the Heisenberg chain, as the tau map says. The ED velocity at that point is within 0.8 percent of the exact 2 pi (H = sum sigma.sigma = 4 sum S.S).
- (1, 0, 0) is the XZ point: the tau chain is the three-colour Kitaev chain whose bond terms form a single anticommuting chain of L Majorana bilinears, leaving L/2 free qubits, which is C143's 2^{L/2 - 1} degeneracy.
- YX and ZY spectra coincide (mirror), as the beta <-> gamma exchange of 3.2 requires.
- Lemma 3 against exact spectra: lsw_check.py, error 0.033 -> 0.008 and 0.20 -> 0.05 as |B| goes from 8 to 32 (linear in 1/|B|).
- Spiral frame: J(1) = J to 6e-16 and the parent's U(1)_m invariance to 2e-16 (spiral.py). At C68's anchor (S = diag(1,1,0.4), DM 0.5 along z) the parent is the anchor itself and the harmonics vanish, recovering C68 section 4.
- Casimir estimator calibration: Heisenberg c_cas = 1.016 at (12, 18) and 1.004 at (18, 24).

## 9. Addendum: open-chain DMRG at L = 96 (Typhon 17594697, typhon_out/dm_L96.jsonl)

[numerical] Centre dimerization at L = 48 and 96, and its ratio D(96)/D(48): eps = 0.05: 0.0521, 0.0407 (0.78); eps = 0.1: 0.0546, 0.0379 (0.69); eps = 0.2: 0.0612, 0.0408 (0.67); (1,1,0): 0.1169, 0.0775 (0.66); Heisenberg: 0.1897, 0.1257 (0.66). In SU(2)_1 the boundary-induced dimerization decays like L^{-1/2} up to logarithms, which the Heisenberg calibration shows as a ratio 0.66; eps = 0.2 and (1,1,0) decay the same way, eps = 0.1 slightly slower, and eps = 0.05 clearly slower (0.78, against 1 for saturated long-range order). This supports section 3.6: eps = 0.05 lies in the dimerizing regime (correlation length still comparable to L), eps >= 0.2 in the SU(2)_1 phase, eps = 0.1 near eps_c. Parity-averaged Calabrese-Cardy fits at L = 96: 0.985 (0.05), 1.09 (0.1), 1.01 (0.2), 0.92 (1,1,0), 0.90 (Heisenberg, whose open-chain log corrections pull the average below 1 at this size). The midpoint estimator 6[S(96/2) - S(48/2)]/ln 2 gives 1.16, 1.30, 1.20, 1.07, 1.05; its excess at small eps decreases with eps and is the A Majorana crossover. No estimator shows c > 1 growing with L. The L = 144 runs are still in job 17594697.
