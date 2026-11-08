# C159 derivations: c <= 1 on the time-reversal strata T and K

This run was done by a Claude Opus agent (2026-11-07), not by Codex.

Standard: theoretical physics. Tags: [proved] (complete argument, or an exact identity checked to machine precision), [derived] (physics-standard derivation with named inputs), [numerical] (sizes, scripts and outputs named), [heuristic], [assumed]. A tag applies to the paragraph or displayed claim it opens.

Conventions. H_L = sum_j [ sigma_j^T J sigma_{j+1} + B . sigma_j ], Pauli matrices, J = S + A(D) with S symmetric and A(D)_ab = eps_abc D_c, as in C68 and C154. Theta = prod_j (i Y_j) K is Kramers time reversal (Theta^2 = -1 on a site). Bosonization follows C68 section 4: for a chain with a U(1) about an axis m, sigma^m ~ (2/pi) d_x phi + a (-1)^j cos 2phi, sigma^+ ~ e^{i vartheta} [ b (-1)^j + c cos 2phi ], dimer (-1)^j sigma_j . sigma_{j+1} ~ sin 2phi, vertex V_{p,n} = e^{i(p vartheta + 2 n phi)} of dimension x_{p,n} = p^2/(4K) + n^2 K. At zero magnetization one-site translation acts as vartheta -> vartheta + pi, phi -> phi + pi/2 (C68). Time reversal acts as phi -> -phi + pi/2, vartheta -> vartheta + pi (section 2.2 checks this against the spin expansion).

## 1. The strata and their normal forms

### 1.1 T

[proved, C68 and C154 Theorem 1] T is the set of bonds with B = 0 whose only on-site symmetry (uniform or site-dependent) is Theta. In the principal frame of S,

    J = diag(J_x, J_y, J_z) + A(D),   D = (D_x, D_y, D_z),   B = 0,      (1.1)

six numbers modulo scale, dimension 5. Generic points have distinct J_a and three nonzero D_a. Boundary strata of T (more symmetry): zero-field XYZ (D = 0, uniform Klein group, dim 2), parity cap T (D along one principal axis, uniform pi rotation, dim 3), U(1) cap T (XXZ plus D along the symmetric axis, dim 2), the cyclic stratum of C154 (dim 2), the free locus cap T (J_z = 0 with D along z, dim 2), and the Heisenberg points.

### 1.2 K

Lemma 1 (normal form of K). [proved] A bond lies in K if and only if, up to a common rotation, it has the form

    J = diag(J_x, J_y, J_z) + D_z A(e_z),   B = (B_x, B_y, 0),   (B_x, B_y) != 0,      (1.2)

with no further symmetry; the antiunitary is A = Theta R_pi(z), which acts on Pauli vectors as the reflection Q = diag(1, 1, -1).

Proof. C68 gives the conditions S n || n, D || n, B . n = 0 for the single improper symmetry with normal n. Take n = e_z. Then S has e_z as an eigenvector, so its xy block is a symmetric 2 x 2 matrix, and a rotation about z (which preserves D || z and the plane of B) diagonalizes it. Conversely Q^T diag Q = diag, Q A(D) Q^T = det(Q) A(Q D) = -A(-D) = A(D) for D || z, and Q B = B for B in the plane, so (1.2) is invariant under the antiunitary. The count is (J_x, J_y, J_z, D_z, B_x, B_y) modulo scale, dimension 5. QED.

So K is the family "XYZ with a DM vector along one principal axis and a field perpendicular to it". This is C155's coordinate set RA. Its exchange part (B = 0) is the parity cap T stratum; the field breaks both the unitary R_pi(z) and Theta, and keeps their product.

[proved] Two sheets. Write the xy block of J as M = [[J_x, D_z], [-D_z, J_y]], det M = J_x J_y + D_z^2. In the polar decomposition J = S O of C154 (O in SO(3)):
- sheet K_I, det M > 0: O is a rotation by theta in (0, pi) about z, the spiral axis is m = e_z, and the parent of C154 section 2 is an XXZ chain about z in the field m . B = 0;
- sheet K_II, det M < 0: M = S_M O_M with O_M a reflection of the plane across a line u, O = diag(O_M, -1) is the pi rotation about u, and the parent is XXZ about the in-plane axis u in the field u . B, generically nonzero.
The free locus does not meet the interior of K (a quadratic Jordan-Wigner form needs the field along the parity axis).

### 1.3 Literature models inside T and K

[proved for membership; literature for the phase diagrams]
- T contains the uniform Kitaev-Gamma chain (in its U6 frame it is the cyclic bond (alpha, beta, gamma) = (-K, -Gamma, -Gamma), C156 section 1), the Kitaev-Heisenberg-Gamma chain, and Kitaev-Gamma chains with DM. For KHG, khg_check.py maps the alternating chain (x, y bonds, K S^x S^x + J S.S + Gamma(S^y S^z + S^z S^y) and its y counterpart) to the uniform bond J' = J0 Q with J0 = -[[K + J, 0, 0], [0, G or J, J or G], [0, J or G, G or J]] and Q the cyclic permutation; the spectra of 6-site rings agree to 2e-14 for three parameter sets, and the only uniform symmetry of J' is Theta (three distinct singular values, so the test is exhaustive). So KHG is a 2-dimensional family in the interior of T, not in the cyclic or parity strata. Yang, Nocera and Affleck (PRR 2, 033268 (2020), arXiv:2004.12954) report nine phases: four Luttinger liquids (c = 1), ferromagnetic, Neel, a distorted spiral and two D_3-breaking ordered phases. Every ordered phase breaks Theta, as LSM requires, and no c > 1.
- K contains XXZ with a DM vector along z in a field perpendicular to it (Garate and Affleck, PRB 81, 144419 (2010); Chan, Jin, Jiang and Starykh, PRB 96, 214441 (2017)): one Luttinger liquid and two Ising-ordered phases, c <= 1. It contains XXZ in a tilted field (J = diag(1, Delta, 1) and the field in the 1-2 plane after relabelling, with reflection normal along the third axis): Neel, paramagnetic and incommensurate Luttinger phases (as summarized in arXiv:1801.10509, found by search, not read in full). Its rank-one boundary contains the Kitaev chain in a [110] field (C156, C158), where C158 finds c = 1.01 +- 0.03 at hp = 0.45 by finite-entanglement scaling.
- The bond-alternating Kitaev-Gamma chain of Luo, Zhao, Wang and Kee (arXiv:2012.03382), which has a c = 1/2 region near the antiferromagnetic Kitaev point, is not in B_2 (two-site cell with integer spin per cell). Section 2.3 explains why such a region cannot occur in T.

## 2. Symmetry, anomaly and the allowed infrared theories

### 2.1 LSM

[proved, literature] On T, Ogata and Tasaki (Commun. Math. Phys. 372, 951 (2019)) forbid a unique gapped ground state that preserves Theta and translation, for a half-odd-integer spin per unit cell. On K the antiunitary squares to +1 on a site and there is no such constraint: the field-polarized product state is a unique gapped ground state at large |B|.

[numerical, side remark] The same obstruction appears in matrix-product form. Damped VUMPS (C153) with a one-site or three-site cell does not converge at the T anchor of C68 (err 0.4 to 0.6 after 300 to 600 iterations, a degenerate dominant transfer eigenvalue, energies that drift), while the two-site cell is the standard choice (C158). An injective MPS with an odd cell cannot be Theta-symmetric for a Kramers site, so a Theta-symmetric critical ground state is approximated by a non-injective MPS with an odd cell. Stage V of the Typhon job therefore uses two-site cells only.

### 2.2 Lemma A: an exact antiunitary makes one field of a floating boson commensurate

Lemma A. [derived; inputs: the automorphisms of the c = 1 compact boson at generic radius, and the identification of lattice symmetries with automorphisms of the infrared theory] Let the infrared of a translation-invariant chain with an exact antiunitary symmetry A (commuting with translation T_1) be a single compact boson (phi, vartheta), with T_1 acting by a shift (a, b) of (phi, vartheta) (that is, in the identity component). Then A reverses exactly one of the two fields, and the shift of the reversed field is 0 or half a period. The other shift is unconstrained.

Proof. At generic radius (at special radii the automorphism group is larger, and the conclusion extends to a whole connected liquid region by continuity of the translation action) the unitary automorphisms are the shifts U(1) x U(1) and the charge conjugation C (phi, vartheta) -> (-phi, -vartheta). Fix the antiunitary T0: phi -> phi, vartheta -> -vartheta, which respects [phi(x), d_y vartheta(y)] = i delta(x - y) because it maps i to -i. Every antiunitary automorphism is T0 times a unitary one, so it is T0 o shift (reverses vartheta) or T0 o C o shift (reverses phi): exactly one field is reversed. Say A reverses vartheta. The shift of vartheta by b is generated by Q_vartheta, proportional to the integral of d_x phi, which A leaves invariant, so A e^{i b Q_vartheta} A^{-1} = e^{-i b Q_vartheta}, the shift by -b; the phi shift is generated by the integral of d_x vartheta, which A reverses, and the antilinearity restores the sign, so a is preserved. Hence A T_1 A^{-1} = T_1 requires shift(a, -b) = shift(a, b) on local operators, that is 2b = 0 modulo the period of vartheta. The case where A reverses phi is the same with the roles exchanged. QED.

Consistency check of the time-reversal action used in the conventions. [proved, algebra] Theta flips every spin component. sigma^m ~ (2/pi) d_x phi + a (-1)^j cos 2phi must change sign, so phi -> -phi + c with cos(2phi - 2c) = -cos 2phi, c = pi/2. Then sin 2phi -> sin(pi - 2phi) = sin 2phi: the dimer field is Theta-even, as it must be (sigma_j . sigma_{j+1} is Theta-even). cos 4phi -> cos 4phi and sin 4phi -> -sin 4phi. For the transverse field Theta sigma^+ Theta^{-1} = -sigma^- forces e^{i vartheta} -> -e^{-i vartheta}, which with antilinearity is vartheta -> vartheta + pi. So Theta reverses phi, as Lemma A allows.

Application to the strata. [derived from Lemma A]
- T. Theta reverses the spin density along any axis, so in the spiral frame (C154 section 2, axis m) it reverses phi. Hence the phi shift is 0 or pi/2, and only vartheta can be incommensurate. A shift 0 would leave a commensurate theory with only an incommensurate vartheta shift and a Theta-even vertex cos(2phi - c) at zero momentum. Its anomaly would be the same at every radius (anomalies are discrete and constant on the conformal manifold), and at small K that vertex is relevant and pins phi with a unique minimum, a symmetric unique gapped state, contradicting section 2.1. So in T the phi shift is pi/2, as for the zero-field XXZ chain. Equivalently: a floating liquid in T has zero magnetization along the spiral axis and a commensurate k_F = pi/2. This is consistent with the spiral frame, where the parent of every T bond is at zero field (m . B = 0).
- K_I. A = Theta R_pi(z) reverses tau^z (the spiral axis), so it reverses phi; the phi shift is 0 or pi/2. On the branch connected to the zero-field parent it is pi/2.
- K_II. The spiral axis u is perpendicular to z, A leaves tau^u invariant (Theta and R_pi(z) both flip it), and A reverses vartheta (it fixes tau^v with v = z x u, flips tau^z, and is antilinear, so e^{i vartheta} -> e^{i vartheta} as an operator identity, that is vartheta -> -vartheta). So the vartheta shift is 0 or pi: the transverse correlations are commensurate, while the magnetization along u, induced by the field, is free and incommensurate.
- N (no antiunitary): both shifts free, translation acts through a 2-torus (C154 section 2.2).

Name. A boson in T or K is "half-floating": one incommensurate angle and one commensurate pair of vertex operators. The closure of {T_1^n} is U(1) x Z_2 (shift pi/2 or pi) or U(1) (shift 0), never the 2-torus of N.

Consequences.
- [derived] In T, the commensurate vertices are e^{+-2 i phi} at momentum pi (the Theta-odd Neel component cos 2phi along m and the Theta-even dimerization sin 2phi) and the umklapp cos 4phi at momentum 0 (Theta-even, dimension 4K). Every other vertex carries an incommensurate momentum p(pi + alpha) + n pi with p != 0. So the only symmetric nonderivative operator of a T-floating liquid with dimension below 2 is cos 4phi, and a T liquid is stable if and only if K > 1/2. At K = 1/2 it undergoes a Kosterlitz-Thouless transition into Neel order along m (Theta broken, period 2) or into a dimerized phase (translation broken, Theta preserved), according to the sign of the umklapp coupling, with a Gaussian c = 1 line between these two ordered phases. This is the zero-field XXZ and J1-J2 scenario, now for every half-floating liquid in T.
- [numerical, check] halffloat_test.py, ring ED at L = 12 (scan of the structure factor S_mm(k) of the spin component along the spiral axis): the peak sits at k = pi exactly at the C68 T anchor (1.1, 0.9, 0.4; D = (0.1, 0.2, 0.5)), at the K anchor (field (0.1, 0.2, 0)) and at zero-field U(1); it moves to 7 pi/6 at the N anchor (field (0.3, 0.4, 0.5), ground-state magnetization along m nonzero) and at the U(1) chain in a field. At this size the check is weak (resolution 2 pi/12) and amounts to "no magnetization along m in T and K_I", which is what Lemma A predicts.
- [derived] In K_I with the phi shift pi/2 the same list holds (cos 4phi, no anomaly). With the phi shift 0 the A-even part of e^{2 i phi} is allowed at zero momentum with dimension K, so such a liquid is stable only for K > 2. In K_II the A-even part of e^{2 i vartheta} (dimension 1/K) is allowed at zero momentum (vartheta shift pi), so a K_II liquid is stable only for K < 1/2. This is the familiar statement that transverse anisotropy gaps the field-induced liquid of XYZ in a longitudinal field unless the liquid is strongly repulsive.

### 2.3 Lemma B: no c < 1 critical point with full symmetry in T

Lemma B. [derived for the Ising, tricritical Ising and Potts models; assumed, standard, for the rest of the unitary c < 1 list] Let h* in T be a non-singular conformal point (A1 of C154: finitely many sectors, no extensive degeneracy) at which neither Theta nor translation is spontaneously broken. Then its infrared CFT has c >= 1.

Proof for Ising. Suppose the infrared is the Ising CFT, with translation and Theta acting as invertible (unitary, respectively antiunitary) automorphisms; they are invertible because the lattice operators are, and they map local operators to local operators. Every such automorphism fixes epsilon: it maps sigma to +-sigma (sigma is the unique primary of dimension 1/8) and must preserve the real OPE coefficient C_{sigma sigma epsilon} = 1/2 in sigma x sigma = 1 + epsilon, so epsilon -> +epsilon (an antiunitary map preserves real coefficients). Hence the uniform epsilon coupling is allowed in T, and the lattice can move along a direction with a nonzero epsilon component (the rank condition of C154 section 4.1, generic). On the disordered side of epsilon, with or without an allowed sigma source, the Ising theory flows to a unique gapped vacuum preserving every symmetry. That lattice point would then have a unique symmetric gapped ground state, which section 2.1 forbids. So the Ising CFT cannot be the full infrared of a symmetric critical point of T. The same argument works whenever the CFT has a relevant operator invariant under all its invertible symmetries whose one sign flows to a unique gapped vacuum: for tricritical Ising and the A-series minimal models this is the Landau-Ginzburg mass phi^2 of the order parameter, for three-state Potts the thermal field (S_3 singlet, appearing with a real coefficient in sigma x sigma^dagger). For the remaining unitary minimal models we use the classification of their invertible unitary symmetries (Z_2 for all but six models, S_3 for the critical and tricritical Potts models, none for four E_7/E_8-type models; Ruelle and Verhoeven, Nucl. Phys. B 535, 650 (1998)) and the statement, quoted as known in arXiv:2207.06433, that these symmetries are anomaly free [assumed for the antiunitary extensions]. This agrees with the bootstrap evidence that the minimal theories carrying spin-1/2 LSM anomalies of Z_2 x Z_2 type are c = 1 bosons (Lanzetta and Fidkowski, PRB 107, 205137 (2023)) and with the c >= 1 bound for transitions between Z_n x Z_n bosonic SPT phases (Tsui, Huang, Jiang and Lee, Nucl. Phys. B 919 (2017), OSTI 1349311). QED.

Corollaries. [derived]
- Every critical point of T is one of: (a) a symmetric CFT carrying the LSM anomaly, c >= 1; (b) a critical theory inside a spontaneously broken background (degenerate vacua, each hosting the CFT), which can have c = 1/2 (an Ising transition between two phases that both break Theta or both break translation); (c) singular (the Kitaev point XZ, c = 1/2 plus a flat band, C143; the ferromagnetic Heisenberg point, z = 2).
- In particular a "disordered" region bounded by Ising lines in T is never trivially gapped. Checks: (i) on the free locus cap T (H = J_x XX + J_y YY + D(XY - YX), field zero) the BdG band 4D sin k +- 2 sqrt(t^2 cos^2 k + delta^2 sin^2 k) (C68 section 3, t = J_x + J_y, delta = J_x - J_y) can have a single Majorana zero only at k = 0 or pi, where the square root equals |t|, so such a zero needs t = 0 and then closes k = 0 and k = pi together (two Majoranas, c = 1); zeros at k != 0, pi occur only in the pocket regime 4D^2 > delta^2, as the two edges of one pocket (c = 1, C68 section 3); so every gapless point of the free locus cap T has c = 1 and never 1/2, and the gapped points (4D^2 < delta^2, t != 0) are Neel or ferromagnetic in spin language, breaking Theta [proved from the free dispersion]. (ii) C144's c = 1/2 line on the Ashkin-Teller plane of XZ (field zero) lies inside a wedge in which the A sector is strongly dimerized (a_stag up to -0.35), so it is of type (b). (iii) Peotta et al. (J. Stat. Mech. P09005 (2014)) find at h_z = 0 (parity cap T) two Ising lines with a "disordered" phase between two ferromagnets; in spin language at zero field the fermionic band insulator is the twofold degenerate ferromagnet along z (empty or full band), which breaks Theta, so both Ising lines are of type (b). (iv) The c = 1/2 region of the bond-alternating Kitaev-Gamma chain is allowed because that chain has no Kramers LSM (integer spin per cell); its uniform B_2 relative, the cyclic stratum, has only SU(2)_1 (c = 1), dimerized and ordered phases (C154).
- If c <= 1 also holds (section 6), every symmetric conformal point of T has c = 1 exactly.
- Near XZ (relevant to C157). At XZ the A Majorana is critical and the B flat band carries the LSM anomaly. On the hypersurface where the A mass is tuned to zero, the B sector must therefore either break a symmetry (dimerize or order, c_tot = 1/2 with degenerate vacua), or combine with A into the anomalous two-Majorana theory (c_tot = 1, as for XZ + ZX), or carry c_B >= 1 by itself, giving c_tot >= 3/2. The last option is the RZFA-type "Ising plus Luttinger" B chain that C157 tests; its kernel analysis places the B chain on the attractive side, where RZFA find Ising (c = 1/2), a tricritical Ising point (c = 7/10) and a dimerized phase, and no c = 3/2.

### 2.4 Anomaly-matching theories with c > 1

[derived; inputs: SU(2)_k data (Affleck and Haldane 1987), operator dimensions] For each candidate, count the relevant spin-zero operators allowed by Theta and translation (k_T) and decide whether it can carry the LSM anomaly.
- SU(2)_2 (c = 3/2). With translation acting as the centre, the centre-SO(3) mixed anomaly is k mod 2, so SU(2)_2 does not carry the Kramers LSM anomaly; it would need a second anomalous sector (c >= 5/2). Its spin-1 primary Phi^{ab} (dimension 1, 9 components, centre-even, Theta-even) gives k_T >= 9 > 5. Not reachable in T.
- SU(2)_3 (c = 9/5) carries the anomaly but has Phi^{ab} of dimension 4/5, again 9 Theta- and translation-even components: k_T >= 9 > 5. Not reachable in the interior of T; on the cyclic stratum A_4 leaves one singlet (k = 1, a curve of the 2-dimensional stratum), which C154's grid did not find.
- SU(2)_1 x Ising (c = 3/2). epsilon, plus sigma times tr g or tr(g sigma^a) (dimension 5/8, momentum 0 when sigma sits at pi), with one or three Theta-allowed components: k_T >= 2.
- Ising x half-floating boson and two half-floating bosons: section 4.
So the only c > 1 candidates that can sit on loci of dimension >= 4 in T are built from half-floating bosons and Ising sectors, and section 4 counts them.

## 3. Controlled regions of T and K

### 3.1 Spiral frame

[derived, C154 section 2 specialized] For a T bond the spiral rotation about m maps H to a zero-field XXZ parent (J_perp = s_perp, J_par = s_m) plus harmonics of U(1)_m charge +-1 (momentum +-theta in the tau frame) and +-2 (momentum +-2 theta). Because the parent is at zero field, the resonance condition (2.3) of C154, q(theta + pi) + n pi = 0 mod 2 pi, reduces to theta in {0, pi/2, pi}: the harmonic with q = +-1 resonates with e^{i vartheta} cos 2phi (dimension 1/(4K) + K) at theta = 0 and with the staggered e^{i vartheta} (dimension 1/(4K)) at theta = pi; q = +-2 resonates with e^{2 i vartheta} (dimension 1/K) at theta = 0, pi and with e^{2 i vartheta} cos 2phi at theta = pi/2. Off these three values every harmonic is irrelevant at small amplitude, and the infrared is the parent's liquid with K > 1/2, half-floating (vartheta shift pi + alpha, phi shift pi/2). Near the resonances the harmonics lock commensurate magnetic orders, all of which break Theta (LSM-allowed). The same holds on K_I (axis z, parent at zero field) with the field as the q = +-1 harmonic (rotating with pitch theta). On K_II the parent is XXZ about u in the field u . B, and the harmonics are the A-allowed transverse anisotropy, which pins vartheta unless K < 1/2 (section 2.2).

### 3.2 Neighbourhood of the antiferromagnetic Heisenberg point (in T)

[derived; inputs: SU(2)_1 with the marginal current-current coupling (Affleck, Gepner, Schulz and Ziman 1989), the chiral-rotation removal of a uniform DM vector (as in Garate and Affleck 2010)] Near (1,1,1) a T bond is SU(2)_1 plus: tr g and tr(g sigma^a) (dimension 1/2), forbidden by translation (centre); the uniform currents, of which Theta allows only D . (J_L - J_R) (the DM vector; D . (J_L + J_R) is a uniform field, Theta-odd); and the 9 marginal couplings g_ab J^a_L J^b_R. The DM term is removed by the chiral rotation that is the continuum form of the spiral frame, after which only the part of g_ab invariant under rotations about D survives (the rest oscillates). The result is an XXZ-type marginal flow about D-hat: an incommensurate Luttinger liquid (c = 1), Neel order along D-hat (Theta broken, gapped) or a dimerized phase (negative isotropic part). For D = 0 this reproduces Baxter: three distinct couplings flow to Neel order along the largest, two equal largest give c = 1. So every T bond in a neighbourhood of the Heisenberg point has c in {0, 1}.

### 3.3 Neighbourhood of the cyclic SU(2)_1 region (in T)

[derived] On the cyclic stratum the infrared is SU(2)_1 with one-site translation acting as the centre times the order-3 rotation P (C154 section 3.3). Moving into T keeps translation exact. tr g is odd under the centre; the Neel vector transforms by -P, which has no eigenvalue 1 (eigenvalues -1, -e^{+-2 pi i/3}); so no relevant operator is allowed. The allowed dimension-1 currents are (J_L - J_R) . n_111 (Theta-even, P-invariant), a DM vector along (1,1,1), removed by the chiral rotation; the marginal couplings are the P-invariant ones (a 3-dimensional space) restricted to those invariant under rotations about n_111. So off the cyclic stratum the SU(2)_1 region continues into c = 1 liquids, Neel order along (1,1,1) or dimerization, with c in {0, 1}.

### 3.4 Neighbourhood of the ferromagnetic Heisenberg point and the onsets of polarized states

[derived, C154 Lemma 3 specialized to B = 0] For B = 0 the uniform product state along n is a classical stationary point when n is an eigenvector of S. At the instability of such a ferromagnet Lemma 3 of C154 gives a single soft species (one incommensurate complex mode, one real mode at 0 or pi, the pair {0, pi}, or a pair +-k_0 bound into one boson). So at the boundary of every ferromagnetic phase of T and K the emerging critical theory has c <= 1. For small perturbations of the ferromagnetic Heisenberg point the classical anisotropy energy n^T S' n is first order in the perturbation and the DM-driven spiral gain is second order, so a generic small T perturbation (distinct eigenvalues of S') gives a gapped Theta-broken ferromagnet. Liquids with large Luttinger parameter (the easy-plane case and its DM-twisted relatives) lie near this region; section 4 explains why they matter.

### 3.5 Summary of controlled sets

[derived and literature] c <= 1 holds at the physics standard on: the U(1) cap T and U(1) cap K neighbourhoods (3.1), zero-field XYZ and its neighbourhood in T (Baxter plus 3.1), the AF Heisenberg neighbourhood (3.2), the cyclic stratum and its T-neighbourhood (C154 and 3.3), ferromagnetic onsets (3.4), the free locus cap T (2.3 (i)), the parity class (C152, pending, A4 of C154), the XZ neighbourhood in the directions studied (C144, C154, C158, C157 pending), and the literature families KHG, Kitaev-Gamma (+DM), Garate-Affleck/Chan et al., tilted-field XXZ, Kitaev [110].

## 4. Codimension of c > 1 candidates in T and K

The count of C154 section 4.1 (k = number of relevant, spin-zero, symmetry-invariant, non-redundant couplings) changes in T and K because the boson is half-floating (Lemma A). Write x for the scaling dimension of the minimal commensurate vertex V_c of a half-floating boson (x = K for V_c = e^{2 i phi} in T and K_I, x = 1/(4K) for V_c = e^{i vartheta} in K_II), and V_c at momentum pi (shift half-period) unless stated.

### 4.1 Table

[derived]

| candidate | c | allowed relevant couplings (T, or K with V_c at pi) | k |
| :-- | :-- | :-- | :-- |
| single half-floating boson | 1 | V_c^2 (dimension 4x) | 0 if x > 1/2 |
| Ising x boson, sigma at momentum pi | 3/2 | epsilon; sigma V_c (the combination with matching Theta or A parity, dimension 1/8 + x) | 2 if x < 15/8, 1 if x > 15/8 |
| Ising x boson, sigma uniform and odd under Theta (T) or A (K) | 3/2 | epsilon; sigma d_x phi (dimension 9/8; d_x phi is the odd density) | 2 |
| Ising x boson, sigma uniform and even | 3/2 | epsilon, sigma, sigma d_x vartheta | 3 |
| two half-floating bosons, both V_c at pi | 2 | V_c1 V_c2 (dimension x1 + x2, momentum 0), V_c1^2, V_c2^2 | 0 if x1 + x2 > 2 and x1, x2 > 1/2 |
| two bosons, one V_c at momentum 0 | 2 | that V_c itself (dimension x2) | 0 only if x2 > 2 |
| SU(2)_2, SU(2)_3 | 3/2, 9/5 | Phi^{ab} | >= 9 (section 2.4) |

Steps behind the table. (i) Ising sigma is a real field, so translation acts on it by +-1 (momentum 0 or pi); there is no incommensurate Ising order. (ii) In N the boson is fully floating and sigma V carries an incommensurate momentum, which is why C154 found k = 1 for P1 at every K. In T and K the commensurate V_c at momentum pi neutralizes sigma at pi, so sigma V_c is allowed. (iii) For the uniform odd sigma, the density d_x phi (sigma^m) is odd under Theta and under A in K_I; the product is even and translation invariant, a spin +-1 relevant coupling counted as in C68 and C154. (iv) In T the anomaly forces at least one boson with V_c at pi (section 2.2), and LSM forbids both bosons from carrying the anomaly (it is Z_2-valued), so in the two-boson case the second boson is anomaly free; its V_c may sit at 0 or pi depending on its translation and Theta action, which gives the last two rows.

### 4.2 The ordered side of P1

[derived] In T or K the ordered side of an Ising x boson point with sigma at pi has <sigma> != 0, which turns sigma V_c into a potential <sigma> V_c at zero momentum, of dimension x. For x < 2 it pins the boson. So a phase with spontaneous period-2 order and a gapless half-floating boson ("floating + SSB", pattern P1 of C154) exists only for x > 2. In N the corresponding argument of C154 section 5.4 used incommensurability instead and needed a second species; in T and K it is an operator statement.

### 4.3 Consequence

Proposition C. [derived, kinematic] In T and in K, every c > 1 candidate from the list above that occupies a locus of codimension <= 1 in the stratum contains a half-floating boson whose commensurate vertex has dimension x > 15/8 (Ising x boson with sigma at pi, which also needs x > 2 for its ordered neighbour), or two half-floating bosons with x1 + x2 > 2 (codimension 0). In XXZ language x = K: the first needs a strongly attractive liquid (K > 15/8, which for the zero-field parent means Delta < cos(11 pi/15) = -0.67), the second needs at least one attractive liquid (K > 1, Delta < 0). All other candidates (SU(2)_k, SU(2)_1 x Ising, Ising x boson with uniform order, and every candidate with x below these thresholds) need k >= 2, i.e. loci of dimension <= 3 in a 5-dimensional stratum, which must lie near strata where extra symmetry removes operators (parity, cyclic, XYZ, U(1)), all of which are controlled (section 3.5) or pending (parity, C152).

So the answer to "can an anomaly-matching c > 1 theory be reached in T, and where" is: symmetry and anomaly do not exclude it, but at codimension <= 1 it can only sit next to attractive half-floating liquids (x > 1 for two bosons, x > 15/8 for Ising x boson). SU(2)_2 is excluded by anomaly and SU(2)_3 by its operator count. In T the attractive liquids come from parents with Delta < 0 (the ferromagnetic side of the XY regime) and in particular from the easy-plane neighbourhood of the ferromagnetic Heisenberg point. The Typhon scans TFM and KFM are aimed there.

[heuristic] The strongest kinematic window, two attractive half-floating bosons, would need a second gapless species at zero magnetization (T) or at a commensurate transverse momentum (K_II). In T the bound-pair mechanism of contrib/d2_two_fluid (a pair sea next to a magnon sea, pinned by U(1) at finite field) has no field to work with: Theta fixes zero magnetization along every spiral axis. The remaining source would be a second minimum of a single-particle band, which Lemma 3 excludes at every ferromagnetic onset.

## 5. Numerics

### 5.1 Done locally (seconds per run; the local machine was loaded)

- [numerical] halffloat_test.py (section 2.2): S_mm(k) peaks at k = pi at the T and K anchors and at zero-field U(1), off pi at the N anchor and at U(1) in a field.
- [numerical] scan_tk.py pilots, 4 T and 4 K random points at rings 12/16 (scan_local/): every point has exactly one nontrivial uniform symmetry (Theta, respectively the reflection), confirming the normal forms. The T pilot gapless points read c_S = 0.99 and 1.00 with c_cas 1.02 and 1.03; the ordered T point has a twofold ground state (ndeg = 2, Theta broken); K points read 0 (gapped) or below 0.92.
- [numerical] scan_tk.py pilot of the ferromagnetic neighbourhood in T (TFM, 6 points, scan_local/TFM_999.jsonl). Parents (AF convention, Delta = s_m/|s_perp|) range from -1.16 (easy-axis ferromagnet) to -0.52. Easy-axis points read c_S 0.87 to 0.91 with crit 0.42 to 0.51 (gapped or slow, not critical by the 1/L test). The attractive point v = (-0.119, -0.629, -0.418, 0.155, 0.498, 0.380), parent Delta = -0.786 (K_parent = 2.35, inside the window K > 15/8 of Proposition C), read c_S = 1.48 and c_fit = 1.30 at rings 12/16 but c_S = 0.71 and c_fit = 1.16 at rings 14/18 (crit 0.11, ground-state momentum changing between sizes). This is the shell effect of an incommensurate liquid at small sizes (C155, d2_classification), not evidence for c > 1; it is the kind of point that stage B refines with DMRG.
- [proved, numerical identity] khg_check.py (section 1.3).

### 5.2 Typhon job 17597033 (pending at the end of this run)

Packed one-node job, partition typhon, 48 cores, 8 h limit, directory /scratch/lustre/jkudlerflam/c159/.
- Stage A (scan/*.jsonl). C155's ring estimator (rings 12, 16; c_S, c_cas, c_fit, crit, degeneracy) on uniform samples of T (16 x 1500 points) and K (16 x 1500), and on the neighbourhoods of the ferromagnetic Heisenberg point in T (TFM, 6 x 1000, radius 0.1 to 0.9) and K (KFM, 4 x 1000) and of the AF Heisenberg point in T (TAF, 6 x 1000, radius 0.2 to 1.2). Each record also carries the spiral-parent data (Delta, theta, harmonic fraction) and the L = 16 levels per momentum. Tests: (i) Lemma B: count of "symmetric Ising-like" points in T (0.3 < c_S < 0.7, crit < 0.15, unique ground state at both sizes), against the same count in K where Ising points are allowed; (ii) distribution of c over gapless points, split by the parent Delta (Proposition C predicts that any c > 1 at codimension <= 1 needs attractive parents); (iii) maxima of c_obj.
- Stage V (fes/*.json). Damped VUMPS FES (two-site cell, chi to 128) at the C68 anchors of T, K and N.
- Stage B (refine/*.json). For the top 8 points per mode with c_obj > 1.1 and crit < 0.15, the best 8 "symmetric Ising-like" points of T, TFM, K, KFM, and 4 random gapless controls per mode: rings 16/20 and open-chain DMRG at L = 32, 64, 128 (chi to 128), estimators c_prof, c_sm, c_two of C155's dmrg_c.py.
- Summary: python3 summarize_tk.py > summary.txt (run automatically at the end of the job).
- Interpretation rules fixed in advance (from the calibrations of C155 and contrib/d2_classification): ring (12, 16) values up to 1.2 to 1.4 occur at c = 1 near ferromagnetic edges and at incommensurate liquids; a c > 1 candidate needs c_two and the profile fits above 1.2 at L = 128 with truncation below 1e-7, and should then go to FES with chi >= 256.

## 6. Conditional theorem for T and K

Assumptions, as in C154 section 7: (A1) gapped, conformal or singular; (A2) the weak c-theorem near conformal bonds; (A3) c <= 1 on U(1), free locus and zero-field XYZ; (A4) c <= 1 on the parity class (C152, pending); (A5) emergent c <= 1 near singular points (C141 to C144, C154, C158; the XZ hyperplane C157 pending); (A5') c <= 1 on the cyclic stratum (C154). New input from this run: the controlled neighbourhoods of section 3.

(H_T), (H_K). [assumed] Every connected critical component of T (respectively K) has a closure point in the controlled set of section 3.5.

Theorem. [derived, conditional on A1 to A5' and H_T, H_K] Every conformal bond of T and of K has c <= 1. In T every conformal point that breaks no symmetry spontaneously has c = 1 exactly (Lemma B gives c >= 1), critical points inside a broken background can have c = 1/2, and nothing else occurs.

Proof. The closure lemma of C154 section 7 (c is constant on a critical component, does not increase towards a conformal frontier point, and is bounded by the emergent value at a singular frontier point) propagates c <= 1 from the controlled set into every component that touches it. QED.

What H_T and H_K contain. [derived, Proposition C] They fail only for a c > 1 island. By Proposition C, an island of codimension <= 1 needs attractive half-floating sectors (x > 1 for two bosons, x > 15/8 for Ising x boson). Islands of codimension >= 2 need two tunings in a 5-dimensional stratum and, by C154 section 5.2, contain lower-c sub-loci of P1/P3 type in their neighbourhood, which again require attractive sectors. So H_T and H_K can be replaced by the single statement

    (A6_TK) [assumed] no point of T or K hosts two independent gapless sectors one of which is an attractive half-floating boson (x > 1),

which is narrower than C154's general one-species assumption (A6), because Lemma A and Lemma B remove the generic (x arbitrary) two-boson and Ising x boson options that N allows.

## 7. Consistency checks

- Lemma 1 against C155's RA coordinates (Jx, Jy, Jz, Dz, fx, fy): identical. Random K points have exactly one nontrivial uniform symmetry (scan_local).
- Lemma A against known liquids: the zero-field XXZ chain with DM (phi shift pi/2, vartheta shift pi + alpha) in T; XXZ in a longitudinal field with transverse anisotropy (K_II type, vartheta commensurate, field-induced incommensurate magnetization, gapped by the A-even e^{2 i vartheta} unless K < 1/2, which matches the Chepiga-Mila picture of XYZ in a field where floating lobes need 1/4 < K < 1/2).
- Lemma B against the free locus cap T (only c = 1 or Theta-broken gapped, proved from the band), the XZ point (c = 1/2 but singular), C144's c = 1/2 line (inside a dimerized background), the cyclic stratum (c = 1 or dimerized/ordered), KHG (c = 1 liquids, all ordered phases break Theta), Peotta et al. at h_z = 0 (Ising lines between Theta-broken phases).
- Section 3.2 at D = 0 reproduces Baxter's XYZ result (c = 1 only when the two largest couplings coincide).
- Proposition C at N: dropping Lemma A restores C154's table (P1 at k = 1 for all K, P3 at k = 0 for all K).
- MPS form of LSM (section 2.1): odd-cell VUMPS fails at the T anchor; two-site cells are used.

## 8. What is not done

- The scans and refinements of section 5.2 have not run yet (job pending). Until they report, the numerical statements about T and K beyond the anchors and the literature families rest on C155's stratum-confined searches (TR, RA, also pending) and on the earlier 103,000-point survey.
- Lemma B for the D- and E-series minimal models uses the assumed anomaly-freeness of their invertible symmetries.
- The candidate list of section 4 is not exhaustive (irreducible c > 1 theories with few relevant operators, for instance N = 1 superconformal minimal models with c between 1 and 3/2, are not counted). For these the operator count needs their embedding, as C68 section 5 already noted.
- H_T and H_K (equivalently A6_TK) remain assumptions.
