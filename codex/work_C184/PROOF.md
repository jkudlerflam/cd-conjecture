This run was done by Claude Opus 5.5 (claude-opus-5-5), agent C184, not by Codex. Astra was not used (out of credits).

# C184 derivations: multicritical points with c > 1 in B_2

Standard: theoretical physics. Tags: [proved], [derived] (physics standard, inputs named), [numerical] (file named), [heuristic], [assumed], [literature] (arXiv number opened in this run or quoted from a RESULT.md that opened it, which is said).

Conventions. Bond h = J_x XX + J_y YY + D(XY − YX) + J_z ZZ + (h_z/2)(Z⊗1 + 1⊗Z) in the parity class; in C152's fermion units w = t e^{iφ} = J_x + J_y + 2iD, Δ = J_x − J_y, V = 4J_z, μ = V + 2h_z. Compact boson normalized as in C165 and C180: for the U(1) liquid of spin-1/2, e^{2iφ} (staggered S^z, dimer) has dimension K and e^{iθ} has dimension 1/(4K), with K = 1 at the XX point and K = 1/2 at the Heisenberg point.

## 1. Framework: what happens where two transition surfaces cross

Let S_1, S_2 be transition surfaces of B_2 (or of a stratum), with infrared CFTs C_1, C_2 and tuning operators O_1, O_2 (the relevant, symmetry-allowed, translation-invariant scalars that the surface tunes to zero). At a transversal crossing there are three possibilities.

(i) Same sector. Both transitions are transitions of one and the same low-energy field. The crossing is then a multicritical point of a single theory, with c no larger than the c of that field. Examples in B_2:
- Ising × Ising crossings on the Ashkin–Teller plane of the parity class (C152 PROOF 8.1 to 8.3). The two Ising surfaces are the two Majorana masses m_± ∝ μ ± 2t cos φ, the coupling ε_1ε_2 has dimension 2 and is exactly marginal, and the crossing is the AT line with c = 1 [derived, C152].
- The zero-field XYZ point J_x = J_y = J_z where the three Gaussian lines J_a = J_b meet. It is SU(2)_1 with c = 1 [exact, Baxter].
- The KT–PT meeting point at Baxter's J_x = −J_z line, eight-vertex class with c = 1 (C152 PROOF line 84) [literature].
- Ising–chiral junctions in the parity class, one-sided cubic Lifshitz points (C152 result 4), and the Ising–PT junctions in N (C154 5.4). These are z ≠ 1 and not conformal.

(ii) Independent sectors, decoupled fixed point. C_1 × C_2 is a fixed point with c = c_1 + c_2. Its codimension in the bond space is

    k = k_1 + k_2 + #{ mixed scalars Φ_1Φ_2 that are relevant (Δ_1 + Δ_2 < 2), invariant under every exact symmetry and under the infrared image of translation, and not redundant }.      (1.1)

If some mixed scalar is relevant and allowed, the generic crossing (that coupling not tuned) flows away from C_1 × C_2. By the c-theorem it goes to a fixed point with c < c_1 + c_2, or to a gapped or first-order crossing. The decoupled point then survives only on a locus of higher codimension, given by (1.1).

(iii) Independent sectors, coupled fixed point. A relevant mixed coupling drives the crossing to a new fixed point with c < c_1 + c_2 (the standard example is the cubic or bicritical fixed point of two order parameters). For c > 1 at the crossing this still needs c_1 + c_2 > 1 at the decoupled point, so the same independent sectors.

So the multicritical route to c > 1 changes the codimension count, and nothing else. It needs two independent soft sectors at one bond, which is the content of C154's (A6) and C168's (R_N). Section 4 makes this a statement about Majorana species.

## 2. Does the closure lemma exclude multicritical c > 1? No.

[derived; this agrees with C174 R1, which I re-checked]
C154 PROOF 7, closure lemma: at a conformal frontier point h′ of a critical component Σ, c(h′) ≥ c(Σ). A multicritical point h* with c(h*) > 1 is a frontier point of the components of its lower-c neighbours (the Ising surfaces with c = 1/2, the Gaussian surfaces or liquids with c = 1). The lemma asks c(h*) ≥ 1/2 and c(h*) ≥ 1 there, which holds. So the lemma is consistent with a multicritical c > 1 point and gives no exclusion. Bounds move from high c to low c, never back.

The multicritical locus Σ* = {bonds flowing to C_1 × C_2} is itself a component, of codimension k. By (A2) its closure contains no conformal point of lower c (C174 R1), so no controlled conformal point of any stratum. C154's (H) then holds for Σ* only if Σ* accumulates on a controlled singular point, which is (H*). Multicritical c > 1 is therefore excluded only by (H*)/(A6) itself, and not by any topological step. A point of Σ* that lies inside a symmetric stratum is not covered by that stratum's generic-point result either, which is C174 R4. The inventory below is meant to fill that gap stratum by stratum.

One more remark on the lemma's inputs. (A2) (C141's weak c-theorem) is used at h* only in the form "bonds near h* have c ≤ c(h*)". That is harmless for this question. It does not say that h* has c ≤ 1.

## 3. Majorana-species counting at weak coupling

Lemma M. [proved for quadratic chains, C152 Prop. F and Lemma P2, C174 R8, C169 §5; derived for weak interactions] A nearest-neighbour qubit chain that is quadratic after Jordan–Wigner (two Majoranas per site, any period p, nearest-neighbour Majorana bilinears) has det H_BdG(K) a polynomial of degree at most 2 in cos K (single bond across a cell boundary), hence at most four Majorana–Weyl zero modes. Its gapless theory is at most one Dirac fermion or two Majoranas, c ≤ 1, at every point and at every crossing of its transition surfaces. Adding interactions of strength g small compared with the gap of all other modes keeps the soft content fixed, since four-fermion operators built from at most four Majorana–Weyl fields have dimension at least 2 (marginal or irrelevant). So every crossing in the perturbative neighbourhood of the free locus, of the XZ point in directions where B stays slow and decoupled (C143, C157), of the AT plane and of the U(1) class has c ≤ 1.

What a c > 1 crossing needs, in Majorana–Weyl components: Ising × Gaussian needs 3 Majoranas per chirality (6 components), TCI × Ising needs a strongly interacting Majorana plus another (c = 6/5 cannot come from free Majoranas at all), Ising² × boson needs 8 components, Potts × boson and Z_5 parafermions are not free-fermion theories. Lemma M gives 4. So every c > 1 crossing needs a sector that is not in the quadratic theory of the chain: a bound state, a blockade-generated band, or a strongly coupled slow sector. In B_2 these are known to exist only at
- (E1) the XZ neighbourhood, where the B Majoranas are a flat band (C143, C144, C157, C158);
- (E2) the Fibonacci point and the cusp, where the blockade generates a second-order band with two valleys (C173 to C180);
- (E3) strong attraction, where bound pairs and clusters form (C152 §7, C175 for the cusp);
- (E4) the neighbourhood of the ferromagnetic Heisenberg point (semiclassical, one soft mode, C160, refereed).

## 4. The inventory

Each entry gives the theory, what the crossing flows to, its codimension, whether a nearest-neighbour bond family has both transitions, and the status.

### 4.1 Ising × Gaussian, c = 3/2

Mixed scalars at the decoupled point, with σ, ε the Ising fields (1/8, 1) and the boson at Luttinger parameter K.
- ε (∂φ)²: dimension 3, irrelevant. ε J (current): spin 1, Lorentz-violating, dimension 2, shifts velocities only.
- ε cos(nφ) or ε cos(mθ) for the cosine O_g that the Gaussian surface tunes: dimension 1 + x_g. It mixes the two tuning parameters and is harmless (redundant with the tunings).
- σ V, with V a vertex of the boson carrying the same quantum numbers as σ: dimension 1/8 + x_V. This decides stability.

Three cases, by where σ lives.

(a) σ at momentum π, boson commensurate (a Gaussian surface at half filling, translation acting as a shift of φ). [derived] For spin 1/2 the boson always has momentum-π vertices, the staggered density and dimer e^{±2iφ} with x = K, and, if σ is odd under the same on-site Z_2 as the in-plane order, also e^{±iθ} at x = 1/(4K). σ e^{2iφ} is relevant for K < 15/8 and σ e^{iθ} for every K > 1/30. A Gaussian surface of a spin-1/2 chain driven by cos 2θ (in-plane anisotropy) needs K > 1/2, and one driven by cos 4φ needs K < 1/2. In both windows a mixed σ·V coupling is relevant and allowed, with one exception. If K > 15/8 and σ is even under the in-plane π rotation Π, then σe^{±iθ} is Π-odd and forbidden, σe^{2iφ} is irrelevant, and the charge-2 momentum-π vertex is e^{2iθ}cos 2φ of dimension 1/K + K > 2. So the decoupled point there has k = 2 (ε and the cos 2θ coefficient), a kinematic window. K > 15/8 is −1 < Δ_eff < −0.67, the same threshold as case (a) of C167's residual (A_W1), and the independent σ it needs is exactly what (A_W1)(a) assumes away. So this window adds no new route, and it is where a targeted search would go. Outside it the decoupled point has codimension at least 3, and the generic crossing locks: the σ·V term orders σ and V together, gapping both, or makes the crossing first order. This is case (iii) with c < 3/2. There is a stronger statement that does not need the counting. At spin 1/2 every momentum-π order (Néel along x, y, z, dimer) is a component of the one SU(2)_1 matrix field g, to all orders in weak anisotropy (C165, refereed by C170, which found the statement stronger than claimed). So near any SU(2)_1 or U(1) liquid there is no independent momentum-π Ising sector to cross with, and crossings of momentum-π transitions are case (i) with c ≤ 1. The XYZ-in-a-field crossings fall here: the transition between x-ordered and y-ordered phases in an x field is an AT/Gaussian transition of the in-plane angle θ, not an Ising × boson point.

(b) σ at momentum π, boson floating (incommensurate). Every σV carries momentum and is forbidden, so the decoupled point is stable with k = 1. This is C154's pattern P1 and the only Ising × boson crossing that is generic on a hypersurface. [derived] It needs an independent period-2 sector inside a floating liquid. C165 and C170 exclude it at weak harmonics, and C172's locking lemma pins it onto the codimension-2 locus near T. At order-one harmonic amplitude it is (R_N)(i),(ii) of C168, open.

(c) σ uniform (Ising at momentum 0, no protecting symmetry). σ itself must be tuned, k ≥ 2, and σJ (dimension 9/8, spin 1) is allowed. This is P2, which needs a first-order surface between two coexisting liquids ending on the Ising × boson line, with an order field independent of the boson's density (C154 5.2).

(d) The known lattice mechanism for Ising on top of a Luttinger liquid is pair unbinding: a pair liquid whose single-particle gap closes, which gives c = 3/2 on the transition between an atomic and a pair liquid. Ruhman and Altman, arXiv:1705.00635 [literature, abstract opened], find it with pair hopping, which in qubit language is a four-site term. Chepiga and collaborators, arXiv:2211.15598 [literature, search snippet; that is the interacting Majorana chain of Rahmani et al., with three-site terms], see an Ising line on top of a floating phase. The spin-1 bilinear–biquadratic chain in a field, arXiv:1012.4518 [literature, search snippet], has the same pair-unbinding c = 3/2 at d = 3. In B_2 a pair liquid needs bound pairs that repel each other. In the U(1) class there are none in the ground state [exact, Bethe ansatz: real roots only for −1 < Δ, ferromagnet below]. With a nearest-neighbour attraction, the same bond that binds a pair binds clusters ((n − 1)V for an n-cluster, C152 §7) [heuristic for Δ ≠ 0], and at the Fibonacci cusp trimers bind (C175) [numerical]. So the pair-unbinding route meets cluster binding (phase separation) before the pair liquid forms. Whether some bond in the parity class or N has pair–pair repulsion with pair binding is not excluded. It needs the pair–pair repulsion to come from a longer-range effective interaction, as the blockade gives at the cusp, which is C181/C182's territory.

### 4.2 Ising × Ising × boson, c = 2

[derived] Mixed scalars: ε_1ε_2 (2, marginal, AT), σ_1σ_2 (1/4) if both orders have the same quantum numbers, plus the case 4.1 terms for each Ising. With translation acting as Z_2 on both σ and the boson floating, k = 3 (C154 table). It needs eight Majorana–Weyl components, twice what Lemma M allows. The two Ising factors would have to be the two Majorana species of the AT plane, and the boson a third, emergent sector. On the AT plane both Majorana species are already critical by translation (C152 PROOF 8.3(a)), so the plane is the natural host. But there is no emergent boson on that plane outside the XZ window (E1), where the B Majoranas are themselves the slow sector, so there is nothing left to form it. Status: no B_2 family has both ingredients.

### 4.3 Tricritical Ising × Ising, c = 6/5

[derived] At TCI × Ising the mixed scalars are ε_T ε_I (dimension 1/5 + 1 = 6/5, relevant, invariant under every symmetry since both factors are Z_2-even, so always allowed) and σ_T σ_I (3/40 + 1/8 = 1/5, relevant; allowed unless a symmetry tells the two Z_2's apart). With the TCI tunings ε_T (1/5), ε′_T (6/5) and ε_I (1), the decoupled point needs k ≥ 4, and k = 5 if σ_Tσ_I is allowed. In the 4-dimensional parity class that is isolated points at best.

The dynamical obstruction in the parity class. Its two Majorana species A and B (C152 PROOF 8.1) have, at the quadratic level, independent masses m_± and velocities t_A, t_B. The only interspecies coupling a nearest-neighbour bond supplies is V ε^A ε^B (the four-Majorana term V n_j n_{j+1}). A TCI point of species A needs an A self-interaction of the O'Brien–Fendley type, ε^A_j ε^A_{j+1}, which is a three-site term (C169 §2.2), comparable to t_A. In B_2 such a term exists only as a virtual-B effect, of order V²/m_B, which needs B gapped, while TCI_A × Ising_B needs m_B = 0. With m_B = 0 the same V is the relevant ε_Tε_I coupling, so the crossing flows away from TCI × Ising. This holds perturbatively in V/t [derived] and is a heuristic at V of order t. The B_2 realization of TCI_A × Ising_B would need V large enough to drive A tricritical directly (as in the interacting Kitaev chain at strong attraction), with B still critical. At V = −1.5t, Δ = 0.5t, φ = 0 and 0.8 the transitions seen in the μ scan of C152's job 17597010 are first-order or gapped (section 6 below) [numerical].

### 4.4 Three-state Potts × boson, c = 9/5

[derived, C169 §9.2, C168, C171] Codimension 2 on T and K_I and at least 3 in N (the chiral Huse–Fisher term is never symmetry-forbidden). At spin 1/2 a period-3 order is either the host's locked spiral, which gaps the host, or a φ-side period-3 density wave of dimension 9K > 9/4, which no nearest-neighbour liquid makes relevant since K > 1/4. So the Potts sector is a second species, under (A_W′) or (R_N). Mixed couplings: σ_P carries momentum 2π/3 and couples only to boson vertices at ±2π/3, absent for a floating boson. So the decoupled point is stable if it exists. The obstruction is the existence of the second species, not the coupling.

### 4.5 Z_5 parafermions, c = 8/7

[derived, C169 §9.1, C168 Lemma P] B_2 has no on-site Z_m (m ≥ 3) beyond U(1), so the Z_5 is translation acting on a period-5 tongue. The chiral term is never forbidden, codimension 3. For m ≥ 5 the clock vertex and any disordering vertex cannot both be relevant in a one-boson description (Δ(O^m)Δ(W) ≥ m²/4), so controlled tongues melt through a c = 1 floating phase. Open only at order-one harmonic amplitude, (R_N)(iv). Not a crossing of two surfaces in the sense of section 1, it is a single-species multicritical point.

### 4.6 SU(2)_2 and other WZW points

[derived] SU(2)_2 has c = 3/2 and is three Majoranas.
- With translation acting trivially on the SU(2) currents and no symmetry, all ten scalars tr g and Φ^{ab} (dimensions 3/8 and 1) are allowed, k = 10 > 8 (C154 table): impossible in N.
- With translation acting through the maximal torus of SO(4), only Φ^{zz} survives, and SU(2)_2 is the special-radius point of Ising × floating boson. That is case 4.1(b).
- An exact SU(2) exists only on the twisted Heisenberg curve (C172), where the chain is the Heisenberg chain, SU(2)_1.
- Anomaly. On T (Kramers Θ, Θ² = −1 per site) and on the cyclic stratum (projective Klein group) the Ogata–Tasaki anomaly must be matched (C154 §1.5). For an SO(3)-symmetric spin-1/2 chain the LSM anomaly is matched by SU(2)_k with k odd [literature, standard; not opened here]. With only Θ or the Klein group the matching condition on k is weaker and I have not derived it, so this is not used.
- SU(2)_3 (c = 9/5): on the cyclic stratum it would be a curve fixed by the A_4-singlet thermal operator (dimension 4/5), and C154 §3.5's 20-point scan found none (C174 R4: grid too coarse to exclude it between points). Elsewhere it has more relevant scalars than SU(2)_2.

Status: SU(2)_2 reduces to 4.1(b). SU(2)_3 on the cyclic stratum is not excluded between grid points. That is a two-parameter search with a clear signature (c = 9/5, two-site entanglement), which nobody has run finely.

### 4.7 The c = 2 point where two Gaussian surfaces meet

[derived] Two cases.
- If both Gaussian surfaces belong to one boson (two different cosines of the same φ, θ), the crossing is c = 1 with two cosines tuned. Example: the three Gaussian lines of zero-field XYZ meet at the Heisenberg point (SU(2)_1).
- If they belong to two bosons, the crossing is P3 with both bosons at commensurate Gaussian points. Mixed vertices cos(m_1φ_1 ± m_2φ_2), and the θ analogues, carry momentum m_1 k_1 ± m_2 k_2, and at a commensurate crossing some of them are invariant and relevant (C180's operator table: cos 2√2φ₋ and the two-pair vertex at the cusp). So k ≥ 3.

Two bosons in B_2 exist only at the cusp (E2), which is C181 to C183's and C186's subject.

### 4.8 Lifshitz-type points

[derived] Every Ising surface in the parity class that meets the chiral liquid G1 does so at a one-sided cubic Lifshitz point (C152 result 4, margin c_3 = t² cos² φ/(4Δ)). Ising surfaces in N end on Pokrovsky–Talapov surfaces (C154 5.4). These have z ≠ 1 and carry no central charge. They are excluded as conformal points, as the spec asks. The protection margin c_3 vanishes as cos φ → 0, so near the AT plane and the XZ window the junction is not protected by the cubic term. That is where the strongly coupled regions E1 and E3 come in, and the target of section 5.

### 4.9 Irreducible c > 1 theories beyond the list (C174 R5)

[heuristic] N = 1 superconformal minimal models with c > 1 (81/70 first) are higher Z_2-multicritical points with supersymmetry. Already the first of them (TCI, c = 7/10) needs a three-site term in the Majorana chain (O'Brien–Fendley). The next ones need more relevant operators tuned and longer-range self-interaction of one Majorana species, which B_2 supplies only virtually (4.3). W_3 minimal models with c > 1 (6/5, 10/7, ...) are Z_3-multicritical and need a period-3 order, which 4.4 says needs a second species at spin 1/2. Neither has a B_2 family with the ingredients.

### 4.10 Summary table

| candidate | c | k (lowest, kinematic) | what the generic crossing flows to | ingredient missing in B_2 | status |
|---|---|---|---|---|---|
| Ising × Gaussian, commensurate | 3/2 | ≥ 3 | locks (σ·V relevant) or AT/Gaussian if same field | independent momentum-π Ising sector (one SU(2)_1 field at spin 1/2) | derived c ≤ 1 near liquids |
| Ising × floating boson (P1) | 3/2 | 1 | decoupled, stable | independent period-2 sector in a liquid | open at order-one harmonics (R_N i, ii) |
| Ising (uniform) × boson (P2) | 3/2 | 2 | decoupled | first-order surface between liquids with independent order | open, (R_N) |
| pair unbinding (Ising on pair liquid) | 3/2 | 1 inside a pair liquid | decoupled | a pair liquid; NN attraction binds clusters | heuristic exclusion, open at the cusp |
| Ising² × boson | 2 | 3 | AT × boson | 8 Majorana–Weyl modes, Lemma M gives 4 | no family |
| TCI × Ising | 6/5 | 4 to 5 | flows away (ε_Tε_I relevant) | A self-interaction needs B gapped | derived perturbatively, numerics first order at V = −1.5t |
| Potts × boson | 9/5 | 2 (T, K_I), 3 (N) | decoupled if floating | period-3 second species | C169, C171, open under (R_N) |
| Z_5 parafermion | 8/7 | 3 | — | order-one clock harmonics | C168, C169, open under (R_N iv) |
| SU(2)_2 | 3/2 | 10 trivially, 1 as P1 | — | reduces to P1 | as P1 |
| SU(2)_3 (cyclic) | 9/5 | 1 in the stratum | — | none kinematically | coarse scan only (C154 §3.5) |
| two Gaussian surfaces, two bosons | 2 | ≥ 3 commensurate, 0 floating | — | a second boson (cusp only) | C181 to C186 |
| two Gaussian surfaces, one boson | 1 | — | — | — | c = 1 (Heisenberg point) |
| Lifshitz | — | — | z ≠ 1 | — | not conformal |

## 5. The crossing chosen for numerics

Choice. The region of the parity class at φ = 1, h_z = 0, V ≈ −2t, Δ ≈ 0.3t to 0.5t. Reasons, in order.
- It is the one uninventoried junction of a symmetric stratum named by C174 R4 that is not near the cusp (which C181 to C183 and C186 cover) or XZ (closed numerically in section 6).
- Two Ising surfaces (the DM-split x-FM/z-FM transition of Peotta et al., arXiv:1403.4568 [literature, abstract opened]) meet the chiral Luttinger region and its first-order boundary there. So three surfaces come together.
- It sits at strong attraction (E3), outside Lemma M's perturbative control, where bound clusters exist. That makes it the place in the parity class where a third, emergent species is least excluded.
- It is inside the kinematic window of 4.1(a). On the U(1) line Δ = 0 of this slice the bond is XXZ + DM, which gauges to XXZ with J′ = t/2 and Δ_eff = (V/4)/(t/2) = V/2. So the Luttinger region −2 < V < −1 has −1 < Δ_eff < −1/2, and K > 15/8 for V < −1.34. There a Π-even momentum-π Ising sector would cross a Gaussian surface at codimension 2, the lowest commensurate value in the inventory.
- The cubic Lifshitz margin c_3 = t² cos²φ/(4Δ) is of order one at φ = 1 (cos²φ = 0.29), so the free-fermion picture predicts a Lifshitz junction (P3 of the registration). A conformal c > 1 reading would therefore be a clean interaction effect.

Predictions are registered in RESULT.md. Results are in section 7.

## 6. Collection of C152's job 17597010 (existing data, read here for the first time)

[numerical, c152_j2_copy/, summary c152j2_summary.out from summ_c152j2.py; S̄ = mean of the four central cuts, c_pair(L_1, L_2) = 6[S̄(L_2) − S̄(L_1)]/ln(L_2/L_1)] The job completed on 2026-11-08 in 2 h 16 min, and its output had not been collected. Truncation is below 3e-9 throughout.
- XZ window (Set 1). The largest C152 reading, φ = 1.45, s = 1.1, V = 0.4, h_z = 0, gives c_pair = 1.46, 1.24, 1.13 for (32, 64), (64, 128), (128, 192), and the window fit falls 2.15 → 1.91 → 1.47 → 1.31 over L = 32 to 192. The other repulsive XZ-window points read c_pair(64, 128) between 1.04 and 1.16, all below their (32, 64) values. At φ = π/2 every point reads 0.86 to 0.99. The attractive points (V = −0.4, −0.8) are gapped or ordered at L = 128 (c_pair 0.32 to 0.44, or entropy jumps of ln 2 from cat states). So the XZ window falls towards c = 1, as C152 and C158 expected.
- Set 2 (fine h_z scans at φ = 1.3, V = −0.4): the largest c_pair is 1.07, at (64, 128).
- Set 3 (strong attraction, V = −1.5, Δ = 0.5, φ = 0 and 0.8, μ from −2.5 to 0.5 in steps of 0.2): every point is gapped at L = 48 and 64 (c_pair ≤ 0.16 apart from one level-crossing jump). The transition is crossed between grid points, which is consistent with a first-order transition at this attraction. No tricritical Ising point was located, and none was resolvable on this grid.
- Set 4 (the Peotta meeting region, φ = 1, h_z = 0, Δ = 0.25, 0.32, 0.4, V = −1.6 to −2.3, L = 48 and 96): the largest c_pair(48, 96) is 1.15 (Δ = 0.4, V = −2.0), with neighbours 0.94 and 0.85. At Δ = 0.25 and 0.32 the larger |V| points are gapped, some with cat-state entropy drops. Nothing reaches 1.2.

## 6b. Junction inventory by stratum (C174 R4), as far as this run takes it

[derived from the cited runs plus sections 1 to 4; numerical where marked]
- U(1) class. The only junctions are the KT point (SU(2)_1), saturation and the empty band (z = 2), and the ferromagnetic point (z = 2). All are exact, c ≤ 1.
- Free locus. Every junction has c ≤ 1 by Lemma M [proved].
- Parity class (C152 table). Ising–Ising crossings: the AT plane, c = 1. Ising–G1 junctions: cubic Lifshitz, not conformal, with margin c_3 = t² cos²φ/(4Δ). Eight-vertex KT–PT points: c = 1. Possible TCI points at attraction: c = 7/10 if they exist. Set 3 at V = −1.5t found the μ-driven transitions first-order on a 0.2 grid [numerical, §6]. TCI × Ising needs k ≥ 4 and flows away (4.3). The XZ window is §6 [numerical]. The Peotta region is §7.
- Staggered U(1) class (2-dimensional). One charge, one boson. Junctions inside it are those of XXZ in a staggered field (Néel plateau edges, PT), and c ≤ 1 by the U(1) argument of C153, except its large-b corner (C153 heuristic, job 17596442 not collected by me).
- Cyclic stratum (2-dimensional). The SU(2)_1 basin and the ordered region meet along a boundary that C154 §3.5 did not resolve. A c > 1 point there would be SU(2)_3 (c = 9/5) or a first-order boundary. Kinematically SU(2)_3 is codimension 1 in the stratum (4.6), so it would be a curve. This is the one place in the inventory where a c > 1 multicritical locus is kinematically cheap and no dynamical argument excludes it. Its natural test is a fine scan of the SU(2)_1/ordered boundary on the cyclic stratum, which I did not run (section 8).
- T and K. Covered by (A_W′) and C171's routes. The Potts × boson and P2 crossings of 4.1 and 4.4 live there.

## 6c. Collection of C154's job 17594765 (the 53-point cyclic grid, never collected)

[numerical, c154_cyc_copy/cyc_grid.out, from C154's cyc_grid.py; momentum-resolved rings L = 12 and 18. h = XZ + b YX + g ZY. x_t is L times the lowest triplet gap, r_s and r_c are the singlet and current gaps over the triplet gap. SU(2)_1 with log corrections has r_s ≳ 1 and r_c ≈ 2.] The job completed on 2026-11-08 in 13 minutes.
- SU(2)_1 region. Every point with g/b ≥ −0.3 has r_c between 1.83 and 2.29, a triplet degenerate to 1e-11 (exact A_4), and c_S between 0.99 and 1.04. c_cas lies between 0.98 and 1.08, except near the XZ corner (b ≤ 0.25), where it reaches 1.97 at (1, 0.1, −0.03). There v is small and the slow sector of C154 §3.6 dominates, and c_S stays at 0.96 to 1.04.
- The boundary. Along g/b = −0.6 the triplet gap collapses as b grows: x_t(18) = 2.03 (b = 0.4), 1.18 (0.55), 0.73 (0.7), 0.107 (0.85), and the point is ordered at b = 1. At (1, 0.85, −0.51), r_s = 9.7 and c_cas = 0.69. All points with g/b = −0.9 are ordered (gap1 = 0, ground-state momentum shifted, r_c ≈ 0.39). A collapsing triplet (the spin-1 level) with a falling c_cas is the signature of a level crossing into the ferromagnet-like ordered state. That is a first-order or z = 2 boundary, as at the XXZ ferromagnetic point. An SU(2)_3 point would instead show r_c near 10/3 (the j = 1/2 primary has Δ = 3/10, against 1/2 at level 1, and the current has Δ = 1), and c_cas rising toward 9/5. Along the approach r_c does rise, 2.18 → 2.44 → 2.88 → 8.1 (b = 0.4 to 0.85), but c_cas falls at the same time (0.93, 0.82, 0.75, 0.69). The rise comes from the triplet collapsing, not from a new tower. No grid point on the boundary shows the SU(2)_3 pattern.
- Reading [numerical, coarse]. No grid point shows c > 1 away from the XZ corner, and the SU(2)_1/ordered boundary looks like a level crossing. This answers the uncollected-job part of C174 R4 for the cyclic stratum. The grid spacing along the boundary is 0.15 in b, so a narrow SU(2)_3 segment is still not excluded. I did not run the finer scan.

## 7. The Peotta region: ring results (the DMRG stage is in Typhon job 17604462)

[numerical, ring_map.py, ring_phi1.jsonl, ring_phi1_table.out (L = 12, 16 on Δ = 0 to 0.7, V = −2.8 to −1.0), ring20.jsonl and ring_cS.out (L = 20 at eight points). Ring c_S(L_1, L_2) = 3[S(L_2) − S(L_1)]/ln(L_2/L_1) from the half-ring entropy of the global ground state.]
- Phases at φ = 1, h_z = 0.
  - z-FM: the ground state is doubly degenerate inside one parity sector, and L·(parity gap) grows linearly. This holds at V ≤ −2 for Δ = 0 and moves to more negative V as Δ grows.
  - A Luttinger region at small Δ for −2 < V < −1. The LL control (Δ = 0, V = −1.6) reads c_S = 1.04 and 0.97, with L·gap ≈ 0.45 to 0.49 flat.
  - An x-ordered region at Δ ≥ 0.5, V ≥ −1.4 (the parity gap closes faster than 1/L).
  - A band along |J_x| ≈ |J_z|, that is V ≈ −2(cos 1 + Δ), between z-FM and x-order. This is where Peotta et al.'s DM-split Ising pair runs.
- On the band (Δ = 0.5, V = −2.2 and −2.0; Δ = 0.6, V = −2.4 and −2.2), c_S(16, 20) = 0.98, 0.96, 0.99, 0.98, and c_S(12, 16) = 0.98 to 0.99. L·g_1 is close to flat (2.65, 2.51, 2.40 at Δ = 0.5, V = −2.2). One c = 1 boson at these sizes, not two. The DM splitting of the band into two Ising lines with a narrow disordered strip is not resolved at L ≤ 20. That is expected if the split is small, since the band is the D = 0 Gaussian line of Baxter's XYZ deformed by DM (C152 PROOF 8.3(b)).
- Toward the meeting point with the LL and z-FM (Δ = 0.3, V = −1.8 and −1.6), the low gaps fall faster than 1/L (L·g_1 = 0.59, 0.56, 0.30). The ring entropy jumps by about 0.2 between L = 12 and 16 at V = −1.6 and then stays flat (c_S = 2.03, then 0.01). That is a level crossing, the first-order LL–FM boundary of Peotta et al., not a c > 1 fixed point.
- Reading [numerical, small rings]. Nothing in the slice reads above c = 1 on rings up to L = 20. The band is c ≈ 1 and the meeting point is a level crossing. This matches P1 and P3 of the registration at ring level. The registered tests are open-chain c_pair at L = 48, 96, 192, which are in the pending job.
- A one-point local preview at L = 32, 64 (local_out/, stopped early to respect the 20-minute local limit) gave c_pair(32, 64) = 1.08 at Δ = 0.5, V = −2.4, inside the band. That is a crossover-size reading and not one of the registered pairs.

## 8. What is left open

1. The registered DMRG test (job 17604462): c_pair at (48, 96) and (96, 192) across the Peotta slice, the φ = 1.2 slice, and the qubit XX and TFI controls.
2. The kinematic window of 4.1(a) (K > 15/8, Π-even momentum-π Ising sector, k = 2) everywhere in the parity class and in T/K_I. That is (A_W1)(a), and a targeted search would locate an independent Π-even period-2 order next to a strongly attractive liquid.
3. A fine scan of the SU(2)_1/ordered boundary on the cyclic stratum (4.6, 6c). That is the only place where c > 1 (SU(2)_3) is codimension 1 inside a stratum and no dynamical argument excludes it, although the coarse data look like a level crossing.
4. The TCI × Ising obstruction (4.3) at V of order t is heuristic. A located TCI point on the attractive Ising surface of the parity class (none found on C152's 0.2 grid at V = −1.5t) would let it be tested directly.
5. Order-one harmonic amplitudes (P1, P2, P4, Potts × boson) remain under (R_N), and the cusp (E2) is with C181 to C183 and C186.
