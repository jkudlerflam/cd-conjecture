# C152 PROOF: the parity-symmetric qubit class

This run was done by Claude Opus agents (2026-11-07), not by Codex. It was interrupted and resumed by a second agent, which added sections 8 to 11. Tags: [proved], [derived] (physics standard, RG, bosonization, controlled perturbation theory), [literature], [numerical], [heuristic], [assumed]. Scripts are in this folder; numerical outputs are named where used.

## 1. The class, its symmetries and its fermion form

### 1.1 Reduction

A nearest-neighbour qubit bond h = σ^T J σ' + B·(σ + σ')/2 commutes with Π = ∏_j σ^z_j iff every term is even under σ^x, σ^y → −σ^x, −σ^y. So J_xz = J_zx = J_yz = J_zy = 0 and B ∥ z. The xy block of J is a real 2 × 2 matrix; its antisymmetric part D(XY − YX) is invariant under a common rotation about z, and that rotation diagonalizes the symmetric part. [proved] Every bond of the class is, after a common z rotation and up to a telescoping term,

    h = J_x XX + J_y YY + D (XY − YX) + J_z ZZ + h_z (Z ⊗ 1 + 1 ⊗ Z)/2.          (1.1)

Five couplings modulo scale: a 4-dimensional stratum of B_2. Any bond invariant under a π rotation about some axis n is in the class after the rotation n → z. In particular the XZ (compass) bond of C143/C144 is in the class with the parity axis along y. Rotating y → z (rotation by −π/2 about x) maps X⊗Z to −X⊗Y = −(XY + YX)/2 − (XY − YX)/2, and after a 45° z-rotation this is J_x = −J_y = ±1/2, D = ±1/2, J_z = h_z = 0, i.e. |J_x − J_y| = 2|D| with J_x + J_y = h_z = 0 (the flat-band point of the free locus, section 2.1).

Symmetries of the parametrization [proved, by direct substitution]: (a) a π/2 rotation about z swaps J_x ↔ J_y and fixes the rest; (b) a π rotation about z on every second site flips (J_x, J_y, D) → −(J_x, J_y, D); (c) a π rotation about x on every site maps (D, h_z) → (−D, −h_z); (d) site inversion maps D → −D. So one may take D ≥ 0, h_z ≥ 0, J_x ≥ |J_y|.

### 1.2 Jordan-Wigner form

With Z_j = 1 − 2n_j and c_j = (∏_{l<j} Z_l) σ^+_j one finds σ^+_j σ^-_{j+1} = c†_{j+1} c_j and σ^+_j σ^+_{j+1} = −c_j c_{j+1}, hence [proved]

    H = Σ_j [ w c†_{j+1} c_j + h.c. − Δ (c_j c_{j+1} + h.c.) + V n_j n_{j+1} − μ n_j ] + const,     (1.2)
    w = J_x + J_y + 2iD ≡ t e^{iφ},   Δ = J_x − J_y,   V = 4 J_z,   μ = 4 J_z + 2 h_z.

This is the interacting Kitaev chain with a hopping phase φ (the inversion-breaking DM term). Half filling and particle-hole symmetry sit at μ = V (h_z = 0). In the notation of Chepiga and Mila (PRB 107, L081106) t_CM = −(J_x + J_y), λ = −Δ (sign removable by c → ic), V = 4J_z, B = −h_z. [numerical, check_free.py] At V = 0 the open-chain ED spectrum coincides with the BdG modes of (1.2) to 1e-5, once the edge sites carry the full field (the bond-split form (1.1) gives them h_z/2; dmrg.py and dmrg_pt.py restore it).

The hopping phase is a pure gauge when Δ = 0: c_j → e^{−ijφ} c_j removes it on the infinite chain (a twist on rings), and moves the pairing to Δ e^{2iφ j}. This is the origin of the asymmetric gapless phases below.

### 1.3 Solvable submanifolds inside the class

- Free locus V = 0 (J_z = 0): 3-dimensional, section 2.
- U(1) locus Δ = 0 (J_x = J_y): XXZ with field and z-DM, Bethe ansatz.
- XYZ at zero field (D = h_z = 0): Baxter; and its free sub-cases with one coupling zero (J_x = 0 or J_y = 0 at D = h_z = 0, free by a Jordan-Wigner transformation along the other axis; Miao et al., PRL 118, 267701).
- Baxter's second integrable line J_x = −J_z inside XYZ (Chepiga-Mila), mapped to XXZ by a sublattice rotation.
- Frustration-free (Peschel-Emery, Katsura-Schuricht-Takahashi) surface, at D = 0: 4Δ² + (μ − V)² = (V + 2t)² in (1.2) units [literature, Chepiga-Mila eq. above their (3)], a disorder line inside the gapped Z2 phase with an exact twofold product ground state.

## 2. The free locus and period-2 free chains

### 2.1 Uniform free chains

At V = 0 the Bogoliubov band is

    E_+(k) = 2 t sin φ sin k + sqrt( (2 t cos φ cos k − μ)² + 4 Δ² sin² k ),   E_−(k) = −E_+(−k).      (2.1)

The many-body excitations are one complex quasiparticle per momentum, with energy E_+(k), k ∈ (−π, π].

Proposition F [proved; C143 Prop. 2.4 in other variables, re-derived here]. (i) The chain is gapless iff either |Δ| ≥ t|sin φ| and μ = ±2t cos φ (Ising surfaces, Majorana at k = 0 or π), or t|sin φ| > |Δ| and μ² ≤ 4(t² − Δ²) (the chiral gapless phase). (ii) c ≤ 1 at every gapless point.

Proof of (ii) by zero counting. E_+(k) E_−(k) = F(k) ≡ 4t² sin²φ sin²k − (2t cos φ cos k − μ)² − 4Δ² sin²k, a polynomial of degree 2 in u = cos k. Unless F ≡ 0, F has at most 2 roots u, hence at most 4 zeros on the k circle with multiplicity. Every zero of E_+ at k gives a zero of E_− at −k, so E_+ has at most 2 zeros with multiplicity. A low-energy theory of a gapless free chain with linear crossings consists of Majorana-Weyl fermions, one per simple zero of det H_BdG(k) = −F(k) (a zero at k ≠ 0, π of E_+ and its partner at −k together form one chiral complex fermion; a Dirac touching at k ∈ {0, π} is one left and one right Majorana-Weyl). Hence c_L + c_R ≤ (1/2) · 4 = 2 and c ≤ 1. F ≡ 0 forces c = μ = 0 and t|sin φ| = |Δ|, which is the XZ point (flat band, C143). QED.

The possible configurations are: two Fermi points (c = 1, the chiral phase); one or two Dirac touchings at k = 0, π (c = 1/2, or c = 1 at μ = 0, cos φ = 0 where the touchings at 0 and π combine into one Dirac fermion); an over-tilted touching plus one Fermi point (c = 1, a measure-zero line inside the chiral phase); quadratic touchings (Pokrovsky-Talapov surfaces |Δ| = t|sin φ| and μ² = 4(t² − Δ²)); and the one-sided cubic junction J where the Ising surface meets |Δ| = t|sin φ| (section 6).

### 2.2 Period-2 free chains (Lemma P2)

Lemma P2 [proved; numerical check period2_det.py]. Let H be any quadratic fermion chain with one orbital per site, nearest-neighbour hopping and pairing and on-site potential, all allowed to alternate with period 2 (arbitrary complex amplitudes). Then det H_BdG(K), K the momentum of the doubled cell, is an even trigonometric polynomial of degree ≤ 2, so c ≤ 1 at every gapless point.

Proof. In the Majorana basis (γ_1, γ_2 on site A, γ_3, γ_4 on site B of a cell) the Bloch matrix is A(K) = [[m_A ε, X(K)], [−X(−K)^T, m_B ε]] up to a factor i, with ε = [[0, 1], [−1, 0]], X(K) = X_0 + X_1 e^{−iK}, X_0, X_1 real 2 × 2 (intra- and inter-cell couplings). For m_A ≠ 0 the Schur complement gives det A = det( m_A m_B ε + X(−K)^T ε X(K) ) (using ε^{−1} = −ε and det(m_A ε) = m_A²). With the 2 × 2 identity M^T ε M = (det M) ε and N ≡ X_1^T ε X_0,

    X(−K)^T ε X(K) = (det X_0 + det X_1) ε + e^{iK} N − e^{−iK} N^T = (d_0 + d_1) ε + 2 cos K n_a ε + 2i sin K n_s,

where N = n_s + n_a ε splits into its symmetric part n_s and antisymmetric part n_a ε. So det A = det( m' ε − 2i sin K n_s ) with m' = m_A m_B − d_0 − d_1 − 2 n_a cos K, and a 2 × 2 computation gives

    det A(K) = (m_A m_B − d_0 − d_1 − 2 n_a cos K)² − 4 sin² K det n_s,

a quadratic polynomial in cos K. The case m_A = 0 follows by continuity. The zero counting of 2.1 then applies (4 Majorana-Weyl zeros at most, c ≤ 1). QED. [numerical] Over 500 random models the Fourier modes |n| ≥ 3 of det H_BdG vanish to 4e-16 relative; random sampling (20,000 models) and an adversarial search (60 starts) find max c_free = 1 (period2_mf.py, period2_search.py).

The uniform free chain is the special case of P2, so 2.1(ii) is a corollary.

## 3. Mean-field closure [derived]

Decouple V n_j n_{j+1} in every channel (Hartree, Fock, Bogoliubov), allowing translation-invariant or period-2 order (CDW, bond order, alternating pairing). Each decoupled term is a bilinear on sites j, j + 1 or on-site, so the mean-field Hamiltonian is a period ≤ 2 nearest-neighbour quadratic chain, and by Lemma P2 every Hartree-Fock-Bogoliubov state of the class has c_MF ≤ 1. This matters because the known c = 3/2 Ising-plus-Luttinger phase of interacting Majorana chains (Rahmani, Zhu, Franz, Affleck, PRB 92, 235123 (2015)) appears already at mean-field level, where the quartic term γ_j γ_{j+1} γ_{j+2} γ_{j+3} generates a range-three Majorana hopping and a second harmonic in the band. In the qubit language that quartic term contains X_{j−1} X_{j+1}, a three-site term. In the parity class the corresponding mechanism is pushed beyond mean field (section 6).

## 4. Bosonization around the U(1) plane [derived]

At Δ = 0 the chain is XXZ in a field with z-DM, gauge-equivalent to XXZ with J' = sqrt((J_x + J_y)² + 4D²)/2 per component and a twist. Its gapless phase is a Luttinger liquid with filling n, parameter K (K = 1 free, K = 1/2 at the isotropic point at half filling, K → 1/4 at the commensurate-incommensurate edge of the period-2 phase for V > 2t) [literature, Bethe ansatz]. With the convention that e^{2iΦ} is the 2k_F density (dimension K) and e^{iθ} the boson creation operator (dimension 1/(4K)), the pairing becomes

    Δ c_j c_{j+1} → −2i Δ sin k_F e^{2iφ x} ψ_R ψ_L ∝ Δ e^{2iφ x} e^{2iθ} + (oscillating terms with e^{2inΦ}, n ≠ 0),     (4.1)

with scaling dimension 1/K + n²K and conformal spin 2n·(…) ≠ 0 for n ≠ 0, carrying momentum 2φ + 2n k_F (in the gauge where the Fermi points are ±k_F).

- φ = 0 (no DM): the n = 0 term is commensurate. For K > 1/2 it is relevant and opens the gapped Z2 (topological, spin-ordered) phase; for 1/4 < K < 1/2 it is irrelevant and the Luttinger liquid survives Δ ≠ 0 as the floating phase (c = 1). Its boundaries are KT (K = 1/2) and Pokrovsky-Talapov (the period-2 edge). This is the Chepiga-Mila phase diagram [literature, confirmed by their DMRG to L = 3001].
- φ ≠ 0, π: the n = 0 term carries momentum 2φ. Shifting θ → θ − φ x turns it into a commensurate cosine plus a twist term (vK φ/π) ∂_x θ. This is a commensurate-incommensurate problem in θ: the Luttinger liquid survives up to Δ_c(φ, K, n) > 0, where the pinning energy of cos 2θ balances the twist, and the transition is Pokrovsky-Talapov. For K < 1/2, Δ_c is not perturbatively small. The free limit is exact: Δ_c = t|sin φ| (Prop. F). Operators with n ≠ 0 have nonzero conformal spin and cannot gap the phase on their own. At second order they only renormalize K and generate the usual umklapp e^{4iΦ} at half filling.
- In every case the low-energy theory is one compact boson. Hence c ≤ 1 in a neighbourhood of the U(1) plane, including all its gapless (chiral, floating) continuations. The asymmetric gapless phase has v_L ≠ v_R (the Fermi points sit at φ ± k_F) but c_L = c_R = 1.

Multicritical points reached from the U(1) side: the KT-PT meeting at the particle-hole line is Baxter's J_x = −J_z line (8-vertex universality, c = 1, ν = π/(2ρ), cos ρ = J_y/J_x) [literature, exact plus DMRG]. Where the Ising lines come down to Δ = 0 they end at the z = 2 band edges of the U(1) chain (saturation, empty band), single species.

## 5. Neighbourhood of the zero-field XYZ critical lines [derived]

At D = h_z = 0, XYZ is critical iff two couplings have equal modulus not smaller than the third (|J_a| = |J_b| ≥ |J_c|, up to the sublattice signs of 1.1(b)), with c = 1 (Gaussian, Baxter's ρ) [literature]. On a line with J_x = ±J_z (U(1) about y after a sublattice rotation), h_z is a transverse field, e^{±iθ_y} of dimension 1/(4K_y) < 2, relevant. It gaps the line except where it is fine-tuned against a second relevant operator. The resulting transitions are Ising (Dmitriev-Krivnov-Ovchinnikov, XXZ in a transverse field), or, inside the particle-hole wedge of section 4, the floating phases. D enters as a twist of the y-U(1) only when J_x = J_y. Otherwise it is a spin-1 operator (∂Φ-type, current) and adds no new sector. One boson plus relevant perturbations gives c ≤ 1 nearby.

## 6. Ising junctions, the rigidity of the free band, and the only weak spot

This is the multicritical locus flagged in the spec: where an Ising surface (Majorana at k = 0) meets a Luttinger region. In the free chain it is the curve J = {μ = 2t cos φ, |Δ| = t|sin φ|} (and its k = π image).

### 6.1 Local form

Near an Ising point at k = 0 write E_+(k) = ε_a(k) + R(k), with ε_a odd (inversion-breaking tilt, from D) and R ≥ 0 even. Particle-hole symmetry and the oddness of p-wave pairing force R(k) = r_1|k| + r_3|k|³ + O(k⁵) (no |k|² term, since ε_s ∝ k² and Δ(k) ∝ k at an Ising point), and ε_a(k) = a_1 k + a_3 k³ + … . On the side where the tilt opposes the cone,

    E_+(k) = s |k| + c_3 |k|³ + O(k⁵),   s = r_1 − |a_1|,   c_3 = r_3 − sgn(a_1) a_3.     (6.1)

[derived] A c = 3/2 Ising × Luttinger point requires E_+ to keep the touching at k = 0 and to change sign twice more. Near J this needs s > 0 small and c_3 < 0: then E_+ crosses zero at k_0 = sqrt(s/|c_3|) and, by periodicity, again at some k_1, giving a Majorana (c = 1/2) plus a chiral Dirac pair (c = 1).

### 6.2 Free value and margin [proved]

From (2.1) at μ = 2t cos φ, expanding sqrt(t² cos²φ k⁴ + 4Δ²(k² − k⁴/3)):

    E_+(k) = 2(Δ − t sin φ)|k| + [ t² cos²φ/(4Δ) − Δ/3 + t sin φ/3 ] |k|³ + O(k⁵),

so at the junction (Δ = t sin φ) the free cubic coefficient is c_3 = t² cos²φ/(4Δ) > 0. It vanishes only for cos φ = 0, where the junction coincides with the XZ flat-band point (μ = 2t cos φ = 0, Δ = t). This is the one-sided cubic onset of C143 Prop. 2.4(iii), and it is consistent with Prop. F, since the touching at 0 plus two more zeros would exceed the bound of 2.

### 6.3 Interactions [derived]

- First order (Hartree-Fock-Bogoliubov): the band stays in the free family with renormalized (t, φ, Δ, μ) (section 3), so c_3 = t_HF² cos²φ_HF/(4Δ_HF) ≥ 0 at the HF junction. No c = 3/2 at first order.
- Beyond first order, interactions add higher harmonics to the quasiparticle dispersion (second-neighbour odd hopping a_2 sin 2k and similar), at order V²/t. They shift c_3 by O(V²/t). A sign change needs V²/t ≳ t_HF² cos²φ_HF/Δ_HF. So the c = 3/2 junction mechanism can only operate within O(|V|) of the curve where the Hartree-Fock band is the XZ flat band (cos φ_HF = 0, μ_HF = 0, Δ_HF = t_HF). That curve starts at the XZ point and is the image of "XZ + V" (in C144's frame XZ + β YY plus the mixing directions). The same rigidity argument covers the other ways to exceed c = 1 (four Fermi points, or a touching plus a Fermi pair at k = π). Each needs E_+ to have more zeros than a degree-2 F allows, so it needs a perturbation of the band at least as large as the margin by which the free band avoids the corresponding degenerate configuration. That margin is bounded below on every compact subset of the free locus avoiding XZ. [derived; the margin function is explicit for the junction, c_3 = t² cos²φ/(4Δ), and is only argued, not computed, for the other degenerations.]
- Near XZ the free band splits into a dispersing half (A, k ∈ (0, π), a critical Majorana) and a flat half (B, k ∈ (−π, 0)). With (cos φ, μ, Δ − t sin φ) of order (V, V, V²) the B half is a strongly interacting single-Majorana-per-site chain of bandwidth O(V²). This is the C144 regime, where C144 finds c_tot ≤ 1 on the Ashkin-Teller plane (cos φ = μ = 0). Its B chain is of Rahmani type (one Majorana per site, quartic interaction with a 1/r tail), and that family does contain a c = 3/2 Ising × Luttinger phase for one sign of the quartic-to-hopping ratio. So the scaling window around XZ + V is the only place where the c ≤ 1 argument is not closed at the physics standard. Numerics there: section 8.

## 7. Strong coupling [derived]

- V → +∞ (J_z ≫ |J_x|, |J_y|, |D|): at density n < 1/2 the low-energy space has no adjacent particles. Pairing creates adjacent pairs at cost V and acts only virtually, generating range-2 pairing t Δ/V and correlated hopping. Hard rods with nearest-neighbour exclusion are a single-species Luttinger liquid with K = (1 − n)² [literature, hard-core exclusion], so the pairing (dimension 1/K) is irrelevant for n > 1 − 1/√2. At n = 1/2 the period-2 CDW (z-Néel) is gapped, and its doped domain walls are again one species of spinless fermions with pairing (the domain-wall theory of commensurate melting, Chepiga-Mila). One species gives c ≤ 1.
- V → −∞ (ferromagnetic J_z): the empty and full states are degenerate at μ = V. An n-particle cluster costs (n − 1)V − nμ, so for μ < V the empty state is stable against every cluster, and at μ = V the transition is first order (phase separation). Single-particle and pair excitations have gaps of order |V|. A pair fluid that coexists with a single-particle fluid would need μ_pair/2 = μ_single over an interval. The cluster energies (n − 1)V make the largest cluster win instead, so no two-fluid phase arises [heuristic; for Δ = 0 this is the Bethe-ansatz statement that the XXZ ground state has no bound strings for J_z/J ≤ −1, where the chain is ferromagnetic]. At intermediate attraction the Ising lines can turn first order through tricritical Ising points (c = 7/10) [literature for interacting Majorana chains; here numerical, section 8].

## 8. The Ashkin-Teller plane, the duality, and the two-Majorana picture (resumed run)

### 8.1 Majorana form on the plane cos φ = 0, h_z = 0 [proved]

Write c_j = (a_j + i b_j)/2. With w = i t (φ = π/2) the hopping is i t (X − X†), X = c†_{j+1} c_j, and a direct expansion gives i t (X − X†) = −(i t/2)(a_j a_{j+1} + b_j b_{j+1}) (the cross terms a_{j+1} b_j + b_j a_{j+1} vanish by anticommutation). The pairing gives −Δ(c_j c_{j+1} + h.c.) = −(i Δ/2)(a_j b_{j+1} + b_j a_{j+1}). In the rotated Majoranas A_j = (a_j + b_j)/√2, B_j = (a_j − b_j)/√2 one has a a' + b b' = A A' + B B' and a b' + b a' = A A' − B B', so the quadratic bond is

    −(i/2)[(t + Δ) A_j A_{j+1} + (t − Δ) B_j B_{j+1}] ≡ −t_A (i A_j A_{j+1}) − t_B (i B_j B_{j+1}),   t_A = (t + Δ)/2, t_B = (t − Δ)/2.     (8.1)

For the interaction, n_j = (1 + i a_j b_j)/2 = (1 − i A_j B_j)/2 (since a_j b_j = −A_j B_j). At μ = V (h_z = 0) the bilinear parts of V n_j n_{j+1} − μ n_j cancel in the bulk, and the quartic part is (V/4)(i A_j B_j)(i A_{j+1} B_{j+1}) = −(V/4)(i A_j A_{j+1})(i B_j B_{j+1}). So

    H_AT = Σ_j [ −t_A ε^A_j − t_B ε^B_j − (V/4) ε^A_j ε^B_j ],   ε^A_j = i A_j A_{j+1},  ε^B_j = i B_j B_{j+1}.     (8.2)

Each species alone is a critical transverse-field Ising chain (its bonds alternate between the σσ and σ^x terms of the Ising spins), and the coupling is the product of the two bond energies on every bond. This is the self-dual quantum Ashkin-Teller chain with unequal Ising couplings for the two species (Kohmoto, den Nijs and Kadanoff, PRB 24, 5229 (1981) for the equal case). In the frame of C144 (parity axis along y) the plane is XZ + α ZX + β YY, the Ashkin-Teller plane of C144, with t_B/t_A ↔ α and V ↔ β.

Symmetries of (8.2) [proved]: the parities P_A, P_B of the two species are separately conserved; translation by one site acts on each species as its Kramers-Wannier duality, so ε^A → −ε^A and ε^B → −ε^B in the continuum (the Ising mass of each species is odd).

### 8.2 Duality [proved, checked numerically]

The automorphism B_j → (−1)^j B_j, A_j → A_j maps t_B → −t_B and V → −V, i.e. (t, Δ, V) → (Δ, t, −V): the point with ratio s = Δ/t and interaction V is unitarily equivalent to the point with ratio 1/s and interaction −V (after rescaling by 1/s). In spin language J_x ↔ D with J_y = −J_x and J_z → −J_z. Check (dual_check.py): the complete ring spectra in both Π sectors agree to 2e-14 at L = 8 for (s, V) = (0.9, −0.8), (0.5, 0.6), (0.3, −1.3); at L = 10 they differ, because on a ring with L ≡ 2 mod 4 the sign (−1)^j changes the boundary condition of the B species. So on the infinite chain the plane is covered by s ∈ [0, 1] with both signs of V, and s = 1 (the XZ line t_B = 0) is self-dual up to V → −V.

### 8.3 Continuum theory [derived]

The A chain has right movers at k = 0 and left movers at k = π (dispersion 2 t_A sin k), so A_j ≈ A_R(x) + (−1)^j A_L(x), and likewise for B when t_B > 0. For t_B < 0 (s > 1) the roles of k = 0 and π are exchanged for B, i.e. B_R and B_L swap.

(a) s < 1 (t_A, t_B > 0). The pair (A, B) is one Dirac fermion ψ = (A + iB)/√2 with unequal velocities, bosonized by a field φ with Luttinger parameter K set by V (K = 1 at V = 0). The operators allowed by P_A, P_B and translation are: ε^A ε^B, which contains (∂φ)² and sets K (marginal); the velocity anisotropy (t_A − t_B)(T_A − T_B) ∝ ψ_R ∂ψ_R + h.c., which bosonizes to a vertex with conformal spin 2 and scaling dimension K + 1/K, hence irrelevant for K ≠ 1 and exactly marginal at the free point, where it is just the free two-velocity theory (c = 1); and the umklapp-type vertices cos 4φ (dimension 4K) and cos 4θ (dimension 4/K). The single-species masses ε^A, ε^B and the combinations ε^A ± ε^B are odd under translation and forbidden. So the plane is c = 1 for 1/2 < K < 2 and gapped beyond, through KT transitions. On the equal-velocity line s = 0 (Δ = 0) this is the XXZ chain with J_z/J = V/(2t), whose critical range −1 < J_z/J ≤ 1 contains 1/2 ≤ K < ∞, and K = 2 at J_z/J = −1/√2, matching the known critical range of the self-dual Ashkin-Teller chain [literature]. The low-energy theory is one compact boson everywhere in this region, so c ≤ 1.

(b) Off the plane. The hopping phase and the field add the lattice bilinears i(A_j B_{j+1} − B_j A_{j+1}) (coefficient ∝ t cos φ) and i A_j B_j (coefficient ∝ h_z, the chemical potential). With the decomposition above, their uniform continuum parts are i(A_R B_R − A_L B_L) and i(A_R B_R + A_L B_L); the parts that mix chiralities carry (−1)^j and average out. For s < 1 these are the current and density of the Dirac fermion (dimension 1, spin 0, exactly marginal): they shift the Fermi points, make the phase incommensurate, and leave c = 1. This is the chiral gapless phase of Prop. F continued to V ≠ 0, and it contains the Ashkin-Teller c = 1 region as its codimension-2 slice. For s > 1 the same bilinears pair A_R with the left-moving B mode, so they are masses m_± ∝ μ ± 2t cos φ for the two Majorana combinations: generically gapped, an Ising surface (c = 1/2) where one mass vanishes, and the c = 1 Ashkin-Teller line only on the codimension-2 plane where both vanish. This reproduces Prop. F at V = 0 and continues it to V ≠ 0. It is also the mechanism behind the splitting of a c = 1 line into two Ising lines by the DM term reported by Peotta et al. [literature, J. Stat. Mech. (2014) P09005].

(c) Near s = 1 (t_B → 0) the B species becomes slow and V is never small compared with its bandwidth. This is the XZ regime of C144, where the B chain is a single strongly interacting Majorana chain. C144 finds on the plane (in its variables) c_tot = 1 outside a wedge, 0 inside, and 1/2 on one line inside the wedge, with no fit above 1.09 at L = 128 [numerical, C144]. By the duality 8.2 the two signs of V are equivalent, so C144's half-plane covers the whole plane.

### 8.4 Consequence for the multicritical loci of the spec

The "Ising line meets Luttinger region" configuration of the spec occurs in this class in two ways. (i) Two Ising surfaces (Majoranas at k = 0 and k = π, or the two combinations m_±) cross: the crossing is the Ashkin-Teller c = 1 line, one boson made of the two Majoranas, not Ising × boson. (ii) An Ising surface meets the chiral c = 1 phase: at V = 0 and at Hartree-Fock-Bogoliubov level this is the one-sided cubic junction of section 6 with c ≤ 1 by the zero count. Beyond first order the rigidity margin c_3 = t² cos² φ/(4Δ) protects it except within O(|V|) of XZ. In the two-Majorana language the reason is that the Ising Majorana and the chiral Dirac fermion are built from the same two Majoranas per site, so a c = 3/2 point would need six Majorana-Weyl components at low energy (three per chirality) from two Majorana species per site, while the quadratic theory allows at most four (Lemma P2) and which interactions can produce only through a strongly coupled slow sector (section 6.3, the XZ window).

### 8.5 The particle-hole-symmetric hyperplane h_z = 0 off the plane [derived, to be confirmed by job 17597010]

At h_z = 0 the composition of the π rotation about x (D, h_z → −D, −h_z) with site inversion (D → −D) is a symmetry, so h_z = 0 is a symmetric hyperplane, like D = 0. Near the Ashkin-Teller plane on the side where the effective ratio s_eff = Δ_eff/(t sin φ)_eff exceeds 1, the hopping-phase bilinear of 8.3(b) is a mass of opposite sign for the two Majorana combinations, i.e. the antisymmetric Ashkin-Teller perturbation, a vertex of dimension 1/K or K depending on the side (the two are exchanged by the duality 8.2). Where that vertex is irrelevant, a c = 1 region survives inside the hyperplane h_z = 0, while h_z (the symmetric combination) gaps it at once. This is the analogue of the floating phase of the D = 0 hyperplane (section 4), one boson, c = 1. The L = 64 data at φ = 1.45, s = 0.9 and 1.0, V = −0.8 fit this pattern: gapless at h_z = 0 (S_mid = 1.07, 1.06), gapped already at h_z = 0.05 (S_mid = 0.64, 0.67, fits 0.03 to 0.12). The fitted c at h_z = 0 (1.13 to 1.17 on the central half, 1.31 to 1.40 on the wide window) is above 1 at L = 64, and the free chain at the same (φ, s) fits 0.99 to 1.01 at L = 64 (free_open.py), so the excess is an interaction effect. The finite-size scaling at L = 32, 64, 128 is in job 17597010; until it returns, c = 1 on this hyperplane near XZ is derived, not numerically confirmed.

## 9. Numerics (calibrated)

All open-chain runs use dmrg_pt.py (C144's NumPy DMRG with the full field restored on the two edge sites), fits of S(l) to (c/6) log[(2L/π) sin(πl/L)] on the central half (c25) and the central three quarters (c12), each parity separately and for the bond-averaged profile. Convergence: truncation below 1e-9 and last-sweep energy change below 1e-6 for every run quoted (most below 1e-9).

9.1 Calibration [numerical]. At L = 64: Baxter 8-vertex point (exact c = 1) 0.98 / 0.97; free chiral phase (c = 1) 0.97 / 0.96, against 0.96 / 0.97 from the exact free-fermion profile at the same L (free_open.py); free Ising (c = 1/2) 0.56 / 0.58, the usual open-boundary excess at L = 64. Free chains at the near-XZ points with h_z = 0 give 0.98 to 1.01 at L = 64 and 0.99 at L = 256, independently of the velocity ratio (free entanglement does not see velocities).

9.2 Candidates from ring ED at φ = 0.9 (c_eff up to 2 to 6 at L = 12 → 16). DMRG (L = 64, 96, χ = 96): averaged fits 0.86 to 1.13, wide-window fits 0.96 to 1.00 at L = 96, pair estimates c(64, 96) = 0.81 to 1.12. All eight are c = 1 Luttinger points. The ring values were level-crossing artifacts at small L [numerical].

9.3 Near-XZ family, L = 64 (job 17594701, 264 runs; φ = 1.3, 1.45, s = Δ/(t sin φ) = 0.9, 1.0, 1.1, V = ±0.4, ±0.8, h_z = 0 to 0.5; table from xzj_table.py). Repulsive V: at φ = 1.3, c = 1 throughout (fits 0.88 to 1.16, scatter from parity oscillations), except the gapped large-h_z end at s = 1.1, V = 0.4. At φ = 1.45 the fits are 1.0 for h_z ≥ 0.15 but rise at small h_z, most for s ≥ 1: 1.55 / 1.54 at (s, V, h_z) = (1.1, 0.4, 0), 1.26 / 1.46 at (1.1, 0.4, 0.05), 1.28 / 1.30 at (1.1, 0.8, 0), 1.23 / 1.26 at (1.0, 0.4, 0). These points lie in the XZ window: at s_eff ≈ 1 the slow half of the free band has width about 2t cos² φ ≈ 0.03, far below V, which is the strongly coupled B regime of C144. There C158 found crossover readings of 2 to 6 at moderate bond dimension that converge to 1 once the slow incommensurate mode dominates. The L = 64 fits there are therefore not evidence for c > 1. They are evidence that L = 64 is not asymptotic, and job 17597010 runs the same points at L = 32, 64, 128 (192 for the two largest). Attractive V: gapped almost everywhere (S_mid 0.2 to 0.6, fits 0.00), with one transition per line in h_z between 0.05 and 0.15. At V = −0.8 its fits are 0.2 to 0.53 (Ising, c = 1/2, as expected for a narrowly missed transition on a 0.05 grid). At V = −0.4, φ = 1.3, s = 1.0 and 1.1, h_z = 0.15 both windows and both parities give 0.87 to 0.90. That is either a narrow c = 1 window (a chiral sliver between two Pokrovsky-Talapov edges, open in the class) or a crossover from a nearby multicritical point. In neither case is it above 1. At φ = 1.45, V = −0.4 the analogous transition sits near h_z = 0.05 (fits 0.88 / 0.87 at s = 1.0, 0.81 / 0.81 at s = 1.1, 1.02 / 1.07 at s = 0.9). The h_z scans at 0.01 resolution with L = 32, 64, 128 are in job 17597010. The h_z = 0 points at φ = 1.45 are discussed in 8.5.

9.4 Small-ring scans (cancelled job 17594696, slices at V = 6, φ = 0 and 0.5, L = 12, 16, 20). The Luttinger region shows c_eff 0.7 to 1.4 with the filling-dependent shell oscillation of incommensurate liquids. The CDW and polarized regions give 0, the Ising edges give about 0.5, and isolated spikes of 3 to 5 sit at single points on the CDW edge (μ = 3.33), where the ring ground state changes sector. No extended region has c_eff above 1.4.

9.5 Pending (job 17597010): finite-size scaling at the h_z = 0 near-XZ points (8.5), the fine h_z scan at the c ≈ 0.87 transitions (9.3), a strong-attraction μ scan (V = −1.5, Δ = 0.5, φ = 0 and 0.8) aimed at tricritical Ising points and at any pair-plus-particle coexistence, and the region at φ = 1, h_z = 0, V ≈ −2t, Δ ≈ 0.3t where the phases of Peotta et al. meet and which they could not resolve. Analysis: python3 fss.py typhon_out/j2/dmrg_*.jsonl after copying /scratch/lustre/jkudlerflam/c152/j2/dmrg_*.jsonl. The decisive numbers are the pair estimates c(L, 2L) = 6[S̄(2L) − S̄(L)]/ln 2, with S̄ the average over the four central cuts.

## 10. Phase diagram of the class: the gapless loci and their central charges

Coordinates: fermion units (t, φ, Δ, V, μ) up to scale, h_z = (μ − V)/2, s = Δ/(t |sin φ|). Symmetries 1.1(a) to (d) reduce to φ ∈ [0, π/2], Δ ≥ 0, h_z ≥ 0. Dimension counts are within the 4-dimensional class.

Open gapless phase.
- G1 [c = 1]. The chiral (inversion-asymmetric) Luttinger liquid, v_L ≠ v_R, Fermi points at φ ± k_F. Contains the U(1) liquid (Δ = 0, Bethe ansatz), the free chiral phase t |sin φ| > Δ, μ² < 4(t² − Δ²) (proved), its continuation to V ≠ 0 (derived, sections 4 and 8.3(b)), and, as its D = 0 slice, the floating phase (1/4 < K < 1/2, Chepiga-Mila). Status: c = 1 proved at V = 0, exact on the U(1) and integrable loci, derived elsewhere, numerical at the points of 9.2 and 9.3.

Codimension-one critical surfaces.
- I [c = 1/2]. Ising surfaces: the topological/trivial (x-ordered/paramagnetic) transitions continuing μ = ±2t cos φ for s > 1 (proved at V = 0, derived for small V by the mass argument of 8.3(b)); the Ising lines between the period-2 CDW or Z2 phases and disordered phases on the repulsive side (Chepiga-Mila, literature); the pair of Ising lines between the z-ferromagnet, a disordered phase and the x-ferromagnet on the attractive side (Peotta et al., literature); the attractive-side transitions of 9.3 (numerical).
- PT [z = 2, not conformal]. Pokrovsky-Talapov and band-edge surfaces bounding G1: band filling and emptying, the commensurate-incommensurate edges of the CDW, and the pairing threshold Δ = Δ_c(φ, K, n) (Δ_c = t |sin φ| at V = 0). The adjacent gapless phase has c = 1.
- KT [c = 1 on the surface]. Kosterlitz-Thouless boundaries of G1 within the symmetric hyperplanes (D = 0: K = 1/2 between floating and Z2; h_z = 0: the analogous edge of 8.5; the Ashkin-Teller plane: K = 1/2 and K = 2).
- First-order surfaces (no conformal point): Luttinger liquid to ferromagnet on the attractive side (Peotta et al.), phase separation at V → −∞ (section 7).

Codimension ≥ 2.
- AT [c = 1]. The crossing of the two Ising surfaces on the plane cos φ = 0 = h_z with s > 1: the Ashkin-Teller line, one boson (8.3, 8.4). By the duality 8.2 it is the image of the s < 1 part of the plane.
- XYZ [c = 1]. Zero-field XYZ critical lines |J_a| = |J_b| ≥ |J_c| (Baxter, exact), including the 8-vertex point J_x = −J_z where the KT and PT boundaries of the floating phase meet (Chepiga-Mila, exact plus DMRG).
- J [Lifshitz, c ≤ 1]. The curves where an Ising surface meets the boundary of G1: one-sided cubic onset at V = 0 (proved, section 6.2), c ≤ 1 at first order in V (proved via Lemma P2) and beyond first order away from XZ (derived, margin c_3 = t² cos² φ/(4Δ)).
- TCI [c = 7/10, possible]. Points where an Ising surface turns first order on the attractive side. They are established in interacting Majorana chains (Rahmani et al.) and in extended Kitaev chains; in this class they are not yet located (job 17597010, Set 3). Either way c = 7/10 < 1.
- XZ [singular point, flat band]. Emergent c_tot ≤ 1 in every direction checked (C143, C144, C154, C158), with the h_z = 0 window of 8.5 pending.
- FM_H [z = 2]. The ferromagnetic Heisenberg point J_x = J_y = −J_z (V = −2t at Δ = D = 0), quadratic magnons, the endpoint where G1 meets the ferromagnets.
- FF. The frustration-free (Peschel-Emery, Katsura-Schuricht-Takahashi) surface is a disorder surface inside the gapped Z2 phase, not critical.

No locus with c > 1 has been found or derived.

## 11. Theorem-style summary for the parity class

Let H be in the parity class, i.e. (1.2) with parameters (t, φ, Δ, V, μ).

Theorem A [proved]. (i) If V = 0, then at every gapless point c ≤ 1, and c = 1 exactly on the chiral phase and the Dirac touchings, 1/2 on the Ising surfaces (Prop. F). (ii) For any period ≤ 2 Hartree-Fock-Bogoliubov state of the class, the mean-field Bogoliubov spectrum has c_MF ≤ 1 (Lemma P2 and section 3). (iii) On the plane cos φ = 0 = h_z the model is the asymmetric self-dual Ashkin-Teller chain (8.2) and the map (t, Δ, V) → (Δ, t, −V) is an exact duality (8.2).

Theorem B [exact results from the literature]. c = 1 at every gapless point of the U(1) locus Δ = 0 (Bethe ansatz) and of the zero-field XYZ locus D = h_z = 0 (Baxter), and in particular at the 8-vertex multicritical point.

Theorem C [derived at the theoretical-physics standard]. c ≤ 1 at every conformal point of the class in each of the following regions: (a) a neighbourhood of the free locus away from its singular set (the flat-band point XZ and the z = 2 surfaces), by the c-theorem applied to the free CFT (c ≤ 1) perturbed by V, with no new low-energy modes from bands gapped by more than |V|; (b) a neighbourhood of the U(1) locus, including the chiral and floating continuations (section 4); (c) the Ashkin-Teller plane and its neighbourhood, where the low-energy theory is one boson built from the two Majoranas per site, or a single Majorana on an Ising surface (8.3); (d) the Ising-chiral junctions J away from O(|V|) of XZ (section 6); (e) strong repulsion, where hard-core exclusion leaves one species of particles or CDW domain walls (section 7); (f) strong attraction, where cluster energies force phase separation and leave no two-fluid phase (section 7, heuristic for Δ ≠ 0).

Numerical status [numerical]. Open-chain DMRG to L = 96 (and calibrated free-fermion profiles to L = 256) at 283 points on the critical path, including all ring-ED candidates and the near-XZ family, finds no c > 1 outside the XZ window. Inside the window (φ = 1.45, small h_z) the L = 64 fits reach 1.40 to 1.55; the free chain there gives 1.0, the derivation gives c = 1, and these values are unresolved until the L = 128, 192 runs return. C144 and C158 cover the Ashkin-Teller plane and the XZ neighbourhood to L = 128 and in the thermodynamic limit (VUMPS), with no c clearly above 1.

What is not closed. (1) The XZ window at h_z = 0 (8.5) until job 17597010 returns; the derivation there rests on the irrelevance of one vertex operator, i.e. on K being on the right side of the threshold, which the numerics must confirm. (2) Intermediate couplings far from all controlled regions (|V| comparable with t and Δ, generic φ), where the statement rests on the one-species heuristic and on the numerics, not on a derivation. A c > 1 island there would have to be a phase with two independent soft sectors built from one fermion per site; within this class the only kinematic route, an Ising Majorana coexisting with a Dirac pair, needs six low-energy Majorana-Weyl components from two Majorana species per site, which is excluded at the free and mean-field levels and in every controlled limit above. Hence, as input A4 of the C154 conditional theorem: c ≤ 1 holds on the parity class at the derived standard in the regions (a) to (f) and at the numerical standard on the critical path sampled, with (1) pending and (2) resting on numerics.
