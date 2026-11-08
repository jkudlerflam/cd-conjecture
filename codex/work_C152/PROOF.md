# C152 PROOF: the parity-symmetric qubit class

This run was done by a Claude Opus agent (2026-11-07), not by Codex. Tags: [proved], [derived] (physics standard, RG, bosonization, controlled perturbation theory), [literature], [numerical], [heuristic], [assumed]. Scripts are in this folder; numerical outputs are named where used.

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
