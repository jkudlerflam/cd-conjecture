This work was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

# C180 derivations: field theory of the Fibonacci-cusp two-valley gas, and the literature

Standard: theoretical physics. Each paragraph opens with its tag, [proved], [derived], [numerical], [heuristic], [assumed] or [literature]. The only script is optable.py (output optable.out), which checks the band facts, enumerates the vertex operators and evaluates the scaling dimensions. Literature items were read in full text unless marked "abstract only".

## 0. Model, conventions, symmetries

C173's effective model (PROOF 3.2), in units |t₂| = Δ²/(2V) = 1, on the blockaded space (no two adjacent particles),

    H_eff = Σ_j [ (w c†_{j+1} c_j + h.c.) − (c†_{j+2} c_j + h.c.) − 2 n_j n_{j+2} − (2/3) n_j n_{j+3} − μ n_j ],   w = |w| e^{iφ}.   (0.1)

Particles never pass, so fermions and hard-core bosons give the same spectrum. Its band is ε(k) = 2|w| cos(k − φ) − 2 cos 2k.

[proved] Gauge map at φ = π/2. Put c_j = i^j d_j. Then c†_{j+1} c_j = −i d†_{j+1} d_j and c†_{j+2} c_j = −d†_{j+2} d_j, so (0.1) becomes real,

    H_d = Σ_j [ |w| (d†_{j+1} d_j + h.c.) + (d†_{j+2} d_j + h.c.) − 2 n_j n_{j+2} − (2/3) n_j n_{j+3} − μ n_j ].   (0.2)

Its band ε_d(q) = 2|w| cos q + 2 cos 2q has ε_d′ = −2 sin q (|w| + 4 cos q), so for |w| < 4 the minima are at q = ±k₀ with cos k₀ = −|w|/4, and ε_min = −2 − |w|²/4. In the c gauge the valleys are k_a = π/2 − k₀ and k_b = π/2 + k₀, so k_a + k_b = π exactly. The mass and cubic term at the minimum are 1/m = ε_d″(k₀) = 8 − |w|²/2 and ε_d‴(k₀) = 2|w| sin k₀ + 16 sin 2k₀. [numerical check, optable.out] At |w| = 0.5, 1, 2 the numerical minima agree with π/2 ∓ k₀ to 6e-6 and with ε_min to 1e-6. At |w| = 1, 1/m = 7.5 and ε‴ = −5.81.

In spin language with d† = S⁺ (magnons over the polarized state), (0.2) is a J1–J2 zigzag chain with in-plane couplings J1 = 2|w| and J2 = 2 (both antiferromagnetic, J1/J2 = |w|), an infinite nearest-neighbour Ising repulsion (the blockade), an attractive next-nearest Ising term J2^z = −J2, and a third-neighbour attraction −2/3. In the c gauge the legs (even and odd sites) are ferromagnetic Heisenberg chains, as C173 PROOF 3.3 states.

Symmetries of the full qubit chain in the parity class, H = Σ [J_x XX + J_y YY + D(XY − YX) + J_z ZZ + (h_z/2)(Z⊗1 + 1⊗Z)]:

- T, translation by one site [exact].
- P = Π_j Z_j, fermion parity [exact in the parity class]. The effective model also conserves N at order Δ²/V. ΔN = 2 terms enter at order |w|Δ/V and Δ³/V², smaller than t₂ by a factor of order Δ/V (C173 3.1). So N is an emergent U(1) broken to Z₂ at higher order.
- Θ = I·K, site inversion times complex conjugation in the Z basis [proved]. K sends Y → −Y and leaves X, Z fixed, so it flips D and fixes everything else. I maps the bond operator A_j B_{j+1} to B_m A_{m+1} with m = −j − 1, so it also flips XY − YX and fixes XX, YY, ZZ, Z. Hence Θ is an antiunitary symmetry of every parity-class bond. It preserves lattice momentum (I and K each send k → −k), so it maps each valley to itself.
- Σ = I × Π_{odd} R_z(π) × Π_j R_z(π/2) [proved]. On the bond, I sends J_x XX → J_x XX, then the odd-site π rotation gives −J_x XX, then R_z(π/2) (X → Y, Y → −X) gives −J_x YY. Likewise J_y YY → −J_y XX. For D, I gives a minus sign, the odd-site rotation a second one, and R_z(π/2) leaves XY − YX invariant, so D → D. ZZ and Z are fixed. So (J_x, J_y, D) → (−J_y, −J_x, D), and Σ is a symmetry exactly when J_x + J_y = 0, which is w = 2iD, φ = ±π/2. In fermion language Σ: c_j → (−1)^j c_{−j} up to the Jordan–Wigner bookkeeping, so c_k → c_{π−k}, and it exchanges the valleys. This agrees with C173 PROOF 3.3 (C179 is refereeing it independently).

## 1. Bosonization of the two-valley gas

[derived] At chemical potential just above ε_min both valleys are occupied with Fermi momenta k_{F,i} = π ρ_i, ρ_a + ρ_b = n. Write the lattice fermion as

    c_j ≈ e^{i k_a x} ψ_a(x) + e^{i k_b x} ψ_b(x),   ψ_i = e^{i k_{F,i} x} R_i + e^{−i k_{F,i} x} L_i,   R_i, L_i ~ U_i e^{i(θ_i ∓ φ_i)},   (1.1)

in Giamarchi's convention (ρ_i = −∂_x φ_i/π, [φ_i(x), ∂_y θ_j(y)] = iπ δ_ij δ(x − y), U_i Klein factors). A vertex operator

    O(n_a, n_b; p_a, p_b) = e^{i(n_a θ_a + n_b θ_b)} e^{−i(p_a φ_a + p_b φ_b)},   p_i ≡ n_i (mod 2),   (1.2)

adds n_i particles to valley i and carries lattice momentum

    P = n_a k_a + n_b k_b + π (p_a ρ_a + p_b ρ_b)   (mod 2π).   (1.3)

Put φ_± = (φ_a ± φ_b)/√2, θ_± = (θ_a ± θ_b)/√2 (canonical).

Remark on statistics [derived]. For ΔN odd one should use the local hard-core boson operators of the spin chain, b_j = σ⁻_j. Their bosonized form has p_i even, and the valley peaks can shift by integer combinations of πρ_a, πρ_b relative to (1.1). This only moves the resonance conditions of section 2.2 by such combinations and does not change which operators are allowed on an open set, because every ΔN-odd momentum contains a single k_a or k_b, which is generic.

### 1.1 Action of the symmetries on the fields

[derived] Translation x → x + 1 multiplies O by e^{iP}. So T acts on the infrared as the internal rotation generated by the valley charges and the densities, and a term is allowed in H only if P ≡ 0 identically on the open parameter set considered.

[derived] Σ exchanges the valleys and inverts space. Inversion sends R ↔ L, ρ(x) → ρ(−x), so φ(x) → −φ(−x), θ(x) → θ(−x). With a ↔ b,

    Σ: φ_a(x) → −φ_b(−x) + c₁,  θ_a(x) → θ_b(−x) + c₂,  and a ↔ b,   (1.4)

so φ_+(x) → −φ_+(−x), θ_+(x) → θ_+(−x), φ_−(x) → φ_−(−x), θ_−(x) → −θ_−(−x), up to constants. The constants are fixed by the lattice and only decide the sign of each allowed cosine [assumed, not computed].

[derived] Θ keeps each valley and, being antiunitary combined with inversion, maps R_i(x) → R_i(−x). Writing this out with R ~ e^{i(θ − φ)}, L ~ e^{i(θ + φ)} and i → −i gives Θ: φ_i(x) → −φ_i(−x), θ_i(x) → −θ_i(−x). Check: ρ_i(x) → ρ_i(−x), and the commutator [φ, ∂θ] = iπδ goes to −iπδ, as an antiunitary map requires.

### 1.2 The Gaussian part

[derived] Under (1.4) the gradients transform as ∂φ_+ even, ∂θ_+ odd, ∂φ_− odd, ∂θ_− even (each evaluated at −x). The Σ-invariant quadratic terms are therefore (∂φ_±)², (∂θ_±)², ∂φ_+ ∂θ_− and ∂θ_+ ∂φ_−, while ∂φ_+∂φ_−, ∂θ_+∂θ_−, ∂φ_+∂θ_+ and ∂φ_−∂θ_− are forbidden. Θ forbids none of these (every gradient is Θ-even). So at the Σ point

    H₀ = (1/2π) ∫ dx [ (u_+/K_+)(∂φ_+)² + u_+K_+(∂θ_+)² + (u_−/K_−)(∂φ_−)² + u_−K_−(∂θ_−)² + 2λ₁ ∂φ_+∂θ_− + 2λ₂ ∂θ_+∂φ_− ].   (1.5)

[derived] Origin of the drift terms. The band is not symmetric about each valley bottom (ε‴ ≠ 0), so in valley i the right and left Fermi velocities differ, v_R − v_L = ε_i‴ k_F². For one valley, (1/4π)[v_R(∂θ − ∂φ)² + v_L(∂θ + ∂φ)²] = (v̄/2π)[(∂θ)² + (∂φ)²] − ((v_R − v_L)/2π) ∂θ∂φ. Σ forces ε_b‴ = −ε_a‴, so the sum over valleys is proportional to ∂φ_a∂θ_a − ∂φ_b∂θ_b = ∂φ_+∂θ_− + ∂φ_−∂θ_+, which gives λ₁ = λ₂ = λ ≈ ε‴ k_F²/2 at weak coupling. They do not change c. They change the scaling dimensions at relative order λ²/(u_+u_−), and they matter for stability (section 3).

[derived] Scaling dimensions. For H₀ = (1/2π)∫ Xᵀ𝓗X with X = (∂φ_+, ∂φ_−, ∂θ_+, ∂θ_−) and symplectic form J, the Gaussian ground-state covariance is G = 𝓗^{−1/2} |𝓗^{1/2}(iJ)𝓗^{1/2}| 𝓗^{−1/2}, and Δ(e^{i(m·φ + n·θ)}) = (1/4) vᵀGv with v = (m_+, m_−, n_+, n_−). For one boson with 𝓗 = u diag(1/K, K) this gives G = diag(K, 1/K) and Δ = (m²K + n²/K)/4, the standard result. With λ = 0, (1.2) has

    Δ(n_a, n_b; p_a, p_b) = [ (p_a + p_b)² K_+ + (p_a − p_b)² K_− + (n_a + n_b)²/K_+ + (n_a − n_b)²/K_− ] / 8.   (1.6)

[numerical check, optable.out] The covariance formula reproduces (1.6) to all printed digits for the eight operators of the table at (K_+, K_−) = (0.35, 1.5).

## 2. Operator table

### 2.1 Operators allowed on an open set

[derived, enumeration checked by optable.py] At the Σ point, ρ_a = ρ_b = n/2 (unless Σ is spontaneously broken) and k_b = π − k_a, so (1.3) reads

    P = (n_a + n_b) π/2 + (n_b − n_a) k₀ + (p_a + p_b) π n/2   (mod 2π).   (2.1)

With k₀ and n independent and generic on an open set, P ≡ 0 identically requires n_a = n_b, p_a = −p_b and (n_a + n_b)π/2 ≡ 0 mod 2π. So n_a = n_b ∈ 2ℤ and p_a = −p_b ∈ 2ℤ. optable.py confirms that these are the only zero-momentum operators with |n_i|, |p_i| ≤ 4 (random parameter points, tolerance 1e-9). Off the Σ point (k_a + k_b ≠ π, ρ_a ≠ ρ_b) the list is empty: no vertex operator is allowed at all.

| operator | (n_a, n_b; p_a, p_b) | sector | momentum | dimension (λ = 0) | allowed where | relevant when |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| cos(2√2 φ_−) (inter-valley backscattering ψ†_{aR}ψ_{aL}ψ†_{bL}ψ_{bR}) | (0,0; 2,−2) | N-conserving | 2π(ρ_a − ρ_b), zero when balanced | 2K_− | Σ point, any class | K_− < 1 |
| cos(4√2 φ_−) | (0,0; 4,−4) | N-conserving | 4π(ρ_a − ρ_b) | 8K_− | Σ point | K_− < 1/4 |
| e^{2√2 iθ_+} + h.c. (two ab pairs, ψ_{aR}ψ_{aL}ψ_{bR}ψ_{bL}) | (2,2; 0,0) | ΔN = 4 | 2(k_a + k_b) = 2π ≡ 0 | 2/K_+ | Σ point, parity class only | K_+ > 1 |
| e^{2√2 iθ_+} cos(2√2 φ_−) | (2,2; 2,−2) | ΔN = 4 | 0 | 2/K_+ + 2K_− | Σ point, parity class | never with K_− > 1 |
| ∂φ_− (valley offset, a pseudospin field) | gradient | N-conserving | 0 | 1 | off Σ only (Σ-odd) | always (it shifts φ_−) |
| ∂φ_+∂φ_−, ∂θ_+∂θ_− (mass imbalance) | gradients | N-conserving | 0 | 2 (marginal) | off Σ only | marginal, changes K, u |

[derived] Reading of the table. In the N-conserving effective model at the Σ point the only nonlinear term is the inter-valley backscattering cos(2√2 φ_−), with dimension 2K_−. This is the g_1⊥ term of the spinful electron gas, written for the valley pseudospin, and it is the operator Hikihara et al. call g_1⊥ cos[√4π(θ_l − θ_s)] in the inner/outer Fermi-point basis (arXiv:1004.0622, their eq. 60, whose irrelevance is their TLL2 condition, eq. 61). In the parity class the two-pair term e^{2√2 iθ_+} is also allowed, because 2(k_a + k_b) = 2π. Its amplitude is of order (|w|Δ/V)² or smaller and it is irrelevant unless K_+ > 1, which needs strong attraction in the charge sector.

### 2.2 Operators allowed only on resonance hypersurfaces

[derived] Every other vertex has momentum that vanishes only on codimension-one sets. The lowest ones, with their resonance condition at the balanced Σ point, are:

| operator | (n_a, n_b; p_a, p_b) | resonance (mod 2π) | dimension (λ = 0) | relevant for |
| :-- | :-- | :-- | :-- | :-- |
| valley transfer ψ†_b ψ_a dressed by φ_+ | (−1,1; ±1,±1) | 2k₀ ≡ ∓πn | 1/(2K_−) + K_+/2 | K_− > 1/(4 − K_+), always in the c = 2 window |
| valley transfer dressed by φ_− | (−1,1; ±1,∓1) | 2k₀ ≡ 0, i.e. k₀ = π, |w| = 4 (valleys merge) | (K_− + 1/K_−)/2 ≥ 1 | K_− ∈ (2 − √3, 2 + √3) |
| pair transfer ψ†_aψ†_aψ_bψ_b, 4k_F-dressed | (−2,2; p, p), p even | 4k₀ + πnp ≡ 0 | 2/K_− + p²K_+/2 | p = 0: K_− > 1. p = ±2: K_− > 1/(1 − K_+) |
| ab pair, ΔN = 2 (parity class) | (1,1; p_a, p_b odd) | π + (p_a + p_b)πn/2 ≡ 0, e.g. n = 1/3 with p_a + p_b = ±6 | large | no |
| aa pair, ΔN = 2 (parity class) | (2,0; p_a even, p_b even) | π − 2k₀ + (p_a + p_b)πn/2 ≡ 0 | 1/(2K_+) + 1/(2K_−) + ... | K_± large |
| single particle, ΔN = 1 (off the parity class, boson form) | (1,0; even, even) | k_a + (integer) πn/2 ≡ 0 | ≥ 1/(8K_+) + 1/(8K_−) | yes, but only on its hypersurface |
| commensurate charge locking | (0,0; p, p), p even | πnp ≡ 0, n = 2j/p | p²K_+/2 | n = 1/3 (p = 6) needs K_+ < 1/9 |

[derived] Each of these, where relevant, opens a tongue of a c ≤ 1 phase around its hypersurface, of width set by its bare amplitude to the power 1/(2 − Δ). Only finitely many are relevant at given K_±, since Δ grows quadratically in the integers. So the complement of the tongues is open. This is the same structure as an ordinary floating phase, here with two incommensurate wavevectors (pattern P3 of C154). The |w| = 0 resonance of the pair transfer (4k₀ ≡ 0 at k₀ = π/2) reproduces C153's resonant conversion vertex at w = 0, where the model is c = 1 (consistency check). Resonances near the candidate [derived, numbers at |w| = 1 where k₀ = 1.8235]. The valley-transfer lines need 2k₀ ≡ ±πn, i.e. n = 0.84 or 1.16, outside the blockaded range n < 1/2. The 4k_F-dressed pair transfer needs 4k₀ + 2πn ≡ 0 with p = ±2, and 4k₀ ≡ 1.011 mod 2π, so it sits at n = 1.011/(2π) ≈ 0.161. That line lies inside the density window of the candidate (μ = −2.2 gives n ≈ 0.17). Its dimension is 2/K_− + 2K_+, which is below 2 when K_− > 1/(1 − K_+), about 1.5 for K_+ ≈ 0.34. So if the pseudospin is strongly easy-plane, a c = 1 tongue (one combination of θ_− and φ_+ locked) is centred at n ≈ 0.161, close to the candidate density, and its width is set by an unknown bare amplitude. If K_− < 1/(1 − K_+) the line is harmless. This is a specific thing for the open-chain and ring numerics to watch for: a density dependence of the c reading with a dip near n ≈ 0.16 at |w| = 1. The commensurate densities n = 1/3 and 1/5 need K_+ < 1/9 and smaller and are harmless in the dilute regime.

## 3. Stability of c = 2

[derived] At the Σ point in the N-conserving effective model, the two-boson Gaussian theory (1.5) is a stable phase (c = 2 on an open set of μ, |w|) when all of the following hold:

- (S1) The quadratic form is positive definite. With the drift terms this is u_±/K_± > 0, u_±K_± > 0, and the two 2 × 2 blocks satisfy λ₂² < (u_+K_+)(u_−/K_−) and λ₁² < (u_+/K_+)(u_−K_−). The first block is the binding one, because u_−/K_− is the small inverse valley susceptibility. Its Schur complement κ_σ⁻¹ = u_−/K_− − λ₂²/(u_+K_+) is the true inverse pseudospin susceptibility, the energy cost of a valley imbalance when the charge sector may carry current. κ_σ⁻¹ → 0 is the valley-ferromagnetic (vector chiral) instability. [numerical illustration, optable.out] At (K_+, K_−, u_+, u_−) = (0.35, 1.5, 4, 0.9) the bound is λ² < 0.84, and λ = 1.2 is unstable. The vector chiral phase has ⟨∂θ_+⟩ ≠ 0 and valley imbalance together, which is what this block describes and what Hikihara et al. find (a nonzero ⟨dθ_+/dx⟩ absorbed by a redefinition of θ_+ in their vector chiral theory).
- (S2) The inter-valley backscattering is irrelevant, K_− > 1, or marginally irrelevant at K_− = 1 with the sign of g that flows to zero (the SU(2)-like separatrix of the Kosterlitz–Thouless flow dK_−/dl ∝ −g², dg/dl = 2(1 − K_−)g).
- (S3) In the parity class, the ΔN = 4 vertex is irrelevant, K_+ < 1. This holds throughout the repulsive dilute regime (section 4, K_+ → 1/2).
- (S4) The parameters avoid the tongues of section 2.2, in particular commensurate densities with relevant charge locking.
- (S5) Both valleys are populated in the ground state. At Σ, with (S1) and (S2), this holds automatically.

[derived] Off Σ, still in the parity class. The Σ-odd term (1/π)√2 h_σ ∂φ_− with h_σ ∝ δ (the valley offset) is a field on the pseudospin. If (S1) and (S2) hold at Σ, the pseudospin sector is an easy-plane liquid, a field polarizes it continuously, the backscattering term acquires momentum 2π(ρ_a − ρ_b) and becomes even less effective, and c = 2 survives until the minority valley empties, at |h_σ| = h_sat. In the dilute limit h_sat = J_⊥(1 + Δ_σ) (the saturation field of the pseudospin XXZ chain, section 4), which is of order n³. Mass imbalance η adds marginal mixing terms and, by C164, a ferromagnetic pseudospin interaction of order η²n² per unlike pair, so easy-plane needs |η| ≲ n^{1/2} at small n. This reproduces C173's wedge |δ| ≲ n³, |η| ≲ n^{1/2} from the field-theory side. The wedge is open in the parity class.

[derived] Off the parity class (in N or T). Parity-odd terms of order ε give ΔN = 1 vertices whose momenta contain one k_a or k_b, so they are allowed only on hypersurfaces (section 2.2). At second order they also shift δ and η, which moves the point inside the wedge. So if the Σ point satisfies (S1)–(S4) with margin, a full neighbourhood in B₂ is c = 2. Nothing in the field theory requires a symmetry to protect the phase. What protects the two-boson liquid is incommensurability: T acts on the infrared as an irrational rotation, which generates an emergent U(1)_a × U(1)_b (U(1)_{a−b} × Z at the Σ point, where k_a + k_b = π).

[heuristic] How wide the region is in practice. The h_sat ∝ n³ estimate, with the C173 ring ED gaps of about 0.07 to 0.1 at N = 4, gives a wedge width in δ of order 0.1|t₂| at n ≈ 0.17. This is consistent with C173's reading persisting to φ = 1.54 (offset 0.12|t₂|) and failing by φ = 1.5.

## 4. The dilute limit: two-body data decide

### 4.1 Impenetrable gas to pseudospin chain

[derived; standard for the strongly repulsive two-component gas] At low density, inter-valley scattering is impenetrable at q → 0 unless a zero-energy resonance occurs (C173 PROOF 4(b)). Order the particles x₁ < ... < x_N and label each by its valley s_i ∈ {a, b} = {↑, ↓}. At leading order (Girardeau) the energy is that of N free fermions of mass m and density n, independent of the label sequence. The first correction comes from the short-distance structure of each neighbouring pair. The zero-energy two-body solution behaves as ψ(r) ∝ (r − 𝒜) away from contact, where 𝒜 is an excluded length that may depend on the labels. The hard-rod energy E = π²N³/(6m(L − Σ_i 𝒜_{i,i+1})²) expands to

    E = E₀ + (π² n³ / 3m) Σ_i 𝒜_{i,i+1} + O(n⁴),   (4.1)

with 𝒜_{i,i+1} an operator on the labels of the pair.

The label structure of 𝒜 [derived]. For a same-valley pair (total momentum 2k_a or 2k_b) it is a number a_s, equal for aa and bb by Σ. For an inter-valley pair (total momentum k_a + k_b = π) the two ordered channels |ab⟩, |ba⟩ (left label first) are degenerate at threshold, and 𝒜 is a 2 × 2 Hermitian matrix. Σ maps (a at x₁, b at x₂) with x₁ < x₂ to (b at −x₁, a at −x₂), i.e. back to the ordered pair |ab⟩, so it does not mix the two channels. Θ keeps valleys and inverts space, so it maps |ab⟩ ↔ |ba⟩ and forces equal diagonal entries. So 𝒜 = [[A_d, A_o], [A_o*, A_d]] on {|ab⟩, |ba⟩}. Its off-diagonal element is the q → 0 transmission (the valley labels exchange positions), and its eigenvalues a_± = A_d ± |A_o| are the excluded lengths of the two inter-valley channels.

Writing P_{↑↑} + P_{↓↓} = 1/2 + 2S^zS^z, P_{↑↓} + P_{↓↑} = 1/2 − 2S^zS^z and |↑↓⟩⟨↓↑| + h.c. = S⁺S⁻ + S⁻S⁺, (4.1) becomes the pseudospin XXZ chain on the particle index,

    H_σ = Σ_i [ J_⊥ (S^x_i S^x_{i+1} + S^y_i S^y_{i+1}) + J_z S^z_i S^z_{i+1} ],
    J_⊥ = (2π²n³/3m) |A_o|,  J_z = (2π²n³/3m)(a_s − A_d),  Δ_σ = (a_s − A_d)/|A_o|.   (4.2)

(The phase of A_o and the sign of J_⊥ can be removed by a rotation of alternate pseudospins.) S^z_tot = (N_a − N_b)/2 is conserved, because a valley flip carries momentum k_b − k_a, and Σ acts as S^z → −S^z.

Consistency check [derived]. If the valleys were an SU(2) doublet (a_s = a_+, the odd channel), (4.2) gives Δ_σ = ±1, an isotropic Heisenberg chain. For the δ-interacting two-component gas, a_odd = 0 and a_even = −2/(mg), so J = (2π²n³/3m)(2/(mg)) = 4π²n³/(3m²g), antiferromagnetic. This is the standard strong-coupling Yang–Gaudin result (Ogata–Shiba factorization), recalled, not re-derived from a source here.

### 4.2 From the pseudospin chain to K_−

[derived] The pseudospin chain (4.2) at zero magnetization has Luttinger parameter K_σ = π/[2(π − arccos Δ_σ)] (K_σ = 1/2 at Δ_σ = 1, 1 at Δ_σ = 0, → ∞ as Δ_σ → −1). Its umklapp cos 4φ_σ has dimension 4K_σ, and in the full system it is the inter-valley backscattering cos 2√2 φ_− of dimension 2K_−, so K_− = 2K_σ. The charge sector is the spinless Fermi sea of total density n, so K_+ → 1/2 as n → 0 (the U = ∞ Hubbard value in this normalization). With excluded lengths of order the lattice spacing, K_+ ≈ (1 − n a_eff)²/2. For blockaded rods (a_eff = 1) at n = 0.17 that gives K_+ ≈ 0.34 [heuristic]. The velocities are u_+ ≈ πn/m and u_− ≈ v_σ/n, with v_σ the pseudospin spin-wave velocity in units of particle index, so u_− is of order n²a/m ≪ u_+.

### 4.3 Criterion

[derived] In the dilute limit, the Σ-point gas is

- c = 2 (two-component liquid) for −1 < Δ_σ < 1, i.e. |a_s − A_d| < |A_o|, which is a_− < a_s < a_+;
- marginal at Δ_σ = 1 (isotropic antiferromagnetic pseudospin, still c = 2 with logarithmic corrections in the pure XXZ chain, but the bare backscattering coupling of the full model decides the sign of the marginal flow);
- c = 1 vector chiral (valley ferromagnet, one valley filled, Σ broken) for Δ_σ < −1, i.e. a_s < a_−;
- c = 1 valley-crystal liquid (pseudospin gapped, section 5) for Δ_σ > 1, i.e. a_s > a_+.

[derived] Relation to Kolezhuk, Heidrich-Meisner, Greschner and Vekua (arXiv:1112.5720). They write the dilute two-valley magnon gas as a two-component Bose gas with couplings g_ab = −2/(m a_ab) extracted from the exact two-magnon problem, and place the vector-chiral/TLL2 boundary at g_11 = g_12, which in the strongly coupled regime they interpret as a crossing of scattering lengths. Their single inter-species length a_12, from the phase shift of the mixed k, k̃ scattering state (their eqs. 28 to 31), corresponds in (4.2) to the inter-valley channel that competes with the same-valley one, so their boundary is Δ_σ = −1, a_s = a_−. The Gross–Pitaevskii description has no Δ_σ > 1 side, so the valley-crystal instability is not in their criterion. I identify their a_12 with a_− by this argument only [heuristic identification]. Their statement that at low density both species form a single Fermi sea, so that the vector chiral to TLL2 transition shows no magnetization cusp, is the impenetrable picture of 4.1 [literature].

[derived] The drift (S1) in the dilute limit. λ ≈ ε‴(πn)²/2 is of order n², while u_+K_+ ~ n and u_−/K_− ~ n². So λ² ~ n⁴ against n³, and the drift instability is subleading as n → 0. At n = 0.17, |w| = 1, the estimates λ ≈ 0.83, u_+K_+ ≈ 1.4, and an inverse susceptibility u_−/K_− of order 0.6 (from the ring ED gaps) give λ²/(u_+K_+) ≈ 0.5, the same order as u_−/K_−. So at the candidate's density the drift correction to the valley susceptibility is not small [heuristic, order-of-magnitude only]. A direct measurement, the energy of the valley-imbalanced ring states at their optimal momentum (C178), includes it automatically.

## 5. The competing c = 1 phases

[derived] All of them are one of the following.

1. Vector chiral (valley-polarized liquid). The pseudospin is ferromagnetic (Δ_σ < −1), or the Schur complement of (S1) vanishes. One valley fills, Σ breaks spontaneously (twofold degenerate), and the remaining boson has c = 1. Its order parameter is the valley imbalance, equivalently the spin current κ = ⟨(S_j × S_{j+1})_z⟩ of the J1–J2 picture, with ⟨∂θ_+⟩ ≠ 0. In the dilute limit the transition at Δ_σ = −1 has a quadratic pseudospin dispersion (z = 2 in the pseudospin sector) and no cusp in the density (KHGV).
2. Valley-crystal liquid (pseudospin gapped). cos 2√2 φ_− is relevant (K_− < 1, Δ_σ > 1) and locks φ_−. The ground state is unique: the two minima of the cosine in a period π/√2 of φ_− are identified with each other by the compactification lattice (π m_a, π m_b) of (φ_a, φ_b) together with a shift of the gapless φ_+. In pseudospin language this is the Néel state of the label chain. The labels alternate along the particle sequence, and the "two Néel states" differ by moving one particle, which the charge sector absorbs. The dominant correlations are the valley-density wave at 2k_F per valley, e^{i√2φ_+} sin √2φ_−, with exponent K_+/2, and the inter-valley pair ψ_aψ_b ~ e^{i√2θ_+} e^{∓i√2φ_−} at momentum k_a + k_b = π, with exponent 1/(2K_+). For K_+ < 1 (repulsive) it reads as a valley crystal, and for K_+ > 1 as an ab pair liquid. These are the same phase with different dominant correlations, the analogue of the Luther–Emery liquid. Whether a gapless-SPT distinction (in the sense of Keselman and Berg, arXiv:1502.02037, abstract only, a spin-gapped c = 1 liquid with protected edge modes) separates the two regimes was not examined.
3. Bound-state liquids. If a two-body bound state lies below 2ε_min, the dilute gas is a gas of pairs. An ab bound state at total momentum π gives phase 2 directly. An aa (and by Σ, bb) bound state at total momentum 2k_a gives a gas of pairs with its own two valleys at 2k_a and 2k_b = 2π − 2k_a, related by Σ. It repeats the whole question one level up: a pair vector-chiral phase (c = 1) or a two-component pair liquid (c = 2). This is the analogue of the nematic and multipolar phases of the ferromagnetic J1–J2 chain near saturation (Hikihara, Kecke, Momoi and Furusaki, arXiv:0807.0858, abstract only) and of the pair superfluid of Mishra, Greschner and Santos (arXiv:1411.7270). Note that in the c gauge each leg of (0.1) is a ferromagnetic Heisenberg chain at w = 0, whose two-magnon bound states exist at every nonzero leg momentum. C173's twobody.py found no bound state below threshold at φ = 1.3, but the lowest two-body level sat only 0.002 to 0.008 above threshold, which may signal a near-zero-energy resonance and a large |a|. At φ = π/2 it was not checked [observation from C173's PROGRESS].
4. Charge-locked phases. At commensurate densities with relevant charge locking (n = 1/3 needs K_+ < 1/9 from the p = 6 vertex, the spinless period-3 condition K < 2/9 in the usual normalization), φ_+ locks. The pseudospin can stay critical (c = 1) or gap (c = 0). C173's period-9 gapped states near n ≈ 0.3 and the first-order jump into the CDW at small |w| belong here.
5. ΔN = 4 locking (parity class only, K_+ > 1). θ_+ locks, a "quartet" condensate of ab pairs with the pseudospin left critical, c = 1. It needs strong attraction and is not expected in the dilute repulsive regime.

## 6. Consistency checks

- [derived] Free-fermion limit. Without the generated range-2 hop, a nearest-neighbour qubit chain has one harmonic in its band and at most two Fermi points (C169, c ≤ r at range r). The cusp evades this only because the interaction (the blockade) kills the range-1 amplitude to leading order while generating range 2 at order Δ²/V. The field theory above is the r = 2 four-Fermi-point liquid (Titvinidze–Japaridze type, C156 2b.3) with interactions.
- [derived] The |w| = 0 resonance reproduces C153 (section 2.2).
- [literature check] Verresen, Vishwanath and Pollmann (arXiv:1903.09179, appendix C) give the second-order Hamiltonian of the Kitaev chain −Σ X_iX_{i+1} + gZ with a soft constraint U n_i n_{i+1}, in C173's convention w = Δ = −1. Their next-nearest hop has coefficient −1/(2U), which matches C173's t₂ = −Δ²/(2V) with V = U. Their n_i n_{i+3} coefficient −1/(3U) matches V₃ = −Δ²/(3V). Their n_i n_{i+2} coefficient −3/U equals C173's −Δ²/V plus −2t²/U from virtual first-order hops next to an occupied site, which are absent at the cusp (t → 0). Their U(1)-breaking term is proportional to tΔ/U, which matches C173's statement that ΔN = 2 enters at order |w|Δ/V. (The text extraction of their eq. A2 leaves the placement of a prefactor 1/2 slightly ambiguous. The reading above is the one consistent with their main text.)
- [literature check] In the J1–J2 language of (0.2) with the blockade but without V₂ and V₃, the model is the hard-core frustrated triangular ladder of Mishra, Greschner and Santos at t = |t′| (their H with t = 1, t′ = −1, after d_j → (−1)^j d_j) at V → ∞ (section 7.1).
- [numerical] optable.py: band facts, operator enumeration, dimension formula (optable.out).

## 7. Literature

### 7.1 The J1–J2 chain near saturation, frustrated ladders, and the closest calibration

- [literature] Kolezhuk, Heidrich-Meisner, Greschner and Vekua, PRB 85, 064420 (2012), arXiv:1112.5720 (read). Dilute two-valley magnon gas, effective couplings from the exact two-magnon problem via a Bethe–Salpeter ladder with an infrared cutoff or off-shell energy E* = π²μ/8, criterion Γ₁₂ < Γ₁₁ for TLL2 (c = 2) against vector chiral. For isotropic antiferromagnetic spin-1/2 chains at saturation they predict vector chiral for 0 < J1/J2 < 1.682 and TLL2 for 1.682 < J1/J2 < 4 (their Table I), and they report that this location disagrees "unexpectedly" with DMRG for spin 1/2 (Okunishi, Hikihara et al.), while the topology agrees. For anisotropic ferromagnetic spin-1/2 chains (J1z/J1 = J2z/J2 = Δ ∈ [0, 1]) the TLL2/VC/bound-state boundaries agree with DMRG. Point relevant here: at S = 1/2 the method works even when |g|m ≫ 1, but its quantitative accuracy on the antiferromagnetic side is poor.
- [literature] Hikihara, Momoi, Furusaki and Kawamura, PRB 81, 224433 (2010), arXiv:1004.0622 (read). Antiferromagnetic J1–J2 chain in a field: TLL1, TLL2 (c = 1 + 1, two Fermi seas of Jordan–Wigner fermions), SDW₂ and vector chiral gapless phases, plus dimer, 1/3 plateau and saturated phases. The high-field TLL2 is equal-density condensation of the ±Q₀ magnons, and unequal densities give the vector chiral phase. The vector chiral phase occurs in two narrow regions and is always adjacent to a TLL2 phase. Their g_1⊥ backscattering, eq. (60), is the operator of section 2.1.
- [literature] Okunishi and Tonegawa, JPSJ 72, 479 (2003), cond-mat/0211458 (abstract only). Strongly frustrated antiferromagnetic zigzag chain in a field: 1/3 plateau, cusps, even–odd effect.
- [literature] Hikihara, Kecke, Momoi and Furusaki, PRB 78, 144404 (2008), arXiv:0807.0858 (abstract only). Ferromagnetic J1: vector chiral at weak field, nematic and higher multipolar (bound-magnon) phases near saturation.
- [literature] Mishra, Greschner and Santos, PRA 91, 043614 (2015), arXiv:1411.7270 (read). Hard-core bosons on the frustrated triangular ladder, H = −t Σ (a†_i a_{i+1} + h.c.) − t′ Σ (a†_i a_{i+2} + h.c.) + V Σ n_i n_{i+1} + W Σ n_i n_{i+1} n_{i+2}, t = 1, t′ < 0. This is the closest published model to (0.2): with V → ∞ and t′ = −1 it is (0.2) at |w| = 1 without V₂ and V₃. Their findings:
  - Dilute limit (their section V, using KHGV): mapping to a J1–J2 spin model with J1z = V, J2z = 0. For J2/J1 > 1/√8 increasing J1z drives chiral superfluid (one valley) to two-component superfluid (c = 2). At J2 = J1 the transition is at J1z ≈ 1.8 J1, i.e. V ≈ 3.6t. The blockade V = ∞ is deep on the c = 2 side.
  - Strong coupling (their section VI): at V = W ≫ t, |t′| and ρ < 1/2 they squeeze out the blockade and obtain a J1–J2 model with correlated next-nearest hopping (their model 8). Its phase diagram has SF and 2SF, and "We do not find any region with finite chirality".
  - DMRG at V = W = 8, t = |t′| = 1, L = 180: c = 2 (2SF) for ρ < ρ₁ ≈ 0.4, then SF (c = 1). They warn that at low fillings finite-size corrections are strong and the open-chain c fits drift upward.
  So the blockaded frustrated ladder at the candidate's ratio t/|t′| = 1 and density is a published c = 2 phase, with a range-2 hop put in by hand.
- [literature] Greschner, Santos and Vekua, PRA 87, 033609 (2013), arXiv:1202.5386 (abstract only): bosons on zigzag lattices, chiral superfluidity, three-body constraint. Mishra, Pai and Mukerjee, PRA 89, 013615 (2014), arXiv:1309.4412 (abstract only): hard-core bosons with frustrated next-nearest hopping and nearest-neighbour V, bond order, CDW and supersolid, no central charge given. Rossini, Lante, Parola and Becca, PRB 83, 155106 (2011), arXiv:1101.3217 (search summary only): unfrustrated zigzag hard-core bosons, one-component TLL at ρ < 1/3.

What this calibration says [derived from the literature]. The question about the candidate narrows to the generated diagonal terms. In the J1–J2 spin language, (0.2) differs from the blockaded MGS model only by J2^z = −J2 (the attraction −2 n_j n_{j+2}, which makes the c-gauge legs ferromagnetic Heisenberg chains) and the third-neighbour attraction −(2/3) n_j n_{j+3}. Neither KHGV's anisotropic diagram (J2z/J2 ∈ [0, 1]) nor MGS (J2z = 0) covers J2z/J2 = −1. [heuristic] A Born-level estimate of g_ab − g_aa ∝ V̂(2k₀) − V̂(0), which is not valid in the hard-core regime and is offered only as a sign indicator, gives at |w| = 1 (cos 2k₀ = −0.875, cos 4k₀ = 0.531, cos 6k₀ = −0.057 in the d gauge): the nearest-neighbour repulsion contributes 2V(−0.875 − 1) < 0 (toward c = 2), the r = 2 attraction contributes −4(0.531 − 1) = +1.9 > 0, and the r = 3 attraction contributes −(4/3)(−0.057 − 1) = +1.4 > 0 (both toward vector chiral). So the generated attractions push in the dangerous direction, and the sign is not settled by any published calculation. The exact two-body computation (C175) decides it. It should reproduce, as a check, the MGS/KHGV prediction that the V₂ = V₃ = 0 limit is two-component.

### 7.2 Nearest-neighbour qubit chains inside B₂ near this corner

- [literature] Verresen, Vishwanath and Pollmann, arXiv:1903.09179 (read). Kitaev chain with a soft or hard Rydberg constraint, a B₂ chain (−Σ XX + gZ + U n n). A stable Luttinger liquid with emergent U(1), c = 1, K = (1 − ρ₀)², critical for −2 ≤ g ≲ 0.0697 at U = ∞. Here the first-order hop equals the pairing (t = Δ), so it sits in C173's strong-rung regime (Lemma S). They do not consider t ≲ Δ²/U, and their hopping is real. Their conclusion lists "central charges c > 1" only as a future direction.
- [literature] Chepiga and Mila, PRB 107, L081106 (2023), arXiv:2206.11754 (read in part). Interacting Kitaev chain with nearest-neighbour pairing λ and repulsion V (a B₂ chain, XYZ in a z field): Z₂, period-2 and floating (c = 1) phases, an eight-vertex multicritical line, c ≈ 1 everywhere it is measured. Real hopping only. With real hopping the two band minima of (0.1) at φ = 0 sit at 0 and π with energies differing by 4|w|, so the Σ degeneracy needs the DM term, which this work does not include.
- [literature, C156] XYZ in a field, XYZ with DM (Peotta et al., arXiv:1403.4568), XXZ with DM and transverse field, Kitaev–Gamma family: all c ≤ 1, none at the cusp. Not repeated here.
- [literature] Corcoran, de Leeuw and Pozsgay, arXiv:2405.15848 (read in part). Classification of time- and space-reflection-symmetric integrable Rydberg-constrained Hamiltonians: at range 3 only the RSOS chains, at range 4 the constrained XXZ chain and a new "double golden chain", and a partial range-5 classification. [derived] In the d gauge (0.2) is real and reflection symmetric. Its next-nearest hop needs the sites on both sides of the hop to be empty, so its density acts on five sites (range 5 in their counting). The range-5 integrable models they list have either a next-nearest hop alone (their H_A, with Hilbert-space fragmentation) or nearest-neighbour hops (H_B, H_C), never both together with the attractions of (0.2). So no integrable or exactly solvable point is known for the candidate. Their range-5 classification is partial, so this is not a proof of non-integrability.
- [literature] Lesanovsky and Katsura, PRA 86, 041601(R) (2012), arXiv:1204.0903 (abstract only): Rydberg lattice gas mapped to interacting Fibonacci anyons (golden-chain physics, c < 1). Chen and Iadecola, PRB 103, 214304 (2021), arXiv:2103.09773 (abstract only): PXP in a strong staggered field, second-order next-nearest hopping in the Fibonacci space, dynamics, two-site cell, not B₂. Sohal, Bernien and Verresen, arXiv:2610.40771 (abstract only): blockaded PXYP ladders, c = 1 gapless topological phase for two legs, a two-component (c = 2) liquid for three legs. Ladders are not B₂.

### 7.3 Is a c = 2 phase of a translation-invariant nearest-neighbour spin-1/2 chain claimed or excluded?

[literature, with the C156 survey] Not claimed anywhere I could find. Every published c = 2 two-component phase in a spin-1/2 chain (J1–J2 TLL2, the frustrated hard-core ladders, Titvinidze–Japaridze, Sudan–Lüscher–Läuchli, Fitouchi–Chepiga) has range 2 or more, or is a ladder. The nearest-neighbour studies closest to the corner (Verresen–Vishwanath–Pollmann, Chepiga–Mila, and C156's list) report c ≤ 1, but none enters the cusp (first-order hop ≲ Δ²/V with the DM term making the valleys degenerate). Not excluded either. No theorem covers it: the rigorous Luttinger-liquid constructions (Benfatto–Falco–Mastropietro) are near free fermions with two Fermi points, frustration-free results need frustration-freeness, and the cusp sits near the extensively degenerate Fibonacci point, outside all of them. Latorre and Sierra state the expectation that no spin-1/2 chain beats the Heisenberg value c = 1 (arXiv:2403.17242, discussed next), but they give it as a conjecture.

## 8. Relation to Latorre and Sierra's conjecture at d = 2

[literature] Latorre and Sierra (arXiv:2403.17242, read) state the conjecture for H_L = Σ h_{i,i+1} with periodic boundary conditions and one-site translation invariance. If the low-energy spectrum is "a collection of unitary CFTs with central charges c_a", then c_T = Σ_a c_a ≤ d − 1. They explicitly allow different velocities v_a (they mention the Hubbard chain's spin and charge velocities). At d = 2 they note that the antiferromagnetic Heisenberg chain saturates the bound and add that, by the conjecture, no other spin model can give a larger central charge.

[derived] Where the candidate sits. Its bond is a translation-invariant nearest-neighbour qubit bond. If the two-component phase is confirmed, its infrared is two Gaussian bosons with u_+ ≠ u_−, coupled by the exactly marginal drift terms of (1.5). That is a collection of unitary CFTs in their sense, with c_T = 1 + 1 = 2 = d. It would violate c ≤ d − 1 at d = 2 by one unit and would saturate c ≤ d. Their tensor-product consistency argument (f(d₁) + f(d₂) ≤ f(d₁d₂)) is untouched. Their remark that no spin-1/2 chain exceeds Heisenberg's c = 1 would be false. The project's replacement bound A = Σ c_a v_a ≤ π M log₂ d is not threatened: A = u_+ + u_− is of order n|t₂| here, while M is of order J_z.

## 9. Open

- The sign of Δ_σ, i.e. (a_s, A_d, |A_o|) from the exact two-body problem at φ = π/2 (C175). Calibration: the V₂ = V₃ = 0 limit should come out two-component (MGS, KHGV).
- Whether aa or ab two-body bound states exist at φ = π/2, and whether the near-threshold level seen by C173 at φ = 1.3 is a resonance (large |a| would break the dilute expansion at n = 0.17).
- The size of the drift correction to the valley susceptibility at n = 0.17 (section 4.3), best measured directly on rings.
- The constants in (1.4), which fix the sign of cos 2√2φ_−. This matters only at the marginal point Δ_σ = 1.
- Whether the pseudospin-gapped phase has a gapless-SPT refinement.
