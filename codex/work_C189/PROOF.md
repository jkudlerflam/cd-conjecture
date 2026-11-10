This was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

# C189 derivations: Luttinger parameters of the two-valley liquid and the infinite-chain test

Standard: theoretical physics (bosonization, boundary CFT, calibrated numerics). Tags [proved], [derived], [numerical], [heuristic] open the paragraph they apply to. Every number names its file in this folder.

## 1. Model, frames and the valley field

C173's effective model (PROOF eq. 3.2), units |t₂| = 1, open or infinite chain, blockade n_j n_{j+1} = 0:

    H = P Σ_j [ (w c†_{j+1} c_j + h.c.) + (t₂ c†_{j+2} c_j + h.c.) + V₂ n_j n_{j+2} + V₃ n_j n_{j+3} ] P,   (1.1)

w = u e^{iφ}, t₂ = −1, V₂ = −2, V₃ = −2/3. [proved] Particles never pass. The NNN hop from j to j+2 needs j+1 empty (blockade), so the Jordan–Wigner string of the hop is +1 and (1.1) is the same matrix in the hard-core boson basis (b_j) and in the fermion basis (c_j = (−1)^{N_{<j}} b_j). Statistics is a choice of fields, not of physics.

Band (c frame). With c_j = L^{-1/2} Σ_k e^{ikj} c_k, Σ_j c†_{j+1}c_j = Σ_k e^{−ik} c†_k c_k, so ε(k) = 2u cos(k − φ) − 2cos 2k. At φ = π/2 the minima sit at sin k = −u/4, i.e. k_A = −0.2527 and k_B = π + 0.2527 at u = 1, with ε_min = −9/4. In the d frame d_j = (−i)^j c_j the momenta shift by −π/2, the valleys sit at ±k₀, cos k₀ = −1/4, and H is real at φ = π/2.

[proved, exact identity] The valley field. Put φ = π/2 − ϑ. Then cos(k − φ) = sin k cos ϑ + cos k sin ϑ, so

    H(ϑ) = H(0) + u sin ϑ P_tot + u (cos ϑ − 1) Q_tot,   P_tot = Σ_k 2 cos k n_k = Σ_j (c†_{j+1}c_j + h.c.),   (1.2)

with Q_tot = Σ_k 2 sin k n_k. P_tot is C176's valley polarization operator (the vector chirality). Near the valleys 2cos k ≈ ±2cos k_A = ±1.94, so P_tot is a pseudospin field of strength δ = 4u cos k_A ϑ (C181 F3), coupling to N_a − N_b. This gives the susceptibility measurement and the operator identification of C_P used in section 4.

## 2. What a hard wall does to the valleys

[derived; my first version was wrong, corrected by the ff2v control] With NN hopping only, open-chain eigenstates are sin(kx), so a wall reflects k into −k. In the d frame that maps valley + into valley −. With NN and NNN hopping a single sine cannot vanish at the two virtual sites x = −1 and −2. At energy E the bulk equation has four real roots ±k₁, ±k₂ (k₁ < k₀ < k₂) plus evanescent ones, and the wall's 2×2 reflection matrix sends the incoming pair (+k₁, −k₂) to the outgoing pair (+k₂, −k₁) with both intra-valley (k₁ → k₂) and inter-valley (k₁ → −k₁) amplitudes. Which one dominates is not fixed by symmetry.

[numerical, kest_controls.out] For free fermions with the candidate's band (ff2v, exact Slater determinant, L = 1024) the wall oscillation is at 2k_F = πn (amplitude 0.24), with 0.03 at 2k₀ and nothing at 2(k₀ ± k_F). The free wall is mostly valley preserving.

[derived] Boundary RG. Write the two-component boson fields so that a valley-preserving wall pins φ_c and φ_s and a valley-flipping wall pins φ_c and θ_s. At the φ_s-pinned boundary the valley-flip reflection is the boundary operator e^{i√2θ_s}, of boundary dimension 1/K_s. At the θ_s-pinned boundary the valley-preserving reflection is e^{i√2φ_s}, of dimension K_s. So for K_s > 1 the valley-flipping boundary is the stable one and for K_s < 1 the valley-preserving one. This matters for the Friedel oscillations and for the zero modes of the parity gap (4.5).

## 3. Two bosonizations and which one the data select

Normalization (Giamarchi). For each valley species σ = a, b: ρ_σ = −∂φ_σ/π + density harmonics e^{2ip(πρ_σ x − φ_σ)}, [φ_σ(x), ∂θ_σ(y)] = iπδ. Then φ_{c,s} = (φ_a ± φ_b)/√2, likewise θ, and H = Σ_{ν=c,s} (v_ν/2π) ∫ [K_ν (∂θ_ν)² + K_ν^{-1}(∂φ_ν)²]. For a vertex operator e^{i(αφ_ν + βθ_ν)} the two-point exponent is (α²K_ν + β²/K_ν)/2. [check] e^{2iφ} gives 2K and e^{iθ} gives 1/(2K), the spinless Luttinger values. Free spinless fermions and hard-core bosons both have K = 1.

### 3.1 Fermion frame (registered, superseded)

Two Fermi seas of one fermion species. ψ_{rσ} ∝ e^{i(θ_σ − rφ_σ)} (r = ±1 for R, L), with fermionic Klein factors between σ = a and b. [derived]
- One-body ⟨c†c⟩: e^{i(θ_c − rφ_c)/√2} e^{i(θ_s − rφ_s)/√2}, exponent (K_c + 1/K_c + K_s + 1/K_s)/4.
- Every charge-1 vertex operator allowed on this lattice of charges has |α_ν| = |β_ν| = 1/√2 in both sectors, so its exponent is at least (x + 1/x)/4 summed over ν, which is ≥ 1. [proved, AM–GM]
- Inter-valley pair at total momentum 0 (d frame), ψ_{R a}ψ_{L b}: e^{i√2θ_c} e^{±i√2φ_s}, exponent K_s + 1/K_c.
- Same-valley pair: e^{i√2θ_c} e^{i√2θ_s}, exponent 1/K_c + 1/K_s.
- The lattice boson b = c e^{iπN} multiplies the charge part by e^{±i√2φ_c}, which only flips the sign of the φ_c coefficient. Its exponent equals the fermion's.

### 3.2 Boson frame (amendment A3)

Two quasi-condensates of hard-core bosons, b ≈ e^{ik_a x} b_a + e^{ik_b x} b_b, b_σ ∝ e^{iθ_σ}, mutually bosonic. [derived]
- Boson one-body at a valley momentum: e^{iθ_a} = e^{i(θ_c + θ_s)/√2}, exponent η_B = 1/(4K_c) + 1/(4K_s).
- Inter-valley pair b_a b_b (total momentum k_a + k_b, π in the c frame): e^{i√2θ_c}, exponent η_P = 1/K_c, with no K_s term. Under the blockade c_j c_{j+d} = −b_j b_{j+d} for d = 2, 3, because the sites between are empty. So C176's D(r) is this boson pair.
- Same-valley pair b_a b_a: e^{i√2(θ_c + θ_s)}, exponent 1/K_c + 1/K_s.
- Jordan–Wigner fermion. c = (−1)^{N_{<j}} b = b_a e^{iπN_a} e^{iπN_b} + (a ↔ b). The own-valley string turns b_a into a fermion e^{i(θ_a − φ_a)}, and the other valley's string is e^{−iφ_b}. Together: e^{i[(θ_c + θ_s)/√2 − √2 φ_c]}, exponent K_c + 1/(4K_c) + 1/(4K_s). [check] For c2 (K = 1) this is 3/2. That is the leg-fermion exponent 1 plus the other leg's string 1/2, the value I derived directly for two decoupled XX legs.
- Inter-valley coherence b_b† b_a at momentum k_a − k_b (2k₀ in the d frame): e^{−i√2θ_s}, exponent 1/K_s.
- Density correlations are the same in both frames: ⟨ρ_ν ρ_ν⟩ = −K_ν/(π² r²) for ν = c, s.

### 3.3 Which frame [numerical, kest_boson.out, kest_run.out]

- The lattice boson decays with η_B = 0.32 to 0.34 in the candidate, at L = 128 to 256 and all three densities. This is impossible in the fermion frame, where every charge-1 exponent is ≥ 1 (3.1). The same fit reads 0.504 on c2 (exact 1/2).
- With η_P = 0.62 the fermion-frame relation gives K_s = η_P − 1/K_c < 0 for any K_c ≤ 1.6. The registered fermion-frame estimator returned −0.21 ± 0.18.
- The boson frame is consistent. η_P = 1/K_c gives K_c = 1.61, the long-wavelength S(q) gives 1.5 to 1.65, η_B gives K_s ≈ 1.5, the predicted fermion exponent K_c + 1/(4K_c) + 1/(4K_s) = 1.93 matches the measured 1.9 to 2.0, and the 2k₀ coherence decays slowly (exponent 1/K_s ≈ 0.7) as the wall profile shows.

### 3.4 The backscattering operator and its relevance

[derived] The non-oscillating inter-valley terms come from products of density harmonics. V ρ_a ρ_b contains cos(2φ_a − 2πρ_a x) cos(2φ_b − 2πρ_b x), whose difference harmonic is cos(2(φ_a − φ_b) − 2π(ρ_a − ρ_b)x). At valley balance (φ = π/2) this is cos(√8 φ_s), with scaling dimension (8K_s)/4 = 2K_s. The operator is a density product, so its form and dimension are the same in both frames. The θ-dual, cos(√8 θ_s), would move two particles between valleys and carries momentum 2(k_a − k_b) = 4k₀ = 1.011 (mod 2π), so translation symmetry forbids it at incommensurate 4k₀/2π. [proved for the momentum count]

[derived, standard] Kosterlitz–Thouless flow for y = g₁/(πv_s): dK_s/dl = −y² K_s²/2, dy/dl = 2(1 − K_s) y. Trajectories are hyperbolas, (K_s − 1)² − y²/4 ≈ const near (1, 0).
- K_s − 1 > |y|/2: y → 0 and K_s → K_s* > 1. Gapless, c = 2.
- Otherwise y grows and a pseudospin gap opens, at Δ_s ∝ exp(−π/√(y² − 4(K_s−1)²)) near the separatrix (exponentially small).
- Signatures in effective, scale-dependent estimators: a marginally irrelevant flow approaches 1 from above like 1/(2 ln(L/L₀)), a relevant one falls through 1.

[heuristic] The fermion intuition does not carry over. There, two Fermi seas start at K_s = 1 and any attraction pushes K_s below 1, which gives an exponentially small gap. In the boson frame the reference point is two quasi-condensates, with K_s ∝ (g − g_ab)^{−1/2} large for weak coupling. An inter-valley attraction g_ab < 0 lowers K_s, and the pair superfluid (K_s < 1, cos(√8φ_s) relevant) needs the attraction to beat the condensates' stiffness. The dilute two-body bound pair of C175 is that regime. Whether the candidate densities are past it is a quantitative question, and the measured K_s ≈ 1.3 to 1.5 says they are not.

## 4. Estimators

### 4.1 K_c from S(q) [derived]

The total density is conserved, so its long-wavelength part is exactly −(√2/π)∂φ_c with no renormalization (Ward identity). ⟨ρ_c ρ_c⟩ = −K_c/(π² r²) gives S(q) = K_c|q|/π for |q| below the lowest oscillation wavevector, and K|q|/(2π) for one component. [checks] Free two-component fermions give |q|/π (K_c = 1). Two XX legs give 2·(L/2)·S_leg(2q)/L = |q|/π, also K_c = 1. In kest.py, S(q) = C(0) + 2 Σ_{r ≤ L/2−4} C(r) cos qr, with C(0) = n(1 − n) and C(r) averaged over the four central rows.

### 4.2 Pair and boson exponents (boson frame)

K_c^P = 1/η_P and 1/K_s^B = 4η_B − η_P (3.2). Errors: σ(1/K_s) = √((4σ_B)² + σ_P²).

### 4.3 Friedel oscillations

For a φ-pinned wall, ⟨e^{i(αφ_c + βφ_s)}(x)⟩ ∝ d(x)^{−(α²K_c + β²K_s)/4}, with d(x) = (2L'/π) sin(πx/L') and L' = L + 1. This is the conformal map of the half-plane one-point function. The 2k_F component e^{i√2(φ_c ± φ_s)} gives (K_c + K_s)/2 (the registered E_F'). At a θ_s-pinned wall the 2(k₀ ± k_F) components e^{i√2φ_c} e^{±i√2θ_s} give (K_c + 1/K_s)/2 (the registered E_F). In the candidate neither is clean: 2k_F is short ranged and the dominant 2k₀ modulation is the coherence e^{i√2θ_s}, which only a valley-flipping wall can seed. That is the stable wall for K_s > 1 (section 2), so this is qualitatively consistent with K_s > 1. I did not extract a number from it.

### 4.4 Pseudospin density: K_s^P, free of the frame and of Z [derived]

By (1.2), P_tot couples to the valley imbalance. Its low-energy projection is P(x) = Z_P ρ_s(x) + (oscillating terms) + (current terms, dropped). The bare weights 2cos k at a valley's two Fermi points are 2.00 and 1.70 (k_F ≈ 0.3). Their mean, 1.85, is Z_P, and their half-difference, 0.15, multiplies a current. That is 8 percent in amplitude and under 1 percent in ⟨PP⟩. Two consequences:
- ⟨P(x)P(0)⟩_smooth = −Z_P² K_s/(π² r²).
- ⟨P⟩/L = Z_P (n_a − n_b) in a polarized state.

The valley weights of n(k) give n_a − n_b = p n, up to a factor Z_N that is 1 for free fermions. So Z_P = ⟨P⟩/(p n) is measured directly on the off-line runs, and K_s^P = π² A/Z_P² with A = −r² C_P(r) averaged over L/8 ≤ r ≤ L/4. [check] For ff2v, Z_P = 1.86 (bare value 2cos k_A = 1.94, reduced by the k-dependence of 2cos k across the Fermi seas), and K_s^P reads 0.94 to 1.04 against exact 1 (kest_cp.out). For the candidate, Z_P = 1.864 to 1.872 over φ = 1.50 to 1.54 and L = 128 to 256, a stable operator weight.

### 4.5 Parity gap, compressibility, susceptibility [derived]

Zero-mode energies of a wall that pins φ_ν (Dirichlet), with φ_c(L) − φ_c(0) = −πN_c/√2:

    E = (π/4L) [ (v_c/K_c) ΔN_c² + (v_s/K_s) ΔN_s² ].

If the wall pins θ_s instead, the second term becomes v_s K_s ΔJ_s², where J_s = N_out − N_in has the parity of N_c. [check] For one free species the zero-mode energy is (πv/2L)ΔN², so E(N+1) + E(N−1) − 2E(N) = πv/L, the open-chain level spacing. With two species and one particle added to species a, ΔN_c = ΔN_s = 1 and the formula above gives πv/(2L) per step, the same thing.

Adding one particle changes N_c by 1 and N_s (or J_s) by ±1, so L·Δ_P = (π/2)(v_c/K_c + v_s/K_s), or (π/2)(v_c/K_c + v_s K_s) for the θ_s wall.

Compressibility: ∂n/∂μ = 2K_c/(πv_c), so v_c/K_c = (2/π) ∂²e/∂n², taken from the second difference of E(N) over n = 1/8, 5/32, 3/16 at fixed L.

Susceptibility: a field μ_s = δ/2 coupling to N_a − N_b shifts ∂φ_s and gives n_a − n_b = (2K_s/πv_s) μ_s, so χ_s = 2pn/δ = 2K_s/(πv_s).
- Valley-preserving wall: the parity gap and χ_s both measure v_s/K_s, so they only test consistency.
- θ_s wall: their product gives K_s².
[numerical, open] The test fails by about 30 percent (RESULT Part 1). Finite-size quantization of N_s at a valley-preserving wall gives a threshold δ_c = π(v_s/K_s)/L, which explains the vanishing p at δ = 0.022 for L = 64 and 128 without a gap.

## 5. Open-chain fit forms and calibration

[derived] For a bulk vertex operator in a strip with a conformal wall, the two-point function is

    d(x−y)^{−η} [ d(x+y)² / (d(2x) d(2y)) ]^{γ},

with γ fixed by how the image charges of φ and θ enter. I fit η and γ freely, with x the fixed reference site and y = x + r. The envelope is the centred RMS of |G| over ±4 samples, which is unbiased to first order for a monotone power law. A one-sided running maximum shifts a monotone function by half its window, which is how the first version of these fits misread c2 as 0.8.

[numerical, kest_boson.out] Calibration at L = 256:
- c2: boson 0.504 (exact 0.5), inter-leg pair 1.020 (1), intra-leg pair 2.041 (2).
- c1: boson 0.918 (exact 1/(2K) = 0.889).
The registered fermion-frame fits (chord only, one-sided envelope) read ff2v's fermion exponent as 1.15 to 1.47 against 1, which is why E_G was excluded.

## 6. The infinite chain

### 6.1 Labels [derived]

Density n = p/q, q-site unit cell, two cells per iDMRG window. At a bond at global position X, a left label is (a, f) with a = q N_left − p X, and a right label is (b, g) with b = q N_right − p m. One site shifts a (or b) by qs − p, so after a cell a returns to its residue and the label sets are the same in every iteration. Charge conservation at density p/q reads a + b = 0 across a bond and a + q(s₁ + s₂) + b = 2p across two sites.

[check] On a finite window of 2q sites with N = 2p and boundary labels (0, 0): 0 + q·2p − p·2q = 0 at the right end. The blockade and N are structural zeros, as in C176. Left-physical labels sit on left-canonical tensors and right-physical labels on right-canonical ones (C176 PROOF 1.1, 1.2).

### 6.2 Why not single-site VUMPS [derived, PROGRESS]

With mixed labels the dimensions d_L(a, f) = rank(P_f ρ_A P_f) and d_R(b, g) = rank(P_g ρ_B P_g) are properties of the state. VUMPS fixes both in advance, and the fixed-point condition A_C = A_L C = C A_R has no solution unless they match. ivumps.py stalls at error 4e-3 for that reason. Two-site updates set both through SVDs.

### 6.3 iDMRG [derived]

McCulloch's prediction with mixed labels is ψ_{t+1} = C_t B_t C_{t−1}⁺ A_t C_t. B_t maps the new mid right basis to the old mid right basis, and C_{t−1}⁺ maps that back to the old mid left basis, which is A_t's left basis. So every contraction is between matching bases. The guess is projected onto the blockade-allowed blocks, then brought to right-canonical form by a no-optimization SVD pass. One iteration does 2q − 1 two-site updates left to right and q − 1 back to the middle, then splits the centre at the middle bond into the mixed bond matrix C and B.

The two-cell unit (A_t, C_t B_t C_{t−1}⁺) closes on the basis of LP_{t−1}. Its transfer matrix T = T[A_t] M[C_t] T[B_t] M[C_{t−1}⁺] gives the correlation lengths ξ = −2q/ln|λ/λ₀|, where λ₀ ≈ 1 checks the closure (0.99999 to 1.00000 in the runs).

### 6.4 String as a gauge [proved]

The fermion string inserts (−1)^s per site, and s = (a_out − a_in + p)/q. So T^{PAR} = e^{iπp/q} G T G⁻¹ with G = diag(e^{iπa/q}). The |eigenvalues| are equal and the momenta shift by πn. Hence ξ₁ᵇ = ξ₁ᶠ, which I observe.

### 6.5 Finite-entanglement scaling [literature, derived interpretation]

At finite χ, S = (c/6) ln ξ + const (Calabrese–Cardy with ξ as the IR length; Pollmann et al. 2009). The level spacings of the transfer matrix scale like the scaling dimensions, ε_X = −ln|λ_X| ∝ Δ_X, at leading order.

In the boson frame the lowest dimensions are:
- ΔQ = 1: Δ₁ = (1/K_c + 1/K_s)/8 (the boson);
- ΔQ = 2: Δ₂ = 1/(2K_c) (the inter-valley pair);
- ΔQ = 0: Δ₀ = min(1/(2K_s), (K_c + K_s)/2, 1).
So ξ₂/ξ₁ = (1 + K_c/K_s)/4 and ξ₁/ξ₀ = Δ₀/Δ₁.

[check] For ff2v the fermion-frame dimensions (1/2, 1, 1) predict ratios 2 and 2. The runs read 2.13 and 2.17.

A gapped pseudospin would make ξ₁ (no ΔQ = 1 operator survives when φ_s is locked) saturate while ξ₂ (e^{i√2θ_c}) keeps growing. That is the registered gap signature.

### 6.6 Validation [numerical, validate_inf.out, local/out]

- The first window, swept by idmrg.py's own code, reproduces fibdmrg on 8, 16 and 32 sites to 1e-11 in energy and 1e-10 in S, and on a complex φ = 1.50 case to 5e-11.
- The infinite energy density converges to the exact value for c1 (error 7.6e-7 at χ = 96) and ff2v (3.5e-5 at χ = 96).
- c_fit = 1.015 (c1) and 1.989 (ff2v).
- The candidate's infinite energy density at n = 1/8 is −0.280922 at χ = 96, still decreasing with χ. C176's open chains give (E(256) − E(128))/128 = −0.280986 at the same N/L. The difference cancels the O(1) edge energy and leaves an O(1/L²) bias. The two agree to 6e-5, with the iDMRG value above, as a finite-χ variational energy should be.

## 7. Limits and controls in one place

- Frame-independent estimators: K_c from S(q) and K_s^P from C_P. Both are calibrated on ff2v (exact K = 1).
- Boson-frame estimators are calibrated on c2 and c1.
- Fermion-frame estimators are calibrated on ff2v but describe the wrong fixed point for this model.
- The infinite-chain code is calibrated on c1 and ff2v. The two-leg control c2 cannot be run in it (hidden U(1) × U(1), RESULT A1).
- What is not settled: the 30 percent parity-gap mismatch, Z_N, and correlation lengths beyond about 100 sites (Typhon pending).
