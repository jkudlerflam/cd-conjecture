# C183 derivations and checks: the full qubit chain near the Fibonacci cusp

Claude Opus 5.5 agent, 2026-11-09. Conventions as in RESULT.md (C173 bond convention plus J_xz, K_xz, h_x), Δ = J_x − J_y = 1, V = 4J_z, r = V/Δ, |t₂| the measured second harmonic of the dressed band, a particle is a down spin.

## 1. Tuning the full chain anywhere near the cusp

[numerical, tune183.py] The dressed one-particle band ε(k) = E_odd(k) − E_vac of the full ring (L = 12, parity split, field far below the band) is fitted with all harmonics to order 4. Its first harmonic is a₁ sin k + a_c1 cos k, and the effective model's 2|w||t₂| cos(k − φ) identifies a₁ = 2|w||t₂| sin φ (set by D, a₁ = 4D at leading order) and a_c1 = 2|w||t₂| cos φ (set by J_x + J_y, a_c1 = 2(J_x + J_y) at leading order). A fixed-point iteration on (J_x + J_y, D) reaches the target |w|, φ to 1e-12. φ and π − φ are equivalent (k → π − k maps 2|w| cos(k − φ) to 2|w| cos(k − (π − φ)) and fixes cos 2k), so φ ≤ π/2 is used.

Consistency checks.
- [numerical] At r = 8, |w| = 1, φ = π/2 it returns D = 0.0290182 and n_vac = 0.0229, C177's values (C177 PROOF 1.2).
- [numerical] The valley offset grows linearly off the stratum, 0.279|t₂| at φ = 1.50 and 0.476|t₂| at φ = 1.45 (r = 8, |w| = 1), i.e. δ ≈ 4|t₂|(π/2 − φ), the derivative of 2|w| cos(k − φ) at the valleys k = π/2 ∓ k₀ with sin k₀ = √15/4 (δ = 4|w| sin k₀ |π/2 − φ| ≈ 3.9|t₂||π/2 − φ|).
- [numerical] At |w| = 3 and φ = 1.45 the upper valley has disappeared (one minimum), as it must once the offset exceeds the barrier.

The dressed chemical potential x is the distance of the field from the band bottom, x = −min_k ε(k)/|t₂|. At |w| = 1, φ = π/2 C173's μ_eff = −2.2 is x = 0.05.

## 2. Stage T, the two-body problem of the full chain on the stratum

### 2.1 Method

[numerical, twob183.py, an extension of C177 twob.py] For each (r, |w|) the full chain is tuned (§1), and the two lowest bosonic two-magnon levels at K = π are computed by momentum-resolved ring ED at L = 14, 16, 18, 20. An effective two-body model with the exact fitted band and short-range interactions V₂, V₃ (model A) or V₂, V₃, V₄ (model B) is fitted to these eight levels and solved at L = 400. C175's threshold scattering code (scat.py, bare band 2|w| cos(k − π/2) − 2 cos 2k) is then run with the model-A V₂, V₃ to give the like-valley length a_s, the Hermitian unlike length matrix with diagonal A_d and off-diagonal |A_o|, and Δ_σ = (a_s − A_d)/|A_o|, plus the lowest two-body level over 49 values of K.

Validation. At r = 8, |w| = 1 the fit reproduces C177 (V₂ = −1.7334, V₃ = −0.4930, E₂ − 2ε_min = +0.00046|t₂|, local_twob_check.json). With the bare V₂ = −2, V₃ = −2/3 the scattering code gives Δ_σ = −1.339, and with V₂ = V₃ = 0 it gives 0.250, C175's calibration values.

### 2.2 Structure through the resonance

[derived] Let the unlike 2 × 2 length matrix have eigenvalues a_reg (regular) and a_div (the channel that passes through the zero-energy resonance). Then A_d = (a_reg + a_div)/2 and |A_o| = |a_div − a_reg|/2. As a_div → ±∞,

    Δ_σ = (a_s − A_d)/|A_o| = [2a_s − a_reg − a_div]/|a_div − a_reg| → −sgn(a_div) + 2(a_s − a_reg)/|a_div| + O(a_div⁻²).

On the bound side the resonant length is large and positive (it is the size of the shallow bound state, a₊ ≈ 21 at |w| = 1, C175), so Δ_σ → −1 from below when a_s < a_reg. On the unbound side a_div is large and negative, so Δ_σ → +1 from below. So crossing the resonance from the bound to the unbound side moves the dilute pseudospin from just inside the ferromagnetic (vector-chiral) side, Δ_σ ≲ −1, to just inside the easy-plane window near its isotropic antiferromagnetic end, Δ_σ ≲ +1. Away from the resonance on the unbound side Δ_σ relaxes toward the hard-rod value (0.25 at |w| = 1, V₂ = V₃ = 0). The bound pair and the easy-plane two-component gas are thus the two sides of one resonance, and the (Δ/V)² corrections of the full chain decide the side.

Check [numerical, twob_map.out]. The stage T data follow this. Every point on the bound side (A_d > 0) has Δ_σ below −1 (from −1.003 at r = 10, |w| = 1, where A_d = 727 and the pair sits on the resonance, to −2.43 at r = 16, |w| = 3), every point on the unbound side (A_d < 0) has Δ_σ between 0.36 and 0.96, and on each side the points closest to the resonance (largest |A_d|) have Δ_σ closest to ∓1. Example, r = 6, |w| = 1.5 (unbound): a_div = A_d − |A_o| = −46.4, a_reg = 0.75, a_s = −0.84 < a_reg, Δ_σ = 0.933. r = 6, |w| = 2 (bound): a_div = 41.7, a_reg = 0.86, a_s = −0.68, Δ_σ = −1.075.

### 2.4 Where the full chain unbinds

[numerical, twob_map.out; the scaling form heuristic] On the stratum the unlike pair is unbound for r < r_b(|w|) with r_b between 16 and 20 at |w| = 0.5 (unbound at 16, A_d = −74, close to resonance), r_b ≈ 10 at |w| = 1 (on the resonance at r = 10, C177 found a sign change near 10), r_b ≈ 7 at |w| = 1.5, r_b ≈ 5.5 at |w| = 2, and r_b < 4 at |w| = 3. So r_b|w| ≈ 10 within about 10 percent. A heuristic reason: C175's resonance function at the Schrieffer–Wolff couplings is D = −2u²/3 (u = |w|), and C177 measured the (Δ/V)² shift of D as about +70/r² at u = 1. If that shift depends weakly on u, D changes sign at r_b = (3c/2)^{1/2}/u ≈ 10/u, which fits the data with c between 60 and 80 across u. In bond language, |w| = 4DV/Δ² (the leading-order map), so the full chain is on the unbound side of the resonance when D V²/Δ³ ≲ 2.5.

This confirms prediction T1 (r_b decreasing in |w|, r_b(1) between 8 and 12, |w| = 3 bound for r ≥ 6). The unbound region is the low-|w|, low-r corner, which is where the first-order hop is smallest compared with the generated interactions' (Δ/V)² corrections.

### 2.3 Scope of the dilute criterion

[derived from C175 PROOF 1 and C180 PROOF 4, with assumed inputs] Δ_σ ∈ (−1, 1] and no bound state at any K give a two-component impenetrable gas (c = 2) at leading order in n, valid for n|A| ≪ 1, with |A| the largest length. Near the resonance |A_d| and |A_o| are tens of sites, so the window shrinks to n ≲ 0.01 to 0.05. Three inputs are not checked here: (i) the three-body problem (trimers bind at the cusp at r = ∞, C175), (ii) the extrapolation of the fitted V's from L ≤ 20 rings to the threshold, which assumes the corrections are short ranged (V₄ is small, below 0.08|t₂|), (iii) at r = 4 and 5 the effective description itself, since Δ/V = 0.25 and 0.2 and the vacuum carries 5 and 3.5 percent virtual particles.

### 2.5 The field dependence of the fitted interactions

[numerical, twob_mu.out, twob_mu_small.out; mechanism derived by C186] The generated interactions come from virtual pairs whose energy denominators contain the field, of the form kV − 2μ in C186's notation, so the (Δ/V)² corrections that decide the side of the resonance depend on where the field sits. Refitting at |w| = 1 with the field raised from μ_eff = −12 to −4 makes V₂ more attractive by 0.03 to 0.05|t₂| per unit of μ_eff and moves the resonance from r ≈ 10 (μ_eff = −12) to ≈ 8 (−8), ≈ 6 (−5) and below 4 (−4). The two-body extraction degrades at μ_eff ≥ −4 (fit rms 2 to 6e-2), where states with four particles come within reach of the two-magnon levels on rings of 14 to 20 sites, so the onset itself (μ_eff ≈ −2.2) is reached only by the trend. Every quantity of §2.1 to §2.4 is therefore a property of the chain far below the band, and does not describe the gas.

## 3. The many-body onset in the full chain, and what the parity gap measures

### 3.1 The dressed vacuum below the onset

[numerical, local/Pr5_32.json, local/Pr5_64.json, local/p5, local/p8, local/p8b] The full chain does not conserve N, only its parity, and its vacuum carries virtual pairs whose density grows as the field approaches the band (0.035 per site far below at r = 5, about 0.06 to 0.10 near the onset). Below the many-body onset the chain is this gapped vacuum. At r = 5, |w| = 1, x = 0.45 the open chain has S̄ = 0.464 (L = 32) and 0.445 (L = 64), a flat bulk S of 0.44, a uniform bulk density of 0.100 with no Friedel oscillation, |⟨s⁺_i s⁻_j⟩| decaying exponentially with a correlation length of about 10 sites (nonzero at odd separations only), and truncation 1e-12 at χ = 160. Its odd-parity partner lies 8.8e-4 = 0.0095|t₂| higher at L = 32, which is the one-quasiparticle gap there. So "n_raw − n_vac" with n_vac taken far below the band is not a gas density near the onset, and the targeting of the Typhon screen was amended to count from the onset (RESULT.md).

The onset is sharp in the L = 24 ladders: at r = 8, |w| = 1, φ = π/2, S̄ goes from 0.26 (x = 0.05) to 0.93 (x = 0.10); at φ = 1.555 from 0.36 (x = 0.10) to 1.05 (x = 0.15). Above it, the gas density grows to about 0.3 at x = 0.3 and then jumps to the dense blockade-packed phase (n ≈ 0.45, S̄ ≈ 0.65) by x = 0.45 (r = 8) or between x = 0.6 and 1.0 (r = 5). So the gas window in the full chain is roughly 0.05 ≲ n ≲ 0.3 above the onset, bounded above by a first-order jump to the packed phase. [heuristic] The jump is the many-body face of the generated attractions V₂, V₃, which favour dense configurations even where the pair is unbound.

### 3.2 The parity gap of the full chain

[derived] With N not conserved, the open-chain ground states in the two parity sectors at fixed field give Δ_p = E_odd − E_even. In a phase with a conserved N at low energy (the gas) and chemical potential μ set by the field, E_odd − E_even = min over the odd neighbours of [E(N ± 1) − E(N) ∓ μ], which is half of E(N+1) + E(N−1) − 2E(N) when μ sits midway between the two thresholds and smaller otherwise. In a gapless one-component or two-component liquid this is O(v/L) and changes sign as μ moves through the finite-size levels. In a pair liquid it stays at Δ_s/2 or above for every μ in the window. So the registered test (L·|Δ_p| growing) detects a pair gap only when Δ_s exceeds the level spacing at the sizes used, and a field scan would be needed to turn Δ_p into the symmetric gap. This is a limitation of the test, not of the criterion.

## 4. Three-body check at unbound points

[numerical, three183.py, three183.out; C175's nbody_box.py with the bare band and the model-A V₂, V₃ of stage T] For three particles on the line in relative coordinates (gaps ≤ R), at the aab and abb total momenta and at K = 0, π, the lowest level E₀ − 3ε_min is positive and falls as 1/R² between R = 150 and 250 (ratio 2.76 against (250/150)² = 2.78), i.e. it is the box scattering level, at all seven unbound points tested: (r, |w|) = (4, 1), (5, 1), (5, 0.5), (6, 1), (6, 1.5), (8, 1), (8, 0.5). The ratio is 2.5 rather than 2.78 at (8, 1) and (6, 1.5), the two points closest to the pair resonance, where the residual attraction slows the approach to threshold, but the level stays positive (+0.0018 and +0.0016 at R = 250). At K = 0 and π (not valley momenta) the lowest level stays far above threshold. Validation: with the bare V₂ = −2, V₃ = −2/3 at |w| = 1 the same code gives −0.020062 at R = 200 (C175: −0.020062). So the trimer of C175 (bound at r = ∞) is absent where the pair is unbound, at these points.
