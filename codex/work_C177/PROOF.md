# C177 derivations and checks: the full qubit chain at the Fibonacci cusp

Claude Opus 5.5 agent, 2026-11-09. Conventions of C173 (PROOF 3): h = J_x XX + J_y YY + D(XY − YX) + J_z ZZ + (h_z/2)(Z⊗1 + 1⊗Z), Δ = J_x − J_y = 1, V = 4J_z, t_2 = −Δ²/(2V), a particle is a down spin, r = V/Δ. Scripts and outputs are named in each item.

## 1. Mapping the full chain onto the effective model at finite V/Δ

### 1.1 The exact one-particle band

[numerical, band.py, band.out, band.json] On a ring of L sites the dressed vacuum is the lowest even-parity state at k = 0, and the dressed one-particle state at momentum k is the lowest odd-parity state in that momentum sector, as long as the chemical potential is far below the band (at μ_eff = −2.2 three-particle states lie lower, which is the dilute gas itself). So ε(k) = E_odd(k) − E_vac is computed at μ_eff = −12 and the chemical potential is mapped separately (1.2). The vacuum is gapped with a correlation length of order one site, and L = 10, 12, 14 agree to 1e-6 relative. Fitting ε(k) to all harmonics up to order 4 gives, at the leading-order D = |t_2|/2 and J_x + J_y = 0:

| r = V/Δ | t_2 / t_2^LO | |w|/|t_2| | a_3/|t_2| (sin 3k) | b_4/|t_2| (cos 4k) | antisymmetric harmonics |
| :-- | :-- | :-- | :-- | :-- | :-- |
| 8 | 0.926 | 1.077 | 0.0055 | −0.021 | < 1e-5 |
| 12 | 0.964 | 1.035 | 0.0030 | −0.011 | < 1e-5 |
| 16 | 0.979 | 1.020 | 0.0018 | −0.0064 | < 1e-5 |
| 24 | 0.991 | 1.009 | 0.0008 | −0.0029 | < 1e-5 |
| 40 | 0.997 | 1.003 | 0.0003 | −0.0011 | < 1e-5 |

Consistency checks.
- [numerical] The renormalization of t_2 scales as the next even order, 1 − t_2/t_2^LO = 0.074, 0.036, 0.021, 0.0094, 0.0034, i.e. (4.7, 5.1, 5.3, 5.4, 5.5) × (Δ/V)², approaching a constant as it must if (3.2) is the leading term of a series in (Δ/V)². This matches C173 PROOF 3.2, whose many-body check fell as V⁻³ in absolute energy, i.e. (Δ/V)² relative to t_2.
- [numerical] The NN hop is almost unrenormalized: |w|/|w^LO| = 1.077 × 0.926 = 0.997 at r = 8. The deviation of |w|/|t_2| from 1 is the t_2 renormalization.
- [numerical, proved in 1.3] The antisymmetric harmonics (cos k, sin 2k, cos 3k, sin 4k) vanish to fitting precision at every r. The band keeps the exact symmetry k → π − k of C173 PROOF 3.3 to all orders, as the lattice symmetry of 1.3 requires.
- [numerical] The extra harmonics allowed by the symmetry (sin 3k, cos 4k) are at most 2 percent of |t_2| at r = 8.

Caveat (added after 3.1). All of this is at μ_eff = −12. The energy denominators of the virtual processes contain the chemical potential, so the renormalizations are μ-dependent at the same order (Δ/V)². At the physical μ the V/Δ = 8 band has |t_2|/t_2^LO = 0.98 to 1.00 and, with the D of 1.2, |w|/|t_2| = 0.93 to 0.95 (twob_phys_*.out).

So at J_x + J_y = 0 the full chain at any r tested has the two-valley band of (3.2) with degenerate, equal-mass valleys, |t_2| renormalized by a few percent, and small symmetric extra harmonics. The interactions V_2, V_3 and the blockade are renormalized at the same relative order (Δ/V)², which is not measured here.

### 1.2 Tuning D and h_z

[numerical, tune.py, tune.json] D is rescaled until the measured |w|/|t_2| = 1 (three fixed-point iterations, residual below 1e-12). The constant c_0 of the fitted band at μ_eff = −12, −9, −6 is linear in h_z to 5e-4 (r = 8) and 6e-7 (r = 40) in absolute units, and h_z is set where −c_0/|t_2| = −2.2, i.e. where the dressed one-particle chemical potential equals C173's μ_eff. The result:

| r | D | h_z | |t_2| | D/D^LO |
| :-- | :-- | :-- | :-- | :-- |
| 8 | 0.0290182 | −3.8791275 | 0.0578845 | 0.9286 |
| 12 | 0.0201223 | −5.9199091 | 0.0401863 | 0.9659 |
| 16 | 0.0153165 | −7.9401807 | 0.0306059 | 0.9803 |
| 24 | 0.0103233 | −11.9602723 | 0.0206380 | 0.9910 |
| 40 | 0.0062296 | −19.9762173 | 0.0124573 | 0.9967 |

At r = 40 this reproduces C173's candidate bond to 3e-5 in h_z and 0.3 percent in D, as it should.

Density. The raw density ⟨(1 − Z)/2⟩ contains a background of virtual pairs. [numerical, ring ED at L = 12, μ_eff = −12] The vacuum density is 0.023, 0.012, 0.0072, 0.0033, 0.0012 per site at r = 8, 12, 16, 24, 40, close to the second-order value 2(Δ/V)² (0.031, 0.014, 0.0078, 0.0035, 0.0013). At finite density part of this background is blocked by the particles. So at r = 8 the dressed particle density is the raw density minus between about 0.01 and 0.023 [heuristic estimate of the blocking].

### 1.3 Symmetry, gauge and the valley order parameter

[proved, algebra] Let S be site inversion followed by a π/2 rotation about z on every site and a π rotation about z on odd sites. On a bond, inversion swaps the two factors, the π/2 rotation maps X → Y, Y → −X, and the π rotation on one site flips the sign of its X and Y. The signs picked up by each term are

| term | inversion | π/2 on all sites | π on one site | total |
| :-- | :-- | :-- | :-- | :-- |
| XY − YX (D) | − | + | − | + |
| XX − YY (Δ) | + | − | − | + |
| XX + YY (J_x + J_y) | + | + | − | − |
| ZZ, Z | + | + | + | + |

So S is a symmetry of the full bond exactly when J_x + J_y = 0, to all orders in Δ/V (this is C173 PROOF 3.3). The bond operator P_j = X_j X_{j+1} + Y_j Y_{j+1} is odd under S. In the effective model S acts as k → π − k and exchanges the valleys, and P_j is (twice) the real part of the NN hopping expectation, which is the valley imbalance. [derived] So P is the order parameter of a single-valley (valley-polarized) state, and ⟨P_j⟩ = 0 in any S-symmetric state.
- [numerical, fc.py checked against ED at L = 10] On the open chain at J_x + J_y = 0 the ground state has ⟨P_j⟩ = 0 to 1e-16. With J_x + J_y = 0.01 (about 0.17|t_2|) and a transverse field, ⟨P_j⟩ is of order 0.2: the valley imbalance responds strongly to the symmetry-breaking term, as expected for two degenerate valleys.
- The XZ + ZX coupling and a transverse field are not S-symmetric, so the generic point G breaks S, parity and J_x + J_y = 0 at once.

[proved, algebra; checked numerically] The site rotation U = Π_j exp(−iθ_j Z_j/2) with θ_j = (−1)^j π/4 maps σ⁺_j → e^{−iθ_j} σ⁺_j. The hop D(XY − YX) = 2iD(σ⁺_jσ⁻_{j+1} − σ⁻_jσ⁺_{j+1}) picks up e^{∓i(θ_j − θ_{j+1})} = e^{∓iπ/2}, which makes it real with an alternating sign, and the pairing Δ(σ⁺σ⁺ + σ⁻σ⁻) picks up e^{∓i(θ_j + θ_{j+1})} = 1. ZZ and Z are invariant. So on the J_x + J_y = 0 stratum the bond is real in this gauge (the imaginary part of the rotated 4 × 4 bond is 4e-17), the open-chain DMRG can run in real arithmetic with a period-2 MPO, and all observables are mapped back with the phases e^{−iθ_i + iθ_j}. Entanglement is gauge invariant.

## 2. Numerical method

- [numerical, ed_check.out, ed_check2.out] fc.py (open chain of the bond with the full on-site field on every site) reproduces sparse ED to 1e-11 at L = 10, 12, 14 for F_8, F_16 and the generic G_8, first with C166's dmrg.py, then with dmrg2.py.
- dmrg2.py differs from C166's dmrg.py only in the local solver and the per-bond MPO. ARPACK needed about 530 matvecs per two-site solve at r = 8, because the two-site problem has a spectral spread of order V and gaps of order t_2/L, so Lanczos converges at a rate set by (gap/spread)^{1/2}. Davidson with the diagonal of H_eff as preconditioner takes the V n_j n_{j+1} term exactly, which removes most of that stiffness, and needs 18 to 50 matvecs per solve in converged sweeps. At L = 24 both codes give the same energy (to 1e-14), S̄ and density.
- vumps_fc.py runs C166's vumps.py (unchanged) with the same MPO and a constant shift. In damped mode it converges (err 5e-6 at χ = 32), in plain mode it stalls at err 1e-3, so production uses damped mode.

## 3. The two-magnon bound state at K = π in the full chain (coordinator's check of C175)

C175 found that in (3.2) at φ = π/2 one particle from each valley binds at total momentum π, with E_b = 0.0195|t_2| at |w| = 1, because the unlike pair sits near an exact zero-energy resonance (the resonance function D = 2V_2V_3 + 4V_2 + V_3u² + 4V_3 + 8, units |t_2| = 1, u = |w|, vanishes at V_2 = 2t_2, V_3 = 0, and the generated V_3 = (2/3)t_2 moves it to D = −2/3, the bound side). The binding is small and the resonance is a fine-tuned feature, so corrections of relative order (Δ/V)² to V_2, V_3 can move D through zero. This section measures them.

Method [numerical, twob.py, twob_fit.py; outputs twob_full.out, twob_full.json, twob_fit.out, twob_fit.json].
- Full chain: momentum-resolved ring ED at L = 14, 16, 18, 20, even parity, K = π, at the band-fit couplings (μ_eff = −12, so four-particle states are far above two-particle ones). E_2(K = π) − E_vac for the two lowest levels.
- Effective two-body model on the same rings: two non-passing hard-core bosons with the blockade, the exact one-particle band of 1.1 (hops t_1..t_4 from the fitted harmonics), and short-range interactions V_2, V_3 (and V_4). Bosonic states are selected by f(d) = e^{iKd} f(L − d). Check: with the bare band and V_2 = −2, V_3 = −2/3 it gives E − 2ε_min = −0.01953 at L = 240 and 400, C175's value.
- With the bare V_2 = 2t_2, V_3 = (2/3)t_2 (and the exact band), the effective model lies below the full chain by a nearly L-independent amount, 0.075 to 0.08|t_2| at r = 8, 0.03 at r = 16, 0.005 at r = 40, for the lowest level. The full chain is less attractive.
- Fitting V_2, V_3 (model A) or V_2, V_3, V_4 (model B) to the 8 levels (two per L) reproduces them to rms 1e-3 and 2e-4|t_2| at r = 8, and 3e-5 and 2e-5 at r = 40. Then the two-body problem is solved at L = 400 with the fitted interactions.

| r = V/Δ | V_2/|t_2| (A) | V_3/|t_2| (A) | V_4/|t_2| (B) | D (A) | (E_2 − 2ε_min)/|t_2| at L = 400 (A; B) |
| :-- | :-- | :-- | :-- | :-- | :-- |
| 8 | −1.733 | −0.493 | −0.015 | +0.31 | +0.0005; +0.0005 (no bound state, lowest scattering level) |
| 12 | −1.867 | −0.576 | −0.008 | −0.20 | −0.0011; −0.0010 |
| 16 | −1.922 | −0.613 | −0.005 | −0.40 | −0.0060; −0.0059 |
| 24 | −1.964 | −0.642 | −0.003 | −0.54 | −0.0124; −0.0124 |
| 40 | −1.987 | −0.657 | −0.001 | −0.62 | −0.0168; −0.0167 |
| ∞ (C175) | −2 | −0.667 | 0 | −0.667 | −0.0195 |

Consistency checks.
- [numerical] The interaction corrections scale as (Δ/V)²: D + 2/3 = 0.98, 0.47, 0.27, 0.12, 0.045 at r = 8 to 40, i.e. (62, 67, 70, 71, 72) × (Δ/V)², approaching a constant. This is the order of the next Schrieffer–Wolff terms, as in 1.1.
- [numerical] At r = 40 the fitted V_2, V_3 differ from the bare values by 0.7 and 1.4 percent, and the binding (0.0168) differs from C175's 0.0195 by the band renormalization and these small shifts. Models A and B agree to 1e-4 in the binding at every r, and the fitted V_4 is small, so the result does not hinge on the interaction range.
- [numerical] Fixing V_3 at its bare value and fitting only V_2 (model C) gives the same signs (unbound at r = 8, bound at r ≥ 12) with a worse fit. Fixing V_2 and fitting only V_3 (model D) fits badly (rms 0.066 at r = 8) and is discarded.

Result [numerical, matched effective two-body model; the extrapolation from L ≤ 20 to L = 400 assumes the corrections are short-ranged, which the small V_4 supports]. The bound pair of C175 is present in the full qubit chain at V/Δ ≥ 12, with a binding that grows with V/Δ toward C175's value (0.001, 0.006, 0.012, 0.017|t_2| at r = 12, 16, 24, 40). It is absent at V/Δ = 8. With D(r) ≈ −2/3 + 70 (Δ/V)², D changes sign near V/Δ ≈ 10. The (Δ/V)² corrections are repulsive and push the unlike pair through the resonance, so at small V/Δ the pairing mechanism is switched off at the two-body level.

What this does and does not imply.
- [derived from C175's criterion] At r = 8 the dilute gas has no bound pair. Whether it is two-component then depends on C175's scattering criterion Δ_σ ∈ (−1, 1] evaluated with the r = 8 interactions, which is not computed here.
- At r ≥ 12 the pair is bound, and C175's dilute conclusion (pair liquid, c = 1) applies in the dilute limit there. The pair binding at r = 12 to 16 (0.001 to 0.006|t_2|) is smaller than at r = 40, so the pair gap at n ≈ 0.17 should be even smaller than C175's estimate [heuristic], and a c ≈ 2 reading at accessible L would not distinguish pairing from a two-component liquid there.
- (Withdrawn, see 3.1.) The window of V/Δ around 8 to 10 without a bound pair was an artefact of matching at μ_eff = −12.

### 3.1 Revision after the C186 conflict (physical chemical potential)

[numerical, twob_phys.py; outputs twob_phys_mu-2.4.out, twob_phys_mu-2.4_retune.out, twob_phys_mu-2.6.out] The matching above was done at μ_eff = −12. At the physical μ (μ_eff = −2.4, −2.6, just below the band bottom) the excitations carry ⟨N⟩ = 1.3 to 1.8 (one excitation, V/Δ = 8), so near-resonant ΔN = 2 processes dress them strongly, and the effective N-conserving couplings depend on μ. Measured directly, B(L) = E_2(π) − 2 min E_1 turns negative at L = 18 to 20 for V/Δ = 8 and 12, less attractive than (3.2) at the same L but with the same trend, and an N-conserving fit at the physical μ (poor, rms of order 1e-2|t_2|) extrapolates to weak binding. The table above therefore describes the chain at μ_eff = −12 only, and the conclusion "unbound at V/Δ = 8, D crosses zero near V/Δ ≈ 10" is withdrawn. The (Δ/V)² scaling of the μ_eff = −12 shifts is still correct as a statement about that chemical potential. See RESULT R2.

## 4. Why the parity-breaking terms are dangerous (heuristic)

[heuristic] A transverse field h_x and an XZ + ZX coupling create or remove one particle at momentum 0. In the dilute cusp gas at μ_eff = −2.2 the one-particle cost at k = 0 is ε(0) − μ = −2 + 2.2 = 0.2|t_2| (band ε = −2 sin k − 2 cos 2k in the physical frame, minimum −2.25 at sin k = 1/4), only 0.25|t_2| above the band bottom, and the valley momentum k_a = arcsin(1/4) = 0.253 is small. So the vertex is close to resonance in both energy and momentum. In bosonization the single-particle operator e^{iθ_a} has a small scaling dimension (of order 1/(4K)), so even a momentum-mismatched source competes with the mismatch energy v k_a only through h^{1/(2 − 1/(4K))}, which falls slowly with h. A commensurate-incommensurate (floating) window must exist for small enough h if the mismatch is finite, but it may be very narrow. The numerics of RESULT R6 (order-one change of the L = 48 state at h = 0.02|t_2|, a gapped state with density locked at 1/5 at 0.25|t_2|) fit a narrow or absent window. None of this is derived, and the width of the window is the open question for openness in B₂.

## 5. Results

See RESULT.md (R1 to R7).
