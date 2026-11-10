This run was done by a Claude Opus 5.5 agent (C190, 2026-11-09), not by Codex. Astra was not used.

# C190 derivations: the qubit-to-effective-model map at the candidate bonds, and the solver

Tags [proved], [derived], [numerical], [heuristic] open the paragraph they apply to.

## 1. The bonds

All four bonds are in the parity class, with no transverse DM term. On an open chain of L sites,

    H = Σ_{j=0}^{L-2} [ J_x X_j X_{j+1} + J_y Y_j Y_{j+1} + D (X_j Y_{j+1} − Y_j X_{j+1}) + J_z Z_j Z_{j+1} ] + h_z Σ_{j=0}^{L-1} Z_j,   (1.1)

with the field on every site, endpoints included (the collaborator's convention, contrib/d2_dmrg_priority_20261109). The bonds are listed below.

| name | J_x | J_y | D | J_z | reference h_z |
| :-- | :-- | :-- | :-- | :-- | :-- |
| p153 (collaborator, unrounded) | (1 + cos 1.53/80)/2 | (−1 + cos 1.53/80)/2 | sin(1.53)/160 | 10 | −19.97589125 |
| rat (collaborator's rational simplification) | 2001/4000 | −1999/4000 | 1/160 | 10 | −199759/10000 |
| sym40 (symmetric line, C173's candidate bond) | 1/2 | −1/2 | 1/160 | 10 | tuned |
| sym80 (symmetric line, V/Δ = 80) | 1/2 | −1/2 | 1/320 | 20 | tuned |

rat is the bond 10 ZZ + (XX − YY)/2 + (XX + YY)/4000 + (XY − YX)/160 with uniform field −199759 Z/10000, as in the spec.

## 2. The map to C173's effective model

[derived; C173 PROOF 0.1 and 3.1, restated with the numbers of these bonds] With n_j = (1 − Z_j)/2 (a down spin is a particle) and the Jordan–Wigner fermions c_j = (Π_{k<j} Z_k) σ⁺_j, the bond (1.1) is C152's fermion form

    H = Σ_j [ w c†_{j+1} c_j + h.c. − Δ (c_j c_{j+1} + h.c.) + V n_j n_{j+1} ] − μ Σ_j n_j + const,
    w = J_x + J_y + 2iD,   Δ = J_x − J_y,   V = 4 J_z,   μ = 4 J_z + 2 h_z.                                  (2.1)

Check of the ZZ and field terms: Z_j Z_{j+1} = 1 − 2n_j − 2n_{j+1} + 4 n_j n_{j+1}, so J_z Σ ZZ gives V = 4J_z and −2J_z per bond on each end of the bond (−4J_z per bulk site), and h_z Z_j = h_z − 2 h_z n_j. The bulk single-particle term is then −(4J_z + 2h_z) n_j = −μ n_j. At the two edge sites only one bond contributes, so the edge chemical potential is μ − 2J_z, lower by V/2: an open end repels particles. The effective model (C176) has no such term, which is one source of the constant offset in S̄ seen by C177.

Second order in Δ/V, inside the blockaded space (no two adjacent particles), C173 PROOF 3.1 gives

    t₂ = −Δ²/(2V)  (next-nearest-neighbour hop, no Jordan–Wigner sign),
    V₂ = −Δ²/V = 2t₂,  V₃ = −Δ²/(3V) = (2/3) t₂,
    a per-particle shift +3Δ²/V (two blocked bonds and two raised bonds), absorbed into μ.                    (2.2)

In units of |t₂| this is C173 eq. (3.2), with

    w/|t₂| = (J_x + J_y + 2iD) · 2V/Δ²,   μ_eff = (μ − 3Δ²/V)/|t₂| = (4J_z + 2h_z − 3Δ²/V) · 2V/Δ².          (2.3)

Numbers [numerical, models.eff_params]:

| bond | Δ | V | |t₂| | |w|/|t₂| | φ = arg w | μ_eff at the reference h_z |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| p153 | 1 | 40 | 0.0125 | 1 | 1.53 | −2.1426 (h_z = −19.97589125) |
| rat | 1 | 40 | 0.0125 | 1.0008 | 1.5308 (arctan 25) | −2.1440 (h_z = −19.9759) |
| sym40 | 1 | 40 | 0.0125 | 1 | π/2 | −2.2 (h_z = −19.97625, C173) |
| sym80 | 1 | 80 | 0.00625 | 1 | π/2 | −2.2 (h_z = −39.988125) |

Inverting (2.3), h_z = (μ_eff |t₂| + 3Δ²/V − 4J_z)/2 (models.hz_for_mu_eff).

Corrections [literature in this project, C173 PROOF 3.2 and C177 R1]. The next even order changes energies by relative (Δ/V)², 6e-4 at V/Δ = 40 and 1.6e-4 at 80. C177 measured |t₂| renormalized by 1 − 5.5(Δ/V)² (0.997 at V/Δ = 40), and the ΔN = 2 terms inside the blockaded space are of relative order Δ/V. The virtual-pair background adds at most about 2(Δ/V)² to the raw density, see 3 below. So at V/Δ = 40 the map is faithful to a few parts in a thousand in the couplings, but the absolute chemical potential carries an O((Δ/V)²|t₂|) shift that is not negligible against the width of the density window, which is why h_z is tuned on the full chain rather than taken from (2.3).

## 3. Which density to match

[derived] The qubit magnetization maps to the particle density by n_j = (1 − ⟨Z_j⟩)/2 exactly. This raw density counts the virtual pairs of the dressed vacuum. To first order in Δ/V each free bond (both ends empty and both outer neighbours empty) carries a pair with amplitude Δ/(V − 2μ) ≈ Δ/V, so the vacuum has n_virt ≈ 2·(Δ/V)² per site in the dilute limit, 1.25e-3 at V/Δ = 40 and 3.1e-4 at 80. In the gas the fraction of free bonds is reduced (about 1 − 3n at low density), so the background in the gas at n ≈ 0.2 is below 1e-3. This is under 0.5 percent of the target density, and C176's S̄ changes with density as dS̄/dn ≈ 2.4 at L = 256 (C176 density_corrected table), so it moves S̄ by at most 0.002, below every tolerance used here. The registered matching is therefore on the raw middle-half density n_mid (the C176 definition, mean of n_j over j ∈ [L/4, 3L/4)), and the background estimate is reported alongside.

Targets. The spec fixes C176's L = 256 values: n_mid = 0.1957 (N/L = 3/16) and 0.1326 (N/L = 1/8). C176's effective-model chemical potentials there, from E(N ± 1) at the same χ, are μ_eff = −2.1654 (3/16, L = 256), −2.1584 (3/16, L = 128), −2.2147 (1/8, L = 256). The full chain is run grand canonically at fixed h_z, so its n_mid is nearly independent of L (C176's fixed-N series instead has n_mid falling with L, 0.220, 0.204, 0.196 at L = 64, 128, 256 for 3/16). The comparison with C176 is therefore made at equal n_mid at each L, by interpolating C176's S̄ across its density families at that L (C176 analyze.density_corrected), and not at equal N/L.

## 4. Parity gap in the grand canonical chain

[derived] The full chain conserves only P = Π_j Z_j = (−1)^N. At fixed h_z, let E_p be the lowest energy in sector p. If the chain were exactly N-conserving with convex E(N) and N₀ the even minimizer of E(N) − μN, then

    E_odd − E_even = min_± [E(N₀ ± 1) − E(N₀) ∓ μ] ≤ ½ [E(N₀ + 1) + E(N₀ − 1) − 2E(N₀)] = D_P/2,

with equality when μ sits midway. So for an unpaired liquid L·(E_odd − E_even) is bounded by L·D_P/2 (C176: L·D_P ≈ 2.2 to 2.4 |t₂|), and its sign alternates with the position of μ between levels. A pair gap Δ_s adds Δ_s to every odd-sector state, so L·(E_odd − E_even) grows like L·Δ_s once L·Δ_s exceeds the level spacing. This is the reading registered in RESULT.md. The near-conservation of N in the full chain (violated only through virtual pairs) does not change the bound at the level of accuracy needed here.

## 5. The solver

[numerical, z2dmrg.py, validated in validate.py] Two-site DMRG with the parity built into every tensor: bond labels are the parity of the particle number to the left, sorted on each bond, so each tensor is a set of blocks and the parity sector is exact without a penalty (C177 used λΠZ). Environments and the two heavy contractions of the two-site matvec are done block by block (about a quarter of the dense cost); the SVD splits into two blocks.

The local eigensolver is the point where earlier full-chain runs (C177 at V/Δ ≥ 16, and the collaborator) lost convergence: local residuals stalled at 1e-6 to 1e-5. The two-site problem has a spread set by V (adjacent pairs, three-in-a-row intermediates) and low-lying gaps of order 0.04|t₂| = 5e-4 (measured at L = 32, diag_local.py). With the diagonal preconditioner Davidson needed 1200 to 1500 matvecs per solve to reach 1e-10 from a converged start, whatever the floor on the denominator (diag_sigma.py), and a Kronecker-sum preconditioner without the central-bond ZZ was worse because it puts doubly occupied central bonds at low energy.

The preconditioner used here is block diagonal in the two physical indices (s₁, s₂). On the block (q_a, s₁, s₂, q_c) of the two-site wavefunction,

    H₀ = H_L^{(s₁)} ⊗ 1 + 1 ⊗ H_R^{(s₂)} + c(s₁, s₂),
    H_L^{(s₁)}[a', a] = Σ_{w even} Le[a', w, a] W_j[w, end, s₁, s₁],   H_R^{(s₂)}[c', c] = Σ_{u even} W_{j+1}[start, u, s₂, s₂] Re[c', u, c],
    c(s₁, s₂) = W_j[start, Z, s₁, s₁] W_{j+1}[Z, end, s₂, s₂]   (the central ZZ),                               (5.1)

so H₀ is the exact Hamiltonian of everything except the terms that flip s₁ or s₂ (pairing and hopping on the three bonds touching the two sites). Each H_L^{(s₁)}, H_R^{(s₂)} is diagonalized once per local problem (8 Hermitian eigenproblems of size about D/2), and the preconditioner (H₀ − E₀ + σ)^{-1}, with E₀ = min(θ, λ_min(H₀)) and σ = 0.01, is applied in that eigenbasis at a cost of about 40 percent of a matvec. On the same local problem as above it reaches 1e-9 in 95 matvecs from a perturbed converged start and in 150 from a random start (diag_sigma.py), against more than 1500 for the diagonal. [heuristic] This works because H₀ places the V-cost states correctly (the doubly occupied central bond through c, an occupied site next to an occupied block edge through the ZZ channel of H_L), so the preconditioned spectrum of the high states is O(1), and on the low states H₀ differs from H only by the second-order shifts of a few |t₂|.

White's density-matrix mixer (perturbation from the Z, Sp, Sm channels, normalized, amplitude 1e-3 to 1e-5) is used in the early sweeps. It is needed because N is not conserved but changes only through virtual pairs, so a local update can move particles across a cut only if the environment basis already holds the neighbouring particle-number states.

Validation [numerical, validate.py, L = 10, both parity sectors, random and product initial states]: ground energies agree with sparse ED of the unrotated Pauli Hamiltonian to 3e-12, entanglement at every cut to 5e-10, ⟨Z_j⟩ to 1.5e-10, the bond operator P_j = X_jX_{j+1} + Y_jY_{j+1} to 6e-11, the boson and fermion one-body matrices on a 6-site window to 3.5e-10, the valley-polarization correlator to 1.2e-8, ⟨ΠZ⟩ to 2e-15, and the energy variance of the final state is below 6e-11. This holds for p153, rat, sym40 without and with the staggered real gauge (C177 PROOF 1.2; observables rotated back), and for the XX control. The two-leg XX control matches ED in the even sector; its odd sector at L = 10 is exactly degenerate (one extra particle on either leg), and a product initial state cannot change the separately conserved leg numbers, so it is run only in the even sector from random states. Against the collaborator's L = 32 rational-bond results, both sectors are lower in energy by 6e-8 and 9e-8, with the same parity splitting (7.24e-5) and S̄ within 3e-4.
