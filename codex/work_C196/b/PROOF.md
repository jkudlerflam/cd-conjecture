This was written by a Claude Opus 5.5 agent (C196b, 2026-11-10), not by Codex.

# C196b methods: embedding, labels, iteration, estimator

Tags [derived], [numerical], [heuristic] open the paragraph they apply to.

## 1. The embedding of the effective model

[derived, checked numerically] In z2dmrg's basis the local state s = 1 has Z = −1 and is the flipped spin, the particle (C190's product starts and density profiles use n = (1 − Z)/2). The bond matrix of models.candidate, in the basis |s₁s₂⟩, has
- ⟨01|h|10⟩ = J_x + J_y + 2iD = w, the amplitude for a particle to move one site to the right (w = 0.0125i at sym40, e^{1.53i}/80 at p153),
- ⟨11|h|00⟩ = J_x − J_y = Δ = 1, real,
- a cost V = 4J_z = 40 per adjacent occupied pair (relative to two separated particles).
A particle at j moves to j + 2 at second order only through |1 1 1⟩ (pair created on j + 1, j + 2, then the pair j, j + 1 removed), whose cost is 2V. The amplitude is Δ·Δ/(−2V) = −|t₂| with |t₂| = Δ²/(2V) = 1/80. In spin language there are no string signs (the blockade keeps site j + 1 empty when the particle at j hops to j + 2). So the effective model of C173 eq. (3.2) in spin language is C176's fibdmrg model in gauge 0 (hard-core bosons), with w_eff = w/|t₂| (= i at sym40, e^{1.53i} at p153), t₂ = −1, V₂ = −2, V₃ = −2/3. No gauge transformation is needed.

The start state of a candidate run is the ground state of this model on an open chain of 2q sites with exactly 2p particles (fibdmrg, χ up to 64), converted to an MPS on the qubit chain. fibdmrg's right-canonical tensors carry right-physical labels (N_right, n) on bond i. With total N = 2p even, the parity of the particles to the left of bond i is (2p − N_right) mod 2 = N_right mod 2, which is the Z2 label of z2dmrg. The conversion orders the blocks of each bond by (Z2 label, N_right, n) with label 0 first and copies each block. The tensors stay right-orthonormal because the conversion only relabels and concatenates blocks. [numerical] The converted 32-site window at sym40 reproduces N = 6 exactly and the effective-model density profile.

The same routine with blockade=False, t₂ = 0, w = −1 (XX chain) or w = 0, t₂ = −1 (two decoupled XX legs, hard-core bosons) gives N-definite starts for the controls.

The embedded state has no virtual pairs. They are local (one adjacent pair, amplitude Δ/V ≈ 1/40) and the qubit-chain iterations add them.

## 2. The q-site-cell infinite DMRG

[derived] McCulloch's iDMRG with a window of two cells, as in C189's idmrg.py, on the C190 qubit-chain solver. At the start of iteration t the window is LP · M[0] ⋯ M[2q−1] · RP, right-canonical with the centre at site 0. Bond 0 and bond 2q carry the basis of Λ_{t−1}, the left Schmidt basis at bond 0 and the right Schmidt basis at bond 2q. The total parity of every finite window is even, so the parity to the left of a bond equals the parity to its right, and the left and right Schmidt bases of a bond carry the same labels index by index (first attempt, PROOF 1). Λ is therefore diagonal and its inverse is taken elementwise (entries below 1e-10 of the largest are dropped).

One iteration: two-site updates to the right over the whole window (2q − 1), then to the left down to the pair (q − 1, q) (q updates). The split of that pair gives Λ_t at bond q. The left cell A[0..q−1] (left-orthonormal) is absorbed into LP, the right cell B[q..2q−1] (right-orthonormal) into RP, and the completed channel of LP (channel 2 of C190's lower-triangular MPO) is shifted by −E·1 with E the last window eigenvalue. If the shifts accumulate to the total energy T_t of the chain before iteration t + 1, the next eigenvalue is T_{t+1} − T_t, the energy of the 2q added sites, so e_t = E/(2q) + e_off with e_off the constant C190 adds per site. The next window is McCulloch's prediction Λ_t B[q] ⋯ B[2q−1] Λ_{t−1}^{-1} A[0] ⋯ A[q−1] Λ_t, right-canonicalized by RQ steps per Z2 block without truncation and normalized.

Each local problem is C190's TwoSite with the (s₁, s₂)-block Kronecker preconditioner and Davidson, residual target 1e-9, at most 300 iterations, σ = 0.01, the settings of the first attempt. Truncation is z2dmrg's rule (global χ, relative cutoff 1e-14 on s²). An optional density-matrix mixer (C190's perturbation through the first attempt's mixer_split) is used only in warm-up iterations.

[numerical, G1] Without truncation and with three extra complete window sweeps per iteration the accumulated eigenvalues equal the open-chain ED ground energies (P = 0) for q = 2, 3, 4 and L ≤ 16 at sym40_n316, p153c and the TFIM, worst deviation 2.7e-11. One pass per iteration is not exact in the first growth steps, because the first prediction is a product across the middle bond and the bonds left of the middle are not revisited after it grows.

## 3. Estimator

The transfer map of two cells, X ↦ Σ M†XM with M = A[0..q−1] Λ_t B[q..2q−1] Λ_{t−1}^{-1}, closes on the left basis of bond 0 (it is the first attempt's two-site map with cells). It conserves r = q_bra ⊕ q_ket, r = 0 for parity-even and r = 1 for parity-odd operators. ξ_r = −2q/ln(|λ_r|/λ₀) in sites. The number per cell comes from the left cell with ρ = Λ_t² at bond q, with generating environments for 1, N and N².

## 4. VUMPS with a multi-site cell (route 2)

[derived] The cell of Lc = 2q sites closing the last iDMRG iteration, M = A[0..q−1] (Λ_t B[q]) B[q+1] ⋯ (B[2q−1] Λ_{t−1}^{-1}), is a uniform MPS of period 2q on the bond-0 basis. Fixed-point QR sweeps around the cell give Y_j M_j = A_L^j Y_{j+1}, LQ sweeps give M_j Z_{j+1} = Z_j A_R^j, and C[j] = Y_j Z_j satisfies A_L^j C[j+1] = C[j] A_R^j. The environments at bond 0 are solved channel by channel (the identity channel exactly, the intermediate channels by two passes around the cell, the completed channel by GMRES with the fixed-point projector and the energy per cell subtracted) and propagated to every bond with C190's env_left and env_right. A_C^j and C[j] are the lowest eigenvectors of the one-site and bond problems (the first attempt's preconditioners, generalized to different left and right bonds), and A_L, A_R follow from the polar rule. [numerical] Converges on the XX control, collapses to the vacuum at the candidate (RESULT).

## 5. Why no exact label can pin the real-particle number

[derived] H changes the spin number by 0 or ±2. A label that equals the number of real particles must give a run of ℓ adjacent flipped spins the value ℓ mod 2 (an isolated virtual pair 11 carries no particle, the intermediate 111 of the second-order hop carries one). The first-order hop w turns 1 0 1 into 1 1 0, which changes that value by 2. The process is virtual (cost V) but it is the first step of the third-order real pair creation, so forbidding it changes the model at order wΔ/V. Any exact label is therefore either the parity (what the code uses) or a change of the Hamiltonian.

## 6. Diagnostic of U(1) breaking

[numerical] On the window of a converged iteration, ⟨b_i b_j⟩ with b = |0⟩⟨1| (two annihilations) and ⟨b_i† b_j⟩ are computed by transfer matrices (local/u1diag.py). In a state symmetric under the real-particle U(1) the first is short-ranged (virtual pairs only); in a pair condensate it is as long-ranged as the second.
