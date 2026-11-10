This was written by a Claude Opus 5.5 agent (C196, 2026-11-10), not by Codex.

# C196 methods and derivations

Tags [derived], [numerical], [heuristic] open the paragraph they apply to.

## 1. Model, labels, unit cell

The chain is C190 PROOF eq. (1.1) without edges. The bulk MPO tensor is site 4 of C190's 8-site MPO (models.candidate, z2dmrg.MPO), with the boundary vectors start (channel 0) on the left and end (channel 2) on the right. [derived] On any open chain this reproduces z2dmrg.MPO exactly. The edge rows and columns that MPO zeroes are never reached from those boundary vectors, and the Pauli bonds have no Z⊗1 part that would change the edge field. Energies: e_phys = e_MPO − shift, with shift = −(J_z + h_z), as in C190 (E_off = const − L·shift).

Bond labels are the parity of the particle number to the left (C190 PROOF 5). In an infinite or half-infinite setting only relative labels matter. With total parity P = 0 on every finite window, the parity to the left of a bond equals the parity to its right, so the left and right Schmidt bases of a bond carry the same label sets.

Unit cell and valleys [derived]. The two valleys sit at k and π − k. A one-site cell imposes translation invariance, which the Hamiltonian has, and identifies no momenta. A two-site cell folds k onto k + π, so the valley at π − k folds onto −k. The two valleys coincide only when k ≡ π/2 and differ by π only when 2k ≡ 0, and neither holds at the candidates (cos k₀ ≈ −1/4 at |w| = |t₂|). So both cells keep the valleys distinct. The ground state is not expected to break translation. A two-site cell can represent a period-two state, a one-site cell cannot.

## 2. Infinite DMRG on the C190 solver, and why it has no fixed point here

[derived] Two-site McCulloch iDMRG. After step t the state reads … A_t Λ_t B_t …, and the next two sites go between A_t and B_t with the guess θ_{t+1}(s₁, s₂) = Λ_t B_t(s₁) Λ_{t−1}^{-1} A_t(s₂) Λ_t. Its left index is the left basis of bond t, B_t maps that to the right basis of bond t − 1, Λ_{t−1}^{-1} maps that to the left basis of bond t − 1, and A_t returns to the left basis of bond t. So every contraction is between matching bases, and the labels satisfy q_c = q_a ⊕ s₁ ⊕ s₂, the block structure of C190's TwoSite. The local problem is C190's TwoSite with the kron4 preconditioner and Davidson (residual 1e-9, at most 300 iterations, σ = 0.01). The cell M(s₁, s₂) = A_t(s₁) Λ_t B_t(s₂) Λ_{t−1}^{-1} closes on the left basis of bond t − 1, and its transfer map X ↦ Σ M†XM conserves r = q_bra ⊕ q_ket. That gives ξ_r = −2/ln(|λ_r|/|λ₀|) per parity sector. [numerical] Untruncated, it reproduces open-chain ED to 5e-14 (L ≤ 10), and on the TFIM it reaches the exact energy density to 1e-13 with λ₀ = 1 to 1e-15.

[derived] Why it fails at incommensurate filling. Suppose H conserves N (XX) or conserves it up to virtual pairs of amplitude Δ/V (the candidate). Then the ground state of each finite window is, up to O((Δ/V)²), an eigenstate of N, and so are its Schmidt vectors. Each step adds two sites, so N must grow by 2n per step on average, and exact parity with P fixed forces even increments. At 2n = 0.39 (XX control) or 0.43 (p153c) the sequence of windows has no fixed point, only a quasi-periodic orbit of particle insertions. [numerical] This is what the runs show. For XX at χ = 32 the energy density wanders by ±5e-3 and the density between 0.184 and 0.199 over 400 steps. For p153c at χ = 32 the cell develops site densities 0.08 and 0.34 that swap every step, and e drifts away from the C190 bulk value. C189 avoided this with U(1) labels at a commensurate density p/q. The qubit chain has no exact U(1), so that cure is not available.

## 3. One-site VUMPS with the C190 preconditioner

[derived] Uniform MPS (A_L, C, A_R, A_C) with one bond basis, block-diagonal C. Left environment L[a′, w, a] (bra, channel, ket), computed channel by channel in the order 0, 1, 3, 4, 5, 6, 2 of C190's lower-triangular MPO. L₀ = 1. Every channel except 2 is a finite sum of single transfers T_O(X) = Σ O_{s′s} A_{s′}^† X A_s. Channel 2 solves (1 − T + |1)(ρ|) L₂ = Y − tr(Yρ) 1 with ρ = CC† (Σ_s A_L ρ A_L† = ρ), energy per site e = tr(Yρ), by restarted GMRES. The right environment follows in the order 2, 1, 5, 6, 3, 4, 0. For the map X ↦ Σ O_{s′s} B̄_{s′} X B_sᵀ the invariant functional is X ↦ Σ_{c′c} X_{c′c}(C†C)_{c′c}, since Σ_s B_sᵀ (C†C)ᵀ B̄_s = (C†C)ᵀ.

A_C is the lowest eigenvector of H_AC, found with C190's Davidson. The preconditioner is the single-site analogue of C190 PROOF eq. (5.1). On the block (q_a, s, q_c),

    H₀ = H_L^{(s)} ⊗ 1 + 1 ⊗ H_R^{(s)} − h₁(s),
    H_L^{(s)} = Σ_{w ∈ {0,1,2}} L[:, w, :] W[w, end, s, s],   H_R^{(s)} = Σ_{u ∈ {0,1,2}} W[start, u, s, s] R[:, u, :].

It keeps every term diagonal in s (the field, the ZZ couplings to both neighbours, both half-chain Hamiltonians), drops only the terms that flip s, and subtracts h₁ once because both sums contain it. C is the lowest eigenvector of H_C. With converged environments H_C = P†H_AC P for the isometry P: C ↦ A_L C. So H_C is preconditioned by P†(H₀ − θ + σ)^{-1}P, which places the doubly occupied bond at its cost V. A Kronecker sum L₂ ⊗ 1 + 1 ⊗ R₀ does not, and with it the C solve hit its iteration cap. The update is A_L = U_{A_C} U_C†, A_R = U_C† U_{A_C}ʳ (polar factors). The gradient norm is ε = max(|A_C − A_L C|, |A_C − C A_R|).

[numerical] TFIM (D = 16): ε = 3e-11, e is 4e-10 above exact. XX chain (D = 16 to 64): ε ≤ 1e-9, e error 5e-5 to 1e-6.

[derived, numerical] Symmetry breaking at finite D. Without a U(1) label, the variational uniform MPS of the XX chain breaks U(1) at finite D. With exact Z2 labels it appears as an equal-weight cat of the two components related by the π rotation. The transfer matrix is then two copies of the same spectrum, one in each parity sector, with |λ| = 1 in the odd sector, and S carries an extra ln 2. Within a cat the "per sector" correlation lengths are the same numbers, so a sector-resolved statement needs the spectrum of one component. [numerical] c from S against ln ξ_even, between consecutive D = 16, 24, 32, 48, 64, reads 1.015, 1.037, 0.975, 0.985. The constant ln 2 drops out of the slope.

## 4. The candidate: what was tried at D = 32 (local, numpy)

[numerical] The fixed field of p153c (h_z = −19.97589125) was used, with a random Z2 start unless stated.
- Plain VUMPS. ε stays at 3e-3 to 1.6e-2 over 200 iterations, e jumps by up to 5e-3, the C solve with the Kronecker-sum preconditioner hits 300 iterations, and the state sits at the dressed vacuum (n ≈ 0.002).
- C from the projection A_L†A_C (no C eigenproblem). Unstable, ε ≈ 1.
- Two-site centre (C190 TwoSite on the infinite environments, A_C and C read off by projection). Unstable, ε ≈ 1.
- C preconditioned through H₀ as above. Every C solve converges (35 to 117 matvecs), but ε stays near 1e-2.
- Damped updates (0.3). ε falls to 3e-4 to 1e-3 and e stays within ±3e-4, but the state is still near the vacuum (n = 0.0027) and does not converge.
- Annealing from a dense liquid (μ_eff = −1.0) down to the target: see RESULT.md.

[heuristic] Reading. The A_C eigenproblem is solved accurately at every step (20 to 30 matvecs), and the environments to 1e-12. What fails is the outer fixed-point iteration. The vacuum and the liquid differ by about 1e-4 per site, which is |t₂|·n, against level spacings of order V in the same local problem. The lowest A_C eigenvector then changes character between steps. This is the scale separation the spec names. It is not cured by the local preconditioner, because the local problem is not where the difficulty lies.
