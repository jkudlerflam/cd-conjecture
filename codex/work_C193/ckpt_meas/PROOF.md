This work was done by a Claude Opus 5.5 agent (C193b, 2026-11-10), not by Codex. Astra was not used.

# C193b derivations

Tags [derived], [numerical], [convention] open the paragraph they apply to. Stage one's derivations (`../PROOF.md`) are cited, not repeated.

## 1. What a checkpoint stores and how a correlator is contracted in its gauge

[convention, from `codex/work_C190/z2dmrg.py` and `run190.py`] A checkpoint `outG/ckpt/<tag>.npz` holds M_j[a, s, b] (j = 0..L−1), the even-label counts n0 on every bond, and the sector P. It is written after every sweep by `MPS.save`, and the last write follows the final left sweep, after which the run's final measurements were taken on the same tensors. After a left sweep the tensors j ≥ 1 are right-canonical, Σ_{s,b} M_j[a, s, b] M_j[a', s, b]* = δ_{aa'}, and M_0 carries the norm (the centre). Local basis s = 0 empty (Z = +1), s = 1 occupied, SP = |0⟩⟨1| removes a particle, SM = SP^T adds one. The labels are Z2 quantum numbers only, there is no fermionic sign in the tensors, so contracting the dense arrays is exact.

[convention] Gauge. Runs with gauge 'stag' were done for H' = U H U†, U = Π_j diag(e^{−iθ_j/2}, e^{iθ_j/2}), θ_j = (−1)^j π/4. An operator O of the original Pauli frame has expectation ⟨ψ|O|ψ⟩ = ⟨ψ'|U O U†|ψ'⟩, and U SP_j U† = g_j SP_j, U SM_j U† = g_j* SM_j with g_j = e^{−iθ_j}, while n_j is invariant. The code uses sp_j = g_j SP, sm_j = g_j* SM, as z2dmrg.measure does (g_j = 1 without gauge).

[derived] Contraction. Let E_k[a', a] be the identity-transfer of sites 0..k−1 (bra index first), E_0 = 1. Because every site to the right of a contracted block is right-canonical, the right part closes to δ, and for i < j

⟨O_i X_{i+1} ... X_{j−1} O'_j⟩ = close(T_j, M_j, O'_j),  T_{k+1} = Σ M_k* (X_k M_k) T_k,  T_{i+1} = Σ M_i* (O_i M_i) E_i,

with close(T, M, O) = Σ_{a',s',b,a,s} M[a', s', b]* O[s', s] T[a', a] M[a, s, b]. The same E_k gives the Schmidt weights on bond k: E_k is the reduced density matrix of sites 0..k−1 written in the orthonormal basis of right states on bond k, so its eigenvalues are the Schmidt weights and S_k = −Σ p ln p. All row propagators alive on bond k are stacked and advanced by one batched pair of matrix products per site, X = E·M_k (one GEMM) and E' = M_k†·(O X) (one GEMM), with the close as a dot product with M_k*.

[derived] Pair two-point function with a right block. For the pair at (j, j + d) define on bond j

R_j^{(d)}[a', a] = Σ M_j*...M_{j+d}* (sp_j ⊗ 1 ⊗ ... ⊗ sp_{j+d}) M_j...M_{j+d},  closed on the right by δ,

built right to left in d + 1 steps. Then ⟨A_{left} sp_j sp_{j+d}⟩ = Σ_{a',a} T_j[a', a] R_j^{(d)}[a', a] for any operator string A_left ending before site j, and the one-point piece ⟨sp_j sp_{j+d}⟩ = Σ E_j R_j^{(d)}.

[derived] Energy. E_mpo is the left MPO environment of z2dmrg_xp.env_left (block-sparse, the DMRG's own code) run over all sites; E_phys = E_mpo + E_off with the run's stored E_off. For runs whose output has no variance step (L = 256 and the controls), the stored energy is the last sweep's Davidson value, which equals ⟨H⟩ of the stored state up to the final residual.

## 2. Which qubit operator is the inter-valley pair

[derived] Effective model (C173 eq. 3.2, C176). In the blockaded subspace the qubit chain at second order in Δ/V is the N-conserving model of hard-core particles with no two adjacent, H_eff = P[Σ (w c†_{j+1}c_j + t₂ c†_{j+2}c_j + h.c.) + V₂ n_j n_{j+2} + V₃ n_j n_{j+3}]P. It is obtained from H by a Schrieffer–Wolff rotation, H_eff = P e^{S} H e^{−S} P, with S = (Δ/V) Σ_j (c†_{j+1}c†_j − c_j c_{j+1}) + O((Δ/V)²) (C190 PROOF §2), which removes the pair-creation term Δ(c†c† + h.c.) between the blockaded and the doubly-occupied-bond sectors. An effective-model operator A_eff is represented on the qubit chain by A = e^{−S} A_eff e^{S} = A_eff − [S, A_eff] + O((Δ/V)²).

[derived] C176 measures the pair c_j c_{j+d} at d = 2 and 3 (d = 1 is excluded by the blockade). Inside the blockaded space the sites strictly between the two members are empty when d ≤ 3 and both are occupied, so the Jordan–Wigner strings of c_j and c_{j+d} cancel to a sign, c_j c_{j+d} = −b_j b_{j+d}, and the two-point function of pairs equals ⟨(b_i b_{i+d})† b_j b_{j+d}⟩ (the sign squares away). With b_j = σ⁺_j (SP_j, removes a particle, the hard-core boson of the qubit chain with no string), the qubit representative is

Π^{(d)}_j = SP_j SP_{j+d}  + δΠ,   δΠ = −[S, SP_j SP_{j+d}] = O(Δ/V).

The correction. [S, b_j b_{j+d}] with S ⊃ (Δ/V) c†_{k+1}c†_k is nonzero only for k + 1 or k ∈ {j, j + d}, and each term removes one particle at j or j + d and adds one at a neighbour, so δΠ is a sum of terms (Δ/V) b_{j or j+d} b†_{neighbour} (number-conserving, local). These carry momentum near zero and the single-particle quantum numbers of a one-body bilinear, not of a pair, and their contribution to ⟨Π†Π⟩ at separation r enters either with amplitude (Δ/V)² ≈ 6e-4 (V/Δ = 40) or as a cross term (Δ/V)⟨b b · b†b⟩, which vanishes in the N-conserving effective state (it changes N by 2). So at leading order the two-point function of the effective pair is that of SP_j SP_{j+d}, with the same exponent, and the operator measured is

D_d(i₀, j) = ⟨(b_{i₀} b_{i₀+d})† b_j b_{j+d}⟩ = ⟨sm_{i₀} sm_{i₀+d} sp_j sp_{j+d}⟩,  d = 2, 3.

[derived] Connected part. The qubit chain conserves only parity, so ⟨b_j b_{j+d}⟩ need not vanish. By the argument above it is ⟨−[S, b b]⟩_eff = O(Δ/V) times a local one-body expectation, a constant of order 1e-3 for d = 2 and smaller for d = 3. Its square, about 1e-6, is subtracted: D_conn = D − ⟨sm sm⟩⟨sp sp⟩, both one-point pieces measured in the same pass.

[derived, C191 frame, stage one PROOF §2] Momentum content. With b_j ≈ e^{ik_A j}ψ_A + e^{ik_B j}ψ_B, b_j b_{j+d} ⊃ e^{i(k_A+k_B)j}(e^{ik_B d} + e^{ik_A d}) ψ_Aψ_B + same-valley terms at 2k_A, 2k_B. On the symmetric line k_A + k_B = π, so the inter-valley pair e^{i√2θ_c} (exponent 1/K_c) sits at Q = π with amplitude |e^{ik_B d} + e^{ik_A d}| = 2|cos((k_B − k_A)d/2)|. With k_B − k_A ≈ 2.66 this is 2|cos 2.66| ≈ 1.77 at d = 2 and 2|cos 3.99| ≈ 1.32 at d = 3, both nonzero. The same-valley pairs e^{2iθ_A} have exponent 1/K_c + 1/K_s and are subleading. Off the line the pair momentum moves to Q = k_A + k_B, which EP finds from the data.

[derived] Two-leg control. With the legs on even and odd sites, d = 3 joins the two legs, b_A b_B, exponent η_A + η_B = 1/2 + 1/2 = 1 at momentum 0. d = 2 is a same-leg pair e^{2iθ}, whose exponent is 2² = 4 times that of e^{iθ}, 4 · 1/2 = 2. On the single XX chain (K = 1) every d gives e^{2iθ} with exponent 2. These are the exact values used in the calibration.

## 3. Open-chain forms and estimators

[derived, stage one PROOF §3] For a vertex with Neumann θ at both ends, ⟨e^{iθ(x)}e^{−iθ(y)}⟩ ∝ F(x, y)^{−η}, F = f(|x − y|) f(x + y) / sqrt(f(2x) f(2y)), f(u) = sin(πu/(2(L + 1))), x = site + 1. For pairs symmetric about the centre F = tan(πr/(2(L + 1))). For the pair function the two operators sit at the pair centres x₁ = i₀ + 1 + d/2 and x₂ = x₁ + r.

[convention] EB (η_b). g(r) = mean of ⟨sm_i sp_j⟩ over the pairs with j − i = r and |i + j − (L − 1)| ≤ 8, so the centres are within 4 sites of the chain centre (stage one averaged over the same 9 centres). The one-body function is computed for every i < j with i + j ≤ L + 8, so r runs out to L − 1. Demodulated at each of the two largest Fourier peaks (the two valleys), refined by the phase slope, low-passed by four binomial passes, and ln|h| fit linearly in ln F on ranges given as fractions of L. The effective model's stored window (sites L/4 to 3L/4) gives the same g(r) out to r ≈ L/2 − 9.

[convention] EP (pair). D_conn rows from the four reference sites i₀ = L/2 − 3..L/2 to j = i₀ + d + 1..L − 1 − d (C176's geometry), the momentum Q from the tapered Fourier transform of the mean row, refined by the phase slope, each row demodulated and low-passed, all rows fit jointly to F(x₁, x₂)^{−η}.

[convention] ES (K_c). C(r) = mean over the reference rows of ⟨n_{i₀}n_{i₀+r}⟩ − ⟨n⟩⟨n⟩, S(q) = 2 Σ_{1≤r≤R} cos²(πr/2R) C(r) cos qr, fit a + bq + cq² on 16 values of q in [lo, hi]·π/R, K_c = π b (C192 sq2.py, with R = ρL).

## 4. Normalization of the structure factor and of K_s

[derived] With n_ν = −∂φ_ν/π per component and ⟨∂φ_ν(x)∂φ_ν(0)⟩ = −K_ν/(2x²), a single component has ⟨δn δn⟩ = −K/(2π²x²) and S(q) = K|q|/(2π). Two components with φ_c = (φ_A + φ_B)/√2 have n = −√2∂φ_c/π, ⟨δn δn⟩ = −K_c/(π²x²) and S(q) = K_c|q|/π. So π·slope reads K_c for a two-component liquid (two-leg XX: 1) and K/2 for one component (XX chain: 1/2). On the two-leg control in chain coordinates C(r) vanishes at odd r and is the leg correlator at leg distance r/2 at even r, −2K/(π²r²), whose smooth part −K/(π²r²) gives the same slope (the staggered part sits at q = π).

[derived, stage one PROOF §2 and §4] η_b = 1/(4K_c) + 1/(4K_s) and η_pair = 1/K_c, so K_s = 1/(4η_b − 1/K_c), with 1/K_c from the pair exponent or from the structure factor, and δK_s = K_s² sqrt((4δη_b)² + δ(1/K_c)²).

## 5. Drift and the Kosterlitz–Thouless bound

[derived, stage one PROOF §5] dK_s/dl = −K_s²y₁²/2, dy₁/dl = (2 − 2K_s)y₁; near K_s = 1 with y_∥ = 2(K_s − 1) this is dy_∥/dl = −y₁², dy₁/dl = −y_∥y₁, gapless iff y_∥ ≥ |y₁|. A change ΔK_s over Δl gives y₁² ≈ −2ΔK_s/(K̄_s²Δl). Here the ranges scale with L, so 128 → 256 compares the same r/L and is a change of the RG scale by Δl = ln 2, which stage one could not do. The bound is drift.py's (central, and with ΔK − σ for the upper limit on |y₁|; y_∥ from K_s(256) lowered by its error).

## 6. Validation [numerical]

`val193.py` (numpy, L = 10 and 12; logs val193_L10.log, val193_L12.log). Cases p153 (complex bond), sym40 with the staggered gauge, the XX chain (both sectors) and the two-leg control (even sector, the odd sector is degenerate). Three comparisons per case and sector, for dens, S at every cut, all of Gb, Cnn, Dpair, ⟨sm sm⟩, ⟨sp sp⟩ and the energy.
- measure_state on the DMRG state against the dense vector of the same MPS: worst 3.7e-15 (contraction exact).
- The ED ground vector written as an exact right-canonical Z2 MPS in the gauge frame, through measure_state, against ED with the plain unrotated operators: worst 3.2e-15 in every correlator, 1.9e-13 in E (operators, gauge phases and contraction exact).
- The DMRG state itself against ED: up to 2.6e-8 in Gb, 7e-9 in the pair pieces, 1e-11 in E for the candidate bonds (V/Δ = 40 makes the small-L DMRG state converge only to this level even at residual 1e-14), 1e-15 for the controls. This is the state, not the measurement.

On stored states (each output's `check` field): energy against the stored value, S̄₆ and the whole entropy profile against the run's final entropies, and the stored 48-site Gb against the new one-body function.
