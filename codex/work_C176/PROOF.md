# C176 derivations and code design

This run was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex.

Standard: theoretical physics. Tags [proved], [derived], [numerical], [heuristic], [assumed] open the paragraph they apply to. The model is C173 PROOF eq. (3.2), in units |t₂| = 1,

    H = P Σ_j [ (w c†_{j+1} c_j + h.c.) − (c†_{j+2} c_j + h.c.) − 2 n_j n_{j+2} − (2/3) n_j n_{j+3} ] P − μ N,     (0.1)

with P the projector on configurations with no two adjacent particles and w = |w| e^{iφ}.

## 1. A blockade-exact, N-conserving DMRG

### 1.1 Why the obvious label scheme fails

[proved] The natural idea is to label every MPS bond by (N_left, f), with f the occupation of the site just left of the bond, and to require f s = 0 when a site of occupation s is attached. On the left-canonical side this is fine: a left-block state with label f has n_{i−1} = f, so states with different labels are orthogonal. On the right-canonical side it is not. A right-block state carrying the label f = 1 only promises that its first site is empty, and a right-block state with f = 0 may also have its first site empty, so right states with different f overlap. Then the centre matrix is not a Schmidt matrix, the per-block singular values are not Schmidt values, and the eigenproblem of a two-site update is a generalized one. Concretely, the reduced density matrix of a constrained state does not commute with the occupation of the last site of the block, ρ_L between left states with n_{i−1} = 0 and 1 is Σ_r ψ_{l r} ψ*_{l′ r}, which is nonzero through right states r with n_i = 0.

### 1.2 The scheme used

[proved] fibdmrg.py labels left-canonical tensors A[i] by left-physical labels (N_left, n_{i−1}) on bond i, and right-canonical tensors B[i] by right-physical labels (N_right, n_i). Each label is a quantum number of the block it describes (the particle number of the block and the occupation of its site next to the cut), so each basis is orthonormal across labels, and orthonormal within a label by the canonical condition. The two-site centre θ(j, j+1) is indexed by (q_l, s₁, s₂, q_r), with q_l left-physical on bond j and q_r right-physical on bond j+2. Its allowed blocks obey

    f(q_l) s₁ = 0,  s₁ s₂ = 0,  s₂ g(q_r) = 0,  N(q_l) + s₁ + s₂ + N(q_r) = N_tot,

so the blockade and N hold as structural zeros of every tensor, with no penalty, at every bond and every sweep. The effective Hamiltonian in this orthonormal basis is P H P restricted to the variational space (the bra and the ket both lie in the constrained space, and the MPO is the plain H, so its matrix elements between constrained states are those of P H P). The NNN hop needs no explicit string, since site j+1 is empty whenever j is occupied, and particles never pass, so hard-core bosons and fermions give the same matrix.

Truncation. Moving right, θ is split per new left-physical label (N_left + s₁, s₁), and moving left per new right-physical label (N_right + s₂, s₂). The singular values of all groups are pooled, the χ largest are kept, and the discarded weight is reported. This is the reduced density matrix restricted to its label-diagonal blocks, so it is exact bookkeeping of the discarded weight but slightly suboptimal relative to an unconstrained SVD (the true Schmidt vectors mix n_{i−1} = 0 and 1, so a label-respecting basis needs at most twice the bond dimension for the same error).

Entanglement. After each split the centre is a matrix between an orthonormal left-physical basis and an orthonormal right-physical basis. Its full SVD within each N sector (which mixes the n = 0 and n = 1 labels) gives the exact Schmidt values of the represented state. S(l) is recorded from these.

MPO. Nine channels: start, b and b† placed one and two sites back (NN and NNN hops), n placed one, two and three sites back (V₂ and V₃), and end. The chemical potential is not in the MPO, since N is fixed. A gauge option runs the model in the frame d_j = (−i)^j c_j, where φ = π/2 is real (w → −iw = |w|, t₂ → +1), and measurements are converted back to the c frame.

### 1.3 Validation

[numerical, validate.py, validate2.py] Against C173's exact diagonalization (eff_ed.py, unmodified), open chains, ten cases with L = 12 to 24 including φ = π/2 in both gauges, φ = 1.52, 1.3, 0.4, a dense case N = 9 at L = 22, and both controls:
- ground-state energy to 3e-12 or better,
- S(l) at every cut to 6e-11,
- the density profile to 5e-11,
- the full fermion one-body density matrix ⟨c†_i c_j⟩ (with Jordan-Wigner string) to 3e-9,
- ⟨P_j⟩, the valley-polarization bond operator, to 3e-11,
- the connected ⟨n n⟩ correlation to 2e-10, and ⟨P_i P_j⟩ to 1e-10 for |i − j| ≥ 3. At |i − j| = 2 the ED reference projects the intermediate state P_i|ψ⟩ onto the blockaded basis and the MPS contraction does not, so they differ there by construction, and the analysis uses |i − j| ≥ 3.
- ⟨n_j n_{j+1}⟩ = 0 exactly (structural).
- The degenerate odd-N c = 2 control (two legs with 2 and 3 particles, gap 5e-15) agrees in energy but not in the state, as expected. At even N it agrees in S and ⟨n n⟩ to 1e-13. (Its fermion one-body matrix differs from eff_ed's because the control runs hard-core bosons whose NNN hop passes, while eff_ed's control is fermions. Only its entanglement is used.)
- Pair correlations (validate3.py). D(r) = ⟨(c_i c_{i+d})† c_j c_{j+d}⟩, d = 2, 3, against ED with explicit Jordan-Wigner signs, in both gauges: agreement to 7e-10. The reduction used in the code is [proved]: c_j c_{j+d} = b_j [Π_{j≤l<j+d}(1 − 2n_l)] b_{j+d} = −b_j b_{j+d} on the blockaded space, because b_j (1 − 2n_j) = −b_j and the sites strictly between the two members are empty when both are occupied (d = 2, 3). The strings of the two pairs cancel, and in the d frame c_j c_{j+d} = i^{2j+d} d_j d_{j+d}, so D_c = (−1)^{j−i} D_d.
- Ring ED, bulk energy (ringcheck.py, ringcheck.out, ringcheck_open.out) [numerical]. At n = 1/6 the ring energy per site is −0.3700, −0.3749, −0.3722, −0.3736 at L = 18, 24, 30, 36 (lower fermion boundary condition; it oscillates with the valley shell filling). The open-chain bulk energy density from DMRG, e = [E(72, 12) − E(48, 8)]/24, is −0.37327. These agree at the level of the ring's shell oscillation (about 1e-3), a check of the MPO coefficients independent of the open-chain ED.

### 1.4 Cost and convergence

[numerical, local/t128.log] The block-sparse contraction makes a sweep at L = 128, χ = 256 take about 3 minutes on two cores. The truncation falls roughly as χ^{-3} (2e-6, 1e-6, 4e-7 at χ = 128, 192, 256, L = 128, n = 3/16), so 1e-7 needs χ ≈ 400 and 1e-8 needs χ ≈ 900 at L = 128. Sweep convergence at small χ is slow, because the dilute particles have to redistribute along the chain, so every schedule starts with 16 sweeps at χ ≤ 64 before the bond dimension is raised.

## 2. Densities and the chemical potential

[numerical, scan_L48_phipi2.out] In the exact model the one-particle band (C173 eq. 3.3) at φ = π/2, |w| = 1 is 2 sin k − 2 cos 2k, with minimum −9/4 at sin k = −1/4 [proved, algebra] and curvature 15/2 there. At L = 48 the open-chain μ(N) = E(N) − E(N − 1) runs from −2.27 to −2.06 for N = 2 to 14, with a strong even-odd alternation (even N lower by about 0.05). On open chains the edges are depleted (mid-chain density 0.204 at L = 128, N = 24, against N/L = 0.1875), so N/L is not the bulk density. C175's ring ED gives n = 0.167 at μ = −2.2. μ = −2.4 lies below the band bottom −2.25 and below every two-body bound level C175 reports (E_b = 0.0195), so it is the empty chain in the exact model [derived]. Production therefore fixes even N at n = N/L = 1/8, 5/32, 1/6, 3/16 and 1/4, which brackets the bulk density 0.17 of C173's reading, and the parity gap uses N ± 1.

## 3. Estimators and the density drift

Definitions are the registered ones (RESULT.md). Two points need derivation.

Density drift [derived, with the numerical input δ ≈ 2]. At fixed N/L the open chain has an edge deficit δ of order one particle per edge, so the mid-chain density is n_mid(L) = N/L + δ/L + o(1/L). In a Luttinger liquid the non-universal constant of S depends on the density through the Fermi momenta, S = (c/6) ln[(2L/π) sin(πl/L)] + s₀(n). For two-component fermions s₀ depends on n as (1/6) Σ_valleys ln k_F (the free-fermion form, assumed to hold up to O(1) changes of the coefficient). Then a two-size estimate at fixed N/L is biased by

    δc ≈ 6 [s₀(n_mid(L₂)) − s₀(n_mid(L₁))]/ln(L₂/L₁) ≈ c [ln n_mid(L₂) − ln n_mid(L₁)]/ln(L₂/L₁) < 0.

Because n_mid falls toward N/L as L grows, the bias is downward. With δ ≈ 2 it is about −0.12 at n = 3/16 for 128 → 256, and −0.2 for 48 → 64 at n = 1/8. So a fixed-N/L reading is not inflated by the drift, and a registered pass is conservative on this count. The secondary estimator in analyze.py (density_corrected) removes the drift by interpolating S̄_W(L) linearly in ln n_mid across the φ = π/2 series at each L and forming c_two at a fixed bulk density. It is not registered and is reported only alongside.

Friedel window [heuristic]. The intra-valley 2k_F oscillation of S has period π/k_F = 2/n sites (two valleys, each holding n/2). That is 16 sites at n = 1/8 and 11 at n = 3/16. A window average over 2W + 1 cuts with 2W + 1 well below that period does not average the oscillation out, which is why W = 1 (C173's estimator) swings by ±1 in c at L ≤ 64 and n = 1/8 while W = 6 and 12 agree. In the job's output W = 6 and 12 are the meaningful ones at low density. W = 1 is kept for continuity with C173.

## 4. Preliminary numerics (L = 48, 64)

See RESULT.md claims 5 to 7, with sources local/pre/summary.txt and local/pre/summary_parity_density.txt. Two checks [numerical]. The N + 1 run at n = 1/8 and the N − 1 run at n = 1/6 at L = 48 are the same problem (N = 7) with the same seed, and they give the same energy, −15.543402010857, so the runs are deterministic (a reproducibility check, not an independent one). The final energies lie within 1e-8 of a linear extrapolation in the truncation at L = 48. At L = 64 they lie within 2e-6 to 4e-6 (3.5e-6 for n = 3/16, 1.9e-6 for n = 5/32 and the c = 2 control), and the parity gaps change by less than 1e-6 between χ = 192 and 256 (analyze.py diagnostics).
