# Central charge of the route-A logical chains at n = 3 and 4

## Claim

Setting. The order-λ² logical chains H = Σ_j K_{j,j+1} of two route-A cells, as spin chains of local dimension 2^n with exact U(1) fermion number: the central-switch cell of `contrib/routeA_compile/` (n = 3, 4) and the folded bus of `contrib/routeA_spatial/` (n = 3). K is the exact second-order bond on the whole 4^n code space at finite conversion amplitude (τ or g), so it carries every generated many-body term, plus the exact third-order bond for the switch cell. The quantity is the Virasoro c.

1. [numerical] The switch cell keeps all 2n Fermi points under self-consistent Hartree-Fock at n = 3 (τ = 3e-4 to 2e-3) and n = 4 (τ = 1e-5, 3e-5), provided small fields break particle-hole symmetry and the common Fermi-sea centre. Its dimensionless couplings are ≤ 1.1e-3 (n = 3) and ≤ 1.3e-4 (n = 4), and both estimators agree with their c = n calibrations. c = n survives.
2. [numerical] The folded bus has c = 3 at g = 0.05. At g = 0.1 weak coupling still predicts 3, but the two estimators disagree (2.9 and 2.4), so c is open. At g = 0.2 Hartree-Fock empties a band for every field tried (c = 2).
3. [numerical] Fermi velocities are 1e-3 to 1e-9 of the logical bandwidth and differ by factors up to 100, and the switch chain shows c = n entanglement scaling only beyond about 2000 cells (n = 3). The measurements therefore use a velocity-equalized proxy with the same Fermi points and dimensionless couplings.
4. The folded bus at n = 4 is not done (Section 5).

## Status

Numerical at the physics standard, on the proxies. Upgrades would be larger χ (queued), extended-precision folded-bus bonds at n = 4, and settling g = 0.1.

## How it was established and checked

- Author: Claude (Opus 5.5) subagent for Jonah Kudler-Flam, 2026-10-28. No Codex call.
- New code, except `vumps.py` and the cell sources in `src/`, copied unchanged.
- Checks. The n-body expansion rebuilds every bond to roundoff. Switch bonds have no number-changing part, folded-bus bonds one of 2e-8, taken as their noise. Ring ED reproduces free-fermion rings to 1e-14 with the spin-chain twist π(N − 1). Nobody else has rerun anything.

## Method

### 1. Bonds and weak coupling

Switch cell (w_D = −0.3, left arm ×1.05). K2 is solved per charge sector with the conversion-free scalar removed analytically (R − R0 = −R C R0). K3 = PVRVRVP is exact, and only V³ on one bond contributes, since a cell touched once keeps a carrier. K3 = 1.32 K2 + R with |R| = 2e-4 |K3|, so third order mostly renormalizes λ², but R holds an imaginary one-cell term 3e-4 λ|T| that competes with the narrowest band. I use λ = 0.05 (0.2 still keeps six Fermi points). The n = 2 phase φ_T = 1e-3 empties two bands at n = 3, so I take φ_T = 1e-5 (n = 3), 1e-7 (n = 4) and a logical field μ = 1e-6, 1e-9. The folded bus uses the `routeA_spatial` working point (D_x = 4, t_d = 1e-5, ε = ½ λ_min(R_s)). With max|T| = 1, two-body norms are 2e-6 to 1e-4 (switch, n = 3), 1e-8 to 9e-8 (n = 4), and 0.07, 0.31, 1.6 (folded bus, g = 0.05, 0.1, 0.2).

Hartree-Fock uses the exact Wick gradient of all bodies. At the Fermi points I compute the antisymmetrized vertex Γ, the forward matrix 1 + Γ_pq/(2π√(v_p v_q)), and for every four-fermion channel the mismatch δ and coupling γ = |Γ|/(2π(v₁v₂v₃v₄)^{1/4}). A marginal channel can act only if γ ln(1/δ) reaches order one.

### 2. Velocities and the proxy

Every low-energy mode hybridizes with the deep flavours, so the narrow bands set the Fermi velocities (switch n = 3: 2e-6 to 2e-4, n = 4: 3e-9 to 2e-7, folded bus: 1e-3 to 9e-3). Direct VUMPS on the folded bus finds no scaling regime at χ ≤ 48 (non-monotonic ξ), and UV/IR ratios of 1e6 and 1e9 rule it out for the switch chains. For the switch chain's free part the exact block-entanglement slope 3 dS/d ln ℓ over doublings from ℓ = 8 to 4096 is 4.1, 4.7, 6.4, 8.3, 6.8, 5.1, 3.7, 3.1, 3.0, and at n = 4 it is still 9.5 at ℓ = 512.

The proxy rescales term by term. Gauge c_j → e^{−iΦj}d_j, rotate each cell to the eigenbasis of the Hermitian part of e^{−iΦ}T_HF (eigenvalues s_a), divide one-body terms by √(s_a s_b) and two-body terms by (s_a s_b s_c s_d)^{1/4}, and normal-order the two-body part to the proxy's own Hartree-Fock state. The one-body map is a congruence, so Fermi points are kept to roundoff (shifted by Φ) while all velocities become 1.7 to 2.7. Dimensionless couplings stay within a factor 1.6 of the true ones (folded bus) and come out 0.9 to 13 times them, mostly larger (switch). The proxies' free parts give c = 3.00 and 4.00 from exact block entanglement.

### 3. Results

FES is VUMPS with a two-cell unit cell, c from S against ln ξ over χ ≥ 16. ΔS is the largest entropy difference from matched XX (n decoupled bands with the Hartree-Fock centres, Fermi momenta and velocities) at equal χ. Ring columns compare ED of the proxy with its free part on the same rings. c_S is the ring-entanglement fit (the exact free c_S falls to 3.03 at L = 128), and c_Cas fits E0(L) − E0,free(L) = AL + B/L with c = n − 12B/(πΣ|v_p|), which cancels the shell oscillations of incommensurate fillings (free part against XX gives 3.02, 3.06, 3.00, 4.00).

| cell | n | coupling | HF Fermi pts | max γ (fwd, other) | FES c chain / free / XX (χ ≤) | ΔS | ring c_S chain / free (L) | c_Cas (L) | c |
| :-- | --: | :-- | --: | :-- | :-- | --: | :-- | :-- | --: |
| fold | 3 | g = 0.05 | 6 | 0.048, 0.027 | 3.41 / 3.02 / 3.17 (32, 96, 96) | 0.05 | 3.357 / 3.352 (7) | 2.98 (4-7) | 3 |
| fold | 3 | g = 0.1 | 6 | 0.16, 0.11 | 2.92 / – / 3.21 (32, –, 64) | 0.31 | 3.39 / – (6) | 2.43 (4-6) | open |
| fold | 3 | g = 0.2 | 4 | 0.05, 0.04 | – | – | – | – | 2 |
| switch | 3 | τ = 3e-4 | 6 | 2e-5, 8e-6 | – | – | – | – | 3 |
| switch | 3 | τ = 1e-3 | 6 | 3e-4, 2e-4 | 2.90 / – / 3.21 (64, –, 96) | 0.004 | 3.364 / 3.364 (7) | 3.00 (4-7) | 3 |
| switch | 3 | τ = 2e-3 | 6 | 1.1e-3, 9e-4 | – | – | – | – | 3 |
| switch | 4 | τ = 1e-5 | 8 | 7e-5, 1e-5 | – | – | – | – | 4 |
| switch | 4 | τ = 3e-5 | 8 | 1.3e-4, 2e-5 | 3.09 / – / 2.96 (48, –, 64) | 0.03 | 4.577 / 4.577 (5) | 4.00 (4-5) | 4 |

Rows with dashes rest on weak coupling, which at couplings ≤ 1e-3 leaves no room for another c. At ξ ≤ 7 FES under-resolves three or four flavours, as matched XX shows (2.96 at n = 4), so the comparison with XX is what counts. For the switch cell at τ = 1e-3 chain and XX entropies agree to 0.004 at every χ from 8 to 64, and ring gaps match the free part to 1 percent. At g = 0.1 the FES slope is consistent with 3, but S and ξ lie well below matched XX, the L = 6 ground state sits at N = 9 where the free part has N = 8, and c_Cas rests on three sizes. A phase change near g = 0.1 and a finite-size artifact would both fit, and I could not tell them apart.

### 4. Fermi momenta, resonances, pair transfer and umklapp

Hartree-Fock Fermi points (k_L, k_R) and centres. Switch n = 3, τ = 1e-3: (−1.599, 1.621), (−1.560, 1.560), (−1.063, 1.041), centres 0.011, 0.0002, −0.011. Switch n = 4, τ = 3e-5: (−1.622, 1.657), (−1.582, 1.590), (−1.538, 1.538), (−1.213, 1.169), centres 0.017, 0.004, 0.0002, −0.022. Folded bus, g = 0.05: (−2.238, 2.233), (−1.560, 1.561), (−0.756, 0.761), centres −0.0027, 0.0001, 0.0026.

Closest distinct-point channels: switch n = 3, mixed exchange at δ = 1.7e-2 (γ = 5e-6) and pair transfer at 2.1e-2. Switch n = 4, a mixed channel at 3e-5 (γ = 2e-5) and pair transfer at 7.9e-3. Folded bus at g = 0.05, pair transfer at 5.0e-3, 5.5e-3, 1.1e-2 (γ ≤ 0.027, γ ln(1/δ) ≤ 0.12), and at g = 0.1 a mixed channel at 1.3e-2 and pair transfer at 2.8e-2 to 5.8e-2 (γ ≤ 0.11, γ ln(1/δ) ≤ 0.31).

Two exact resonances appear when a breaking field is off. At μ = 0 a particle-hole symmetry pairs the outer bands (ν_a + ν_b = 1) and half-fills the middle one, so the mixed umklapp 2k_a + 2k_b = 2π is resonant to 6e-8 (n = 3) and 2e-7 (n = 4). At φ_T = 0 all centres coincide and pair transfer is resonant. With γ near 1e-4 a gap would scale like e^{−1/γ}, invisible but formally fatal, so μ and φ_T are required. Same-point channels such as the 4k_F umklapp of the nearly half-filled middle band have dimension ≥ 3 at K ≈ 1. No pair-transfer or umklapp flow shows up in any row with c = n.

### 5. Not done

At n = 4 the folded bus needs t_d far below 1e-5, since its loop link bypasses the express suppression and makes R_s indefinite (lowest eigenvalue −4e-4). Its scalar part is 2.5e12 times its hopping, so double precision leaves noise of 1e-4 to 8e-4 of |T| against a narrowest band of 1.3e-5. The true chains were not followed past their crossover length.

## Reproducing

numpy and scipy only. Bonds: `make_bonds_compile.py`, `bonds34.py spatial` (`data/bond_*.npz`). Reports: `scan_report.py` (`data/report_*.json`). Chains are `cell:n:param:λ:μ:level` in `models34.py` (full, free, xx, and proxies pfull, pfree, pxx). VUMPS: `run_fes34.py`, `fes_fit.py`. Rings: `ring_ed.py`, `ring_fit.py`. Gaussian references: `gauss_calib.py`, `crossover.py`. The queue was full, so runs went locally and, niced, on the Typhon login node. SLURM array 17517729 (`slurm/fes_tasks.txt`, χ to 256) is queued, writing to `/scratch/lustre/jkudlerflam/routeA_n34_chains/data/fes/`.

## Relation to the rest of the repository

This addresses STATUS "Where help is most useful" item 2 for the logical chains at n = 3 and 4, not the microscopic cells. It supports c = n for the switch cell and for the folded bus at small g, and adds two requirements, explicit particle-hole and centre-breaking fields, and a crossover length of thousands of cells.
