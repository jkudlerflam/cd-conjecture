# Luttinger parameter of the tuned ten-letter demonstrator

## Claim

Class: the nearest-neighbour ten-letter bond of [`../routeA_demonstrator/cell.py`](../routeA_demonstrator/cell.py), unchanged, with the collaborator's Jp/Jz tuning from [`../decisive_attack_20261028/demonstrator/`](../decisive_attack_20261028/demonstrator/README.md), on the full microscopic chain in its conserved period-3 type sector, at logical filling n = 1/3 (off half filling), for λ from 0.0125 to 0.15. The tuning keeps Jp Jz = −2, so the second-order hopping stays at −9/32 for every anisotropy, and varies ρ = Jp/|Jz| = 2, 1, 1/2, which gives the second-order anisotropies g = 27/182, 27/364, 27/728. No four-site (C53) term is used anywhere.

1. [numerical] The microscopic chain has c = 1. Interval entanglement on open chains of M = 24 and 48 cells (72 and 144 sites) gives c between 0.977 and 1.011, equal to that of its logical chain at the same M to 0.002. The same estimator on the bare XXZ chain gives 0.984, 0.989, 0.993 at M = 48, 96, 192, so this is the calibrated finite-size value of c = 1. Finite-entanglement scaling (VUMPS, ξ up to 35 cells) gives 0.96 to 1.00 with fit errors 0.02 to 0.05, against 0.94 to 0.98 for the logical chain on the same ξ window.
2. [numerical] K and v are those of an XXZ chain. At every (λ, g) the microscopic Luttinger parameter equals the Bethe-ansatz K of the XXZ chain whose anisotropy g_x(λ) and hopping are read off the exact two-cell logical bond, to 2 × 10⁻⁴ for λ ≤ 0.05, 1.5 × 10⁻³ at λ = 0.1, and 9 × 10⁻³ at the worst point (λ = 0.15, ρ = 1/2, where g_x = −0.33). That residual scales as λ², and the velocity behaves the same way.
3. [numerical] Measured against the second-order prediction K_BA(g), the microscopic K deviates linearly in λ, and the deviation goes to zero as λ → 0. At ρ = 1/2 it is +0.0092, +0.0190, +0.0419 at λ = 0.0125, 0.025, 0.05. The slope is dK/dΔ × dg_x/dλ, so the whole O(λ) correction is a local shift of the marginal parameter. This is the numerical content of input S in this family: the infrared of the matched retarded action is the XXZ Luttinger liquid at the locally renormalized anisotropy, with no further renormalization of K visible at the 10⁻³ level.
4. [numerical, exact two-cell algebra] The tuning has a definite order of limits. The third-order shift of the anisotropy is g_x − g ≈ −0.41λ, −0.70λ, −1.39λ for ρ = 2, 1, 1/2 (small λ), and its coefficient grows as ρ decreases. At fixed λ the "nearly XX" chain is therefore not nearly XX. At λ = 0.1, ρ = 1/2 the effective anisotropy is −0.15 and K = 1.09. The anisotropy is small only when λ ≪ g/1.4 at ρ = 1/2, which is the order (ratio first, λ last) that the collaborator's README prescribes.

## Status

Every claim is numerical, on finite open chains (up to 144 sites) and on infinite MPS at finite bond dimension. Nothing here is a theorem. The tests are in-family and at one filling, n = 1/3. They assume the frozen P B P type sector, which the bond conserves exactly and which `../routeA_demonstrator/` selects at these couplings, but selection is not retested here. The thermodynamic K and v of the microscopic chain are obtained by combining the microscopic/logical ratio at M ≤ 48 with the logical chain's extrapolation (the finite-size corrections of the two chains agree to the digits shown). When this was written, the M = 48 microscopic runs at λ = 0.0125 and 0.025 (array 17523039) were still running, and the M = 48 run at λ = 0.1, ρ = 1 had not finished its excited state before its time limit. At those points Table 3 uses M ≤ 24, where the microscopic/logical difference is already below 4 × 10⁻⁵ (λ ≤ 0.025) or stable between M = 12 and 24 (λ = 0.1). Rerunning `python3 analyze.py --readme` after copying their JSON files into `data/obc/` updates the tables.

## How it was established and checked

- Author: Claude (Opus 5.5) subagent for Jonah Kudler-Flam, 2026-10-28/29. No Codex. Typhon SLURM arrays (short and medium QOS) for the microscopic runs, logical chains locally.
- The microscopic chain is built from the blocks of `cell.bond()` itself (`chain.py`). Two cells reproduce `cell.inter_cell` to 2 × 10⁻¹⁴ in spectrum, and every term commutes with the charge.
- DMRG (`dmrg.py`) against exact diagonalization in the fixed-charge sector: logical M = 12 and micro M = 3 (32768 states), ground and first excited energies to 10⁻¹⁵ (`test_small.py`).
- Bethe ansatz (`bethe.py`): Δ = 0 gives K = 1 and v = sin(πn) exactly, B → ∞ reproduces K = π/(2(π − γ)) and v = (π/2) sin γ/γ, quadrature converged to 10⁻¹⁵, and the compressibility identity dn/dh = K/(πv) holds to 10⁻⁸.
- Calibration of every open-chain estimator on the bare XXZ chain against the Bethe ansatz (Table 1): K to 3 × 10⁻⁵, v to 1 × 10⁻³.
- One numerical pitfall is recorded because it matters. The logical scale is 10⁻⁴ to 10⁻² of the cell gap, so the two-site Lanczos needs about 40 Krylov vectors per update (60 to 80 for λ ≤ 0.025). With 10 vectors the energy creeps down by 10⁻⁸ per sweep for tens of sweeps and looks converged when it is not. All microscopic runs here use the larger Krylov spaces, start from the encoded logical ground state, and end with discarded weights between 10⁻¹⁰ (M = 48) and 10⁻¹⁸ (M = 12).
- Not independently rerun by anyone else.

## Method and results

### Model and tuning

Sites are P (ℂ²), B (ℂ⁸), P, repeated, which is the period-3 type sector of the d = 10 chain. Each nearest-neighbour term is the corresponding block of the original bond, with D₀ = D₁ = 1, κ = 0.6, BB repulsion 2, μ = 0.55 and the port coupling λ[Jp(σ⁺σ⁻ + σ⁻σ⁺) + Jz ZZ]. The second-order logical chain is λ²[Jxy(σ⁺σ⁻ + h.c.) + Jzz ZZ], which after the staggered rotation is the XXZ chain J Σ(SˣSˣ + SʸSʸ + Δ SᶻSᶻ) with J₂ = 2|Jxy|λ² and Δ = g = 2Jzz/|Jxy| (C47 and the collaborator's formulas).

| ρ = Jp/\|Jz\| | Jp | Jz | Jxy | g |
|---|---|---|---|---|
| 2 (original) | 2 | −1 | −9/32 | 27/182 = 0.1484 |
| 1 | √2 | −√2 | −9/32 | 27/364 = 0.0742 |
| 1/2 | 1 | −2 | −9/32 | 27/728 = 0.0371 |

The λ-dressed logical chain ("logical_x") uses the exact two-cell des Cloizeaux bond of `../routeA_demonstrator/effective.py`, which is XXZ plus a constant at every λ (other Pauli terms below 10⁻¹⁶). Its anisotropy g_x and hopping J_x are

| λ | g_x (ρ = 2) | g_x (ρ = 1) | g_x (ρ = 1/2) | J_x/J₂ (ρ = 2, 1, 1/2) |
|---|---|---|---|---|
| 0.0125 | +0.1434 | +0.0654 | +0.0197 | 0.961, 0.965, 0.957 |
| 0.025 | +0.1391 | +0.0570 | +0.0015 | 0.925, 0.932, 0.917 |
| 0.05 | +0.1322 | +0.0408 | −0.0391 | 0.860, 0.871, 0.844 |
| 0.10 | +0.1225 | +0.0067 | −0.1510 | 0.749, 0.765, 0.716 |
| 0.15 | +0.1151 | −0.0353 | −0.3315 | 0.658, 0.676, 0.608 |

The charge Q is conserved exactly, and on the code it equals Σ Z_c/2, so the logical filling of the microscopic chain is exact, with no O(λ) dressing. Filling n = 1/3 means Q = M/6 on open chains, selected by a penalty w(Q − M/6)² in the MPO (eigenstates of that sector are unchanged). VUMPS uses a one-site field −h_f Z_c on the B letter instead, with h_f from the Bethe ansatz for the logical_x chain, and the filling reached is listed.

### Estimators

All lengths are in cells.

- K_gap = E₁/Δ_c on open chains, with E₁ the lowest neutral gap (≈ πv/L) and Δ_c = E(Q₀+1) + E(Q₀−1) − 2E(Q₀) (≈ πv/(KL)). The ratio cancels the effective length to leading order. K_inf is a fit K + a/M + b/M² through M = 12, 24, 48.
- v from 1/E₁ = (M + δ)/(πv), linear in M, through the two largest sizes. δ comes out at 0.73 in every chain.
- c_int and K_fl from the entanglement entropy and the charge variance of the left block of l cells, fitted to (c/6 or K/2π²) ln[(2(M+1)/π) sin(πl/(M+1))] + a + (b₁ cos 2k_F l + b₂ sin 2k_F l)/chord^K, with k_F = πN/(M+1), the exponent fixed at K_gap, and 4 ≤ l ≤ M − 4. The open-chain Friedel wavevector is πN/(M+1), not π/3. Using π/3 spoils the fit badly at these sizes.
- c_FES = 6 dS/d ln ξ from VUMPS with a one-cell unit cell, over χ ≥ 24 (ξ ≥ 8 cells). The microscopic chain needs about four times the logical bond dimension for the same ξ, because port dressing splits each logical Schmidt value, and at χ = 12 and 16 it sits on a plateau (ξ ≈ 2.5) that is excluded. The logical_x chain is fitted over the same ξ window.

### Calibration on the bare XXZ chain

Table 1 runs the estimators on the second-order logical chain and compares with the Bethe ansatz at n = 1/3. K_gap converges as 1/M and K_inf(12, 24, 48) agrees with K_BA to 3 × 10⁻⁵. The M = 192 values sit slightly above that 1/M trend and were not used. c_int approaches 1 from below as M grows, and K_fl converges slowly (0.918 against 0.927 at M = 192), so both interval estimators are used only in the calibrated sense, microscopic against logical at the same M.

<!-- T1 -->
| rho | g | K_BA | K_gap M=24 | 48 | 96 | 192 | K_inf (12,24,48) | v_BA/J | v(48,96)/J | c_int M=48 | 96 | 192 | K_fl M=192 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 0.1484 | 0.92659 | 0.92516 | 0.92595 | 0.92629 | 0.92657 | 0.92662 | 0.9335 | 0.9338 | 0.984 | 0.989 | 0.993 | 0.918 |
| 1 | 0.0742 | 0.96131 | 0.96059 | 0.96099 | 0.96116 | 0.96138 | 0.96133 | 0.9006 | 0.9008 | 0.984 | 0.989 | 0.993 | 0.952 |
| 0.5 | 0.0371 | 0.98012 | 0.97976 | 0.97996 | 0.98005 | 0.98023 | 0.98013 | 0.8835 | 0.8837 | 0.983 | 0.989 | 0.994 | 0.971 |
<!-- /T1 -->

### Microscopic chain against its logical chain at equal size

Table 2 is the cleanest comparison, with the same geometry, the same filling and the same estimator. The difference dK = K_micro − K_logical_x grows as λ² (ρ = 2 at M = 24 gives 2.8 × 10⁻⁵, 1.0 × 10⁻⁴, 4.0 × 10⁻⁴, 9.0 × 10⁻⁴ at λ = 0.025, 0.05, 0.1, 0.15), and it is nearly independent of M. The M = 48 microscopic energies are converged to about 10⁻⁹ (ground states) and 10⁻⁸ (excited states), against 10⁻¹² at M ≤ 24, so at λ ≤ 0.05 the M = 48 values of dK carry an uncertainty of about 1 × 10⁻⁴ and differ from M = 24 by that much. The gap ratio E₁(micro)/E₁(logical_x), which is the velocity ratio, also departs from 1 as λ². This is what multi-cell terms of relative size λ do once the two-cell bond has been resummed. The interval estimators agree to 0.002 (c_int) and 0.007 (K_fl).

<!-- T2 -->
| lam | rho | M | K_gap micro | K_gap logical_x | dK | E1 micro / E1 logical_x | c_int micro | c_int logical_x | K_fl micro | K_fl logical_x |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.0125 | 1 | 12 | 0.96413 | 0.96414 | -7.7e-06 | 0.99993 | - | - | - | - |
| 0.0125 | 1 | 24 | 0.96503 | 0.96503 | -7.6e-06 | 0.99992 | 1.002 | 1.002 | 1.114 | 1.114 |
| 0.0125 | 0.5 | 12 | 0.98884 | 0.98885 | -8.5e-06 | 0.99988 | - | - | - | - |
| 0.0125 | 0.5 | 24 | 0.98912 | 0.98913 | -8.4e-06 | 0.99987 | 0.998 | 0.998 | 1.122 | 1.122 |
| 0.025 | 2 | 12 | 0.92753 | 0.92756 | -2.7e-05 | 0.99978 | - | - | - | - |
| 0.025 | 2 | 24 | 0.92934 | 0.92937 | -2.8e-05 | 0.99977 | 1.011 | 1.010 | 1.104 | 1.104 |
| 0.025 | 1 | 12 | 0.96852 | 0.96855 | -3.0e-05 | 0.99971 | - | - | - | - |
| 0.025 | 1 | 24 | 0.96930 | 0.96933 | -3.0e-05 | 0.99970 | 1.002 | 1.002 | 1.115 | 1.115 |
| 0.025 | 0.5 | 12 | 0.99912 | 0.99916 | -3.7e-05 | 0.99953 | - | - | - | - |
| 0.025 | 0.5 | 24 | 0.99914 | 0.99918 | -3.7e-05 | 0.99950 | 0.996 | 0.996 | 1.125 | 1.125 |
| 0.05 | 2 | 12 | 0.93072 | 0.93082 | -1.0e-04 | 0.99918 | - | - | - | - |
| 0.05 | 2 | 24 | 0.93244 | 0.93255 | -1.0e-04 | 0.99914 | 1.010 | 1.010 | 1.104 | 1.105 |
| 0.05 | 2 | 48 | 0.93321 | 0.93325 | -3.9e-05 | 0.99919 | 0.984 | 0.984 | 0.943 | 0.943 |
| 0.05 | 1 | 12 | 0.97711 | 0.97723 | -1.2e-04 | 0.99890 | - | - | - | - |
| 0.05 | 1 | 24 | 0.97768 | 0.97780 | -1.2e-04 | 0.99883 | 1.000 | 1.000 | 1.117 | 1.118 |
| 0.05 | 1 | 48 | 0.97798 | 0.97802 | -4.4e-05 | 0.99888 | 0.983 | 0.983 | 0.980 | 0.980 |
| 0.05 | 0.5 | 12 | 1.02302 | 1.02320 | -1.8e-04 | 0.99818 | - | - | - | - |
| 0.05 | 0.5 | 24 | 1.02243 | 1.02262 | -1.8e-04 | 0.99807 | 0.993 | 0.992 | 1.135 | 1.135 |
| 0.05 | 0.5 | 48 | 1.02231 | 1.02240 | -8.3e-05 | 0.99812 | 0.983 | 0.984 | 1.019 | 1.019 |
| 0.1 | 2 | 12 | 0.93508 | 0.93547 | -3.9e-04 | 0.99705 | - | - | - | - |
| 0.1 | 2 | 24 | 0.93668 | 0.93708 | -4.0e-04 | 0.99688 | 1.010 | 1.009 | 1.104 | 1.106 |
| 0.1 | 2 | 48 | 0.93743 | 0.93774 | -3.1e-04 | 0.99691 | 0.984 | 0.984 | 0.946 | 0.946 |
| 0.1 | 1 | 12 | 0.99565 | 0.99615 | -5.0e-04 | 0.99580 | - | - | - | - |
| 0.1 | 1 | 24 | 0.99574 | 0.99624 | -5.0e-04 | 0.99555 | 0.998 | 0.997 | 1.123 | 1.124 |
| 0.1 | 0.5 | 12 | 1.09686 | 1.09833 | -1.5e-03 | 0.99300 | - | - | - | - |
| 0.1 | 0.5 | 24 | 1.09430 | 1.09575 | -1.4e-03 | 0.99264 | 0.985 | 0.984 | 1.169 | 1.170 |
| 0.1 | 0.5 | 48 | 1.09362 | 1.09483 | -1.2e-03 | 0.99267 | 0.985 | 0.985 | 1.084 | 1.085 |
| 0.15 | 2 | 12 | 0.93819 | 0.93906 | -8.7e-04 | 0.99384 | - | - | - | - |
| 0.15 | 2 | 24 | 0.93969 | 0.94058 | -9.0e-04 | 0.99349 | 1.010 | 1.008 | 1.103 | 1.107 |
| 0.15 | 2 | 48 | 0.94043 | 0.94120 | -7.7e-04 | 0.99348 | 0.984 | 0.984 | 0.948 | 0.949 |
| 0.15 | 1 | 12 | 1.01962 | 1.02088 | -1.3e-03 | 0.99082 | - | - | - | - |
| 0.15 | 1 | 24 | 1.01908 | 1.02035 | -1.3e-03 | 0.99028 | 0.995 | 0.993 | 1.132 | 1.134 |
| 0.15 | 1 | 48 | 1.01908 | 1.02015 | -1.1e-03 | 0.99023 | 0.984 | 0.983 | 1.016 | 1.017 |
| 0.15 | 0.5 | 12 | 1.24879 | 1.25748 | -8.7e-03 | 0.98595 | - | - | - | - |
| 0.15 | 0.5 | 24 | 1.24146 | 1.24991 | -8.5e-03 | 0.98570 | 0.977 | 0.975 | 1.261 | 1.267 |
| 0.15 | 0.5 | 48 | 1.23957 | 1.24745 | -7.9e-03 | 0.98598 | 0.990 | 0.990 | 1.224 | 1.231 |
<!-- /T2 -->

### Thermodynamic K and v against the XXZ predictions

K micro is K_inf of the logical_x chain plus dK at M = 24 (M = 12 where that is all there is), and v micro is the logical_x velocity from its two largest sizes times the gap ratio at the same M. The error in parentheses (last digits) adds the spread between the logical_x extrapolations, the spread of dK or of the gap ratio over the other sizes (and against the direct microscopic extrapolation where it exists), and the calibration offset. Velocities are in units of the second-order hopping J₂ = 2|Jxy|λ², so the second-order prediction is the λ-independent column v_BA(g)/J₂.

<!-- T3 -->
| lam | rho | g | g_x(lam) | K micro | K_BA(g) | K_BA(g_x) | K micro - K_BA(g) | v micro / J2 | v_BA(g)/J2 | v_BA(g_x)/J2 | M max |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.0125 | 1 | 0.0742 | +0.0654 | 0.9657(1) | 0.9613 | 0.9657 | +0.0044 | 0.8656(9) | 0.9006 | 0.8651 | 24 |
| 0.0125 | 0.5 | 0.0371 | +0.0197 | 0.9893(1) | 0.9801 | 0.9893 | +0.0092 | 0.8383(8) | 0.8835 | 0.8379 | 24 |
| 0.025 | 2 | 0.1484 | +0.1391 | 0.9307(1) | 0.9266 | 0.9307 | +0.0041 | 0.8606(9) | 0.9335 | 0.8601 | 24 |
| 0.025 | 1 | 0.0742 | +0.0570 | 0.9699(1) | 0.9613 | 0.9699 | +0.0086 | 0.8321(8) | 0.9006 | 0.8318 | 24 |
| 0.025 | 0.5 | 0.0371 | +0.0015 | 0.9992(1) | 0.9801 | 0.9992 | +0.0190 | 0.7949(8) | 0.8835 | 0.7949 | 24 |
| 0.05 | 2 | 0.1484 | +0.1322 | 0.9337(3) | 0.9266 | 0.9338 | +0.0072 | 0.7960(19) | 0.9335 | 0.7965 | 48 |
| 0.05 | 1 | 0.0742 | +0.0408 | 0.9781(3) | 0.9613 | 0.9782 | +0.0168 | 0.7701(16) | 0.9006 | 0.7709 | 48 |
| 0.05 | 0.5 | 0.0371 | -0.0391 | 1.0220(4) | 0.9801 | 1.0222 | +0.0419 | 0.7133(13) | 0.8835 | 0.7146 | 48 |
| 0.1 | 2 | 0.1484 | +0.1225 | 0.9379(4) | 0.9266 | 0.9383 | +0.0113 | 0.6885(17) | 0.9335 | 0.6905 | 48 |
| 0.1 | 1 | 0.0742 | +0.0067 | 0.9958(1) | 0.9613 | 0.9963 | +0.0345 | 0.6625(11) | 0.9006 | 0.6653 | 24 |
| 0.1 | 0.5 | 0.0371 | -0.1510 | 1.0927(9) | 0.9801 | 1.0942 | +0.1126 | 0.5620(8) | 0.8835 | 0.5662 | 48 |
| 0.15 | 2 | 0.1484 | +0.1151 | 0.9408(5) | 0.9266 | 0.9417 | +0.0142 | 0.6005(15) | 0.9335 | 0.6043 | 48 |
| 0.15 | 1 | 0.0742 | -0.0353 | 1.0187(8) | 0.9613 | 1.0200 | +0.0574 | 0.5682(12) | 0.9006 | 0.5738 | 48 |
| 0.15 | 0.5 | 0.0371 | -0.3315 | 1.2374(22) | 0.9801 | 1.2459 | +0.2573 | 0.4134(8) | 0.8835 | 0.4195 | 48 |
<!-- /T3 -->

Two columns carry the result. K micro − K_BA(g) shrinks linearly as λ decreases at each ρ, and its slope is the one predicted by the two-cell shift of the anisotropy. With dK/dΔ = −0.444, −0.493, −0.521 at n = 1/3 and dg_x/dλ = −0.41, −0.70, −1.39, the predicted slopes are 0.18, 0.35, 0.72, and the λ → 0 values measured are about 0.16 (ρ = 2, λ = 0.025), 0.35 (ρ = 1, λ = 0.0125) and 0.74 (ρ = 1/2, λ = 0.0125). The column K_BA(g_x) then matches K micro to the λ² residual of Table 2. The velocity has a larger first-order correction, about −3λ relative, coming from the renormalized hopping J_x, and v micro equals v_BA(g_x)·J_x to 0.05 to 1.5%.

Where M = 12, 24 and 48 all exist for the microscopic chain, its own extrapolation (Table 4) agrees with the combined estimate of Table 3 within the quoted errors.

<!-- T4 -->
| lam | rho | K_inf micro (12,24,48) | K_inf logical_x (12,24,48) | v micro (24,48)/J2 | v logical_x (24,48)/J2 |
|---|---|---|---|---|---|
| 0.05 | 2 | 0.93392 | 0.93385 | 0.7966 | 0.7972 |
| 0.05 | 1 | 0.97829 | 0.97820 | 0.7706 | 0.7714 |
| 0.05 | 0.5 | 1.02231 | 1.02222 | 0.7136 | 0.7149 |
| 0.1 | 2 | 0.93814 | 0.93829 | 0.6890 | 0.6911 |
| 0.1 | 0.5 | 1.09332 | 1.09416 | 0.5621 | 0.5662 |
| 0.15 | 2 | 0.94117 | 0.94172 | 0.6009 | 0.6048 |
| 0.15 | 1 | 1.01926 | 1.01999 | 0.5684 | 0.5740 |
| 0.15 | 0.5 | 1.23888 | 1.24586 | 0.4134 | 0.4192 |
<!-- /T4 -->

### Finite-entanglement scaling

At equal ξ the microscopic entropy exceeds the logical one by a constant (the port dressing across the inter-cell cut, growing with λ), and the slopes, which are c/6, agree within the fit errors. On this window (ξ ≤ 35 cells) the calibrated logical chain gives 0.94 to 0.98, so the estimator carries a systematic of a few percent at these ξ, consistent with `../routeA_demonstrator/` (0.90 to 0.92 for ξ ≤ 32 at half filling). The VUMPS gradient errors stop at 10⁻⁶ to 10⁻⁵ for the microscopic chain (10⁻⁹ for the logical one), the same behaviour as in the half-filled runs there.

<!-- T5 -->
| lam | rho | chi (micro) | xi micro (cells) | c micro | c logical_x, same xi window | S micro - S logical_x at equal xi | n reached |
|---|---|---|---|---|---|---|---|
| 0.05 | 2 | 24 to 64 | 8.4 to 29.8 | 0.960 +- 0.052 | 0.936 +- 0.035 | 0.021 to 0.028 | 0.3331 to 0.3333 |
| 0.05 | 1 | 24 to 64 | 8.6 to 33.8 | 1.001 +- 0.028 | 0.936 +- 0.035 | 0.028 to 0.036 | 0.3331 to 0.3334 |
| 0.05 | 0.5 | 24 to 48 | 8.7 to 21.6 | 0.967 +- 0.016 | 0.977 +- 0.059 | 0.029 to 0.044 | 0.3334 to 0.3335 |
| 0.1 | 2 | 24 to 64 | 8.4 to 31.6 | 0.986 +- 0.024 | 0.936 +- 0.035 | 0.069 to 0.078 | 0.3332 to 0.3334 |
| 0.1 | 1 | 24 to 64 | 8.6 to 32.7 | 0.986 +- 0.023 | 0.936 +- 0.036 | 0.078 to 0.086 | 0.3336 to 0.3338 |
| 0.1 | 0.5 | 24 to 64 | 8.8 to 33.9 | 0.993 +- 0.021 | 0.982 +- 0.036 | 0.111 to 0.120 | 0.3346 to 0.3348 |
| 0.15 | 2 | 24 to 64 | 8.4 to 31.6 | 0.970 +- 0.044 | 0.936 +- 0.035 | 0.127 to 0.132 | 0.3334 to 0.3335 |
| 0.15 | 1 | 24 to 64 | 8.6 to 33.4 | 0.975 +- 0.049 | 0.982 +- 0.036 | 0.146 to 0.149 | 0.3344 to 0.3345 |
| 0.15 | 0.5 | 24 to 64 | 8.9 to 34.6 | 0.989 +- 0.016 | 0.982 +- 0.036 | 0.203 to 0.211 | 0.3382 to 0.3383 |
<!-- /T5 -->

### What this says about the RG input

The analytic gap identified in C47 and in the collaborator's SINGLE_CHANNEL_RG.md is the infrared flow of the matched retarded action, which is the XXZ reference plus a summable O(λ³) correction. The numerics show what that flow does in this family. The corrections act as local renormalizations of the two marginal data (anisotropy and hopping), at first order in λ through the two-cell bond and at second order through longer-range terms, and the chain then behaves as the XXZ Luttinger liquid at those renormalized values, with c = 1. There is no sign of a relevant perturbation, of a gap, or of a drift of K with system size beyond the XXZ chain's own 1/M approach. This is heuristic as a statement about the RG [heuristic], and numerical as a statement about K, v and c at the listed sizes [numerical]. It also says that a theorem stated "near XX" covers the tuned chain only in the order ρ first, λ last, with λ well below g/|dg_x/dλ|, which is 0.36, 0.11 and 0.027 for ρ = 2, 1, 1/2.

## Reproducing

numpy and scipy only. `chain.py` (micro and logical open chains, MPOs, the encoding of a logical MPS into the micro chain), `dmrg.py` (two-site DMRG with site-dependent dimensions, excited states by projection, block charge moments), `bethe.py` (K and v of XXZ at any filling), `run_obc.py` (one open-chain point, writes `data/obc/NAME.json`), `vumps_fes.py` (one VUMPS FES point, uses `../d2_followup_220_1282/vumps.py`, writes `data/fes/NAME.json`), `analyze.py` (Tables 1 to 5, `data/summary.json`, `python3 analyze.py --md`, or `--readme` to refill the tables here), `test_small.py` (ED checks). Task lists: `tasks_typhon.txt` (Typhon array `array.slurm`, M = 24 and 12 and VUMPS; its M = 48 lines were cancelled and rerun from `tasks_typhon3.txt` with χ ≤ 160), `tasks_typhon2.txt` (λ ≤ 0.025; its M = 48 lines were rerun from `tasks_typhon4.txt` with χ ≤ 128), and the matching `array*.slurm`, `tasks_local*.txt` (logical chains, run locally with 2 threads each). Typhon arrays 17522077 (VUMPS), 17522172 and 17522173 (M = 24, 12), 17522428 (M = 48), 17522224 and 17523039 (λ ≤ 0.025). The microscopic runs take 20 minutes (M = 12) to 8 hours (M = 48) on 4 cores. Logs are in `data/logs/` (local runs) and in the Typhon scratch folder `demonstrator_luttinger` (cluster runs).
