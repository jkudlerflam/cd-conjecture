# Encoded route-B chain: does the logical chain flow to c = r/2?

## Claim

Class: the chain of route-B cells of `contrib/programmable_cell/` §3(B), each cell one site of dimension D = q⁴2ⁿ (pL, sL, n logical qubits, sR, pR), H = Σ_j H_cell(j) + λ Σ_j Σ_a F_a(pR_j) F_a(pL_{j+1}), H_cell = H_L + H_R as there (η_a = iγ_aP on the right), F_a the r = 2q₊q₋ odd port channels, Δ = 1. Translation invariant, nearest neighbour in big sites, Hermitian. At first order it is r Majorana chains, c = r/2. Question: the Virasoro c at finite λ/Δ.

1. q = 2 (q₊ = q₋ = 1, r = 2, n = 1, D = 32). c = 0.98 to 1.05, each within its error of 1, for λ/Δ = 0.05 to 0.8 (g/Δ = 0.5) and 0.05 to 0.2 (g/Δ = 0.2), with z = 1 and no translation breaking. An exact U(1) with half-odd charge per cell, particle-hole symmetry, and an attractive second-order XXZ anisotropy make umklapp irrelevant. Verdict: c = r/2 confirmed. [numerical]
2. q = 3 (q₊ = 1, q₋ = 2, r = 4, n = 2, D = 324). One graded port per end gives r = 4 with exact compressions, so no second port is needed. The logical chain with the bond to all orders in λ gives c = 2.1 ± 0.3 (λ = 0.4, ξ ≤ 11 cells), the second-order chain c = 1.97 ± 0.06 (ξ ≤ 23). The full D = 324 chain is resolved only to ξ ≈ 2.5 cells, where its logical charge matches them. Mechanism [one-loop]: the SU(2) current coupling is marginally irrelevant, and a chemical potential ∝ λ²/Δ from q₊ ≠ q₋ dopes the charge sector before umklapp can open a Mott gap. With the sign of the λ² terms reversed, c → 1, so the estimator does see a relevant coupling when there is one. Verdict: c = 2 at all accessible scales for λ/Δ ≤ 0.4, one-loop beyond. Near λ ≈ 0.75Δ (g = 0.5) two cells have a level crossing and the encoded description ends. [numerical, one-loop]

## Status

Numerical for q = 2. For q = 3, numerical on logical chains derived from the full one up to ξ ≈ 25, and one-loop beyond. The full q = 3 chain at χ of several hundred would upgrade it.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), sub-session of the coordinating session, 2026-10-27. No Codex/GPT was used.
- Code: new, on top of `../programmable_cell/majorana_bank_check.py` (cell) and `../d2_followup_220_1282/vumps.py` (VUMPS).
- Checks: exact cell degeneracy and compressions (10⁻¹⁵). Sparse, dense and factorized MPOs agree to 12 digits. The second-order chain matches two-cell diagonalization with residual ∝ λ³. The continuum projection respects SU(2) × U(1) to 10⁻¹² and is calibrated on the Hubbard chain. All estimators were run on Ising and the first-order chains.
- Not rerun by anyone else.

## Method

VUMPS with the cell as one site (a low-rank null-space expansion avoids two-site tensors at D = 324), in unit cells of 4 cells (q = 2, whose finite-χ state is a period-4 spiral) or 2 (q = 3), so translation breaking shows in cell-resolved entropies, bond energies and U(1) charge. Finite-entanglement scaling: c = 6 dS/d ln ξ over the last three χ (cell-boundary S, ξ in cells), and e(χ) − e∞ ∝ ξ^−p, p = 1 + z. Open chains: two-site DMRG on [pL | sL⊗logical⊗sR | pR] (inner bonds qχ, lossless), Calabrese–Cardy fits with a staggered term.

Logical chains (q = 3). `effective.py` gives the Schrieffer–Wolff chain through O(λ²), a two-cell term T2 and a three-cell term T3 (middle cell excited). `core_bond` replaces the two-cell part by its all-order value from two-cell diagonalization. `marginal.py` expands quartic terms with c_{a,j} → ψ_{aR} + (−1)^j ψ_{aL} and projects on so(4) = su(2)₊ ⊕ su(2)₋ currents (su(2)₊ holds the U(1) charge).

## Results

Calibration: Ising 0.498 (ξ to 2190 cells), first-order r = 2 chain 0.988 (ξ to 491), first-order r = 4 chain 2.08 (χ = 128 to 256, ξ = 26 to 44, pairwise slopes 1.92 to 2.08). The penalty-projected excited state reproduces the exact open-chain gap of the r = 2 chain to six digits.

q = 2 (c ± spread of the last pairwise slopes; leak = weight of the left pair outside |Ω⟩):

| g | λ | λ/α² | χ | ξ | c | p | leak |
|---|---|---|---|---|---|---|---|
| 0.2 | 0.05 | 0.46 | 32–64 | 18–52 | 1.054 ± 0.107 | 1.42 | 0.07 |
| 0.2 | 0.1 | 0.92 | 32–64 | 18–47 | 0.995 ± 0.050 | 1.55 | 0.07 |
| 0.2 | 0.2 | 1.84 | 32–64 | 17–44 | 0.985 ± 0.011 | 1.71 | 0.10 |
| 0.5 | 0.05 | 0.19 | 24–48 | 18–43 | 0.977 ± 0.088 | 1.79 | 0.18 |
| 0.5 | 0.1 | 0.38 | 32–64 | 25–59 | 0.993 ± 0.075 | 2.05 | 0.19 |
| 0.5 | 0.2 | 0.77 | 32–64 | 23–57 | 0.990 ± 0.070 | 2.13 | 0.21 |
| 0.5 | 0.4 | 1.54 | 48–96 | 37–89 | 1.006 ± 0.030 | 1.99 | 0.26 |
| 0.5 | 0.8 | 3.07 | 32–64 | 19–47 | 1.007 ± 0.047 | 1.93 | 0.37 |

All runs have Q = 0 and cell-to-cell differences of S and bond energy below 10⁻⁸. p ≈ 2 at g = 0.5 gives z = 1 (at g = 0.2 the energy differences sit near the convergence threshold). Open chains (g = 0.5, L = 16 at χ = 64, L = 32 at χ = 32) give c = 0.85 to 1.01, against 0.93 for the first-order chain at L = 32, so they agree within that finite-size bias. At second order T2 adds −0.120 λ² Z_jZ_{j+1} (g = 0.5; −0.026 λ² at g = 0.2) and T3 a correlated hopping: XXZ with Δ_xxz = −0.46λ/Δ (−0.24λ/Δ), K > 1, umklapp (dimension 4K) irrelevant. Half-odd charge per cell excludes a trivial gap, and dimerization would need K < 1/2. The two-cell level crossing is near λ ≈ 1.0Δ.

q = 3 couplings per λ², normalized so that the Hubbard chain with U = 1 has +0.5 (charge) and −0.5 (spin), positive meaning marginally relevant ("all orders": exact two-cell bond plus second-order T3):

| g | λ | spin | charge, U(1) axis | charge, umklapp | μ/λ² |
|---|---|---|---|---|---|
| 0.5 | O(λ²) | −0.495 | −0.048 | −0.796 | 0.0423 |
| 0.5 | 0.2, all orders | −0.488 | −0.052 | −0.783 | 0.0386 |
| 0.5 | 0.4, all orders | −0.461 | −0.059 | −0.732 | 0.0293 |
| 0.2 | O(λ²) | −0.170 | −0.044 | −0.253 | 0.0297 |
| 0.2 | 0.4, all orders | −0.150 | −0.016 | −0.237 | 0.0176 |

The spin coupling is irrelevant throughout. With |g_⊥| ≫ |g_∥| the charge sector alone would open a Mott gap ∼ v exp(−π/2y), y ≈ 0.7λ/Δ, about 10⁻³ at λ = 0.4. But T2 contains μλ² i(c₀c₁ + c₂c₃) on each cell of each bond (absent for q₊ = q₋), a fermion chemical potential μ_f = 4μλ² (0.027 at λ = 0.4) that dopes the charge sector by δn = 2μ_f/πv = 0.061 per cell, the value VUMPS finds (−0.062). Keeping only the 68 cell states below 2Δ changes these couplings by under 2% (6% for the U(1)-axis one).

VUMPS on the logical chains (g = 0.5). Second order, λ = 0.4: c = 1.97 ± 0.06 (χ = 32 to 64, ξ = 13 to 23). All-order bond, λ = 0.4: c = 2.1 ± 0.3 (χ = 24 to 48, ξ = 6 to 11, pairwise 1.76 to 2.39). Doping sets in abruptly once v/ξ < μ_f, as one χ step with jumps in Q and S, and pairs straddling it are dropped. With the λ² terms multiplied by −4 the spin coupling turns relevant and c = 1.01 (ξ = 7 to 15, Q = +0.144).

Full chain, D = 324, χ ≤ 32. The cell boundary carries the bond's port fluctuations besides the logical entanglement, so ξ reaches only 1.4 to 2.5 cells, too short for a slope. No translation breaking at λ = 0.2, 0.8 or g = 0.2 (below 10⁻⁹). At g = 0.5, λ = 0.4 the run keeps a period-2 asymmetry (ΔS = 0.04) at χ = 24 that the 68-state truncated cell does not show. The logical charge, −0.014 (λ = 0.2) and −0.022 (λ = 0.4), matches the logical chains at comparable ξ (−0.016, −0.033).

## Translation to the fixed-d₀ construction

In `programmable_cell` §4 a cell is a segment of b sites with d₀ = 4k_c q. The big site is its logical content (port and partner on the register factors ℂ^q of the end sites and their neighbours, logical qubits on ancillas), so D = q⁴2ⁿ costs no extra local dimension. Δ and g are energy terms of the port-idle history cell, λ the bond term of its item 5. The conditions found are dimensionless: λ below the two-cell level crossing (0.75Δ to 1.0Δ here), while λ/(α²Δ) may reach 3. Microscopically Δ is the clock-tower gap ∝ J_c/b², so λ and the bandwidth λα² shrink with b, which the infrared analysis does not see. Two lessons carry over. q₊ ≠ q₋, which route B needs for c > d₀ − 1, generically produces an O(λ²) chemical potential that removes umklapp. The SU(q₊) × SU(q₋) current couplings are marginal and their signs must be computed for the target (q₊, q₋). For q = 3 only SU(2) appears and its sign is benign. For q₊ = 8k_c they are unknown.

## Reproducing

numpy and scipy only. `cells.py`, `effective.py`, `marginal.py`, `core_couplings.py` run in minutes. VUMPS: `run_fes.py {sparse|trunc|effective|core|logical|ising} --q --qp --n --g --lam --ncells --chis ... --out` (`--lanczos --eig_factor 0.1` for D = 324). Open chains: `fdmrg_md.py encoded --q --g --lam --L --chis`. Tables: `analyze.py fes|obc`. Task lists and slurm scripts are beside them. Typhon jobs 17512493 (controls, q = 2), 17512848 (q = 3, `*_arr.json`), 17512851 (q = 2, g = 0.2), 17512979 and 17512980 (open chains). Logical-chain, truncated-cell and some open-chain runs were `srun --overlap` steps inside 17512493 because of queue limits. Data: `data/typhon/` (cluster), `data/core/` and `data/eff/` (local logical-chain runs), `data/dev/` (unit-cell and factorization tests).

## Relation to the rest of the repository

This tests §5 item 2 of `contrib/programmable_cell/` (REPORT §9.6) and complements `contrib/ir_control/`, whose detuning argument holds here for umklapp but not for the SU(2) coupling between degenerate Fermi points, whose sign decides. "Where things stand" could say that the encoded chain flows to c = r/2 for q = 2, and for q = 3 at all accessible scales, with the non-abelian coupling signs the open issue at large q.
