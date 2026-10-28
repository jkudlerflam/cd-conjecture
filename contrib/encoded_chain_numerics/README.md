# Encoded route-B chain: does the logical chain really flow to c = r/2?

## Claim

Class: the chain of route-B cells of `contrib/programmable_cell/` §3(B), each cell treated as one site of dimension D = q⁴2ⁿ (port pL, partner sL, n logical qubits, sR, pR). H = Σ_j H_cell(j) + λ Σ_j Σ_a F_a(pR_j) F_a(pL_{j+1}), with H_cell = H_L + H_R as in that README (η_a = iγ_aP at the right end), F_a the r = 2q₊q₋ odd port channels, Δ = 1. This is translation invariant, nearest neighbour in the big sites and Hermitian. At first order it is r uniform Majorana chains, c = r/2. The question is the Virasoro c of the full chain at finite λ/Δ.

1. q = 2 (q₊ = q₋ = 1, r = 2, n = 1, D = 32). c = 1.00 ± 0.03 for λ/Δ = 0.05, 0.1, 0.2, 0.4, 0.8 at g/Δ = 0.5 and for λ/Δ = 0.05 to 0.4 at g/Δ = 0.2, with z = 1 and no translation breaking. The chain has an exact U(1) with half-odd charge per cell and is particle-hole symmetric, and at second order it is an XXZ chain on the attractive side (anisotropy −0.46 λ/Δ at g = 0.5), so umklapp is irrelevant. Verdict: c = r/2 confirmed. [numerical]
2. q = 3 (q₊ = 1, q₋ = 2, r = 4, n = 2, D = 324). One graded port per end already gives r = 4 with exact compressions, so no second port is needed. RESULT-Q3. Verdict: VERDICT-Q3.
3. Mechanism for q = 3 [one-loop, heuristic]. The chain has an exact U(2) = (U(1) × SU(2))/Z₂. At second order the spin current coupling has the sign of the repulsive Hubbard chain (marginally irrelevant). At commensurate filling the charge sector would flow to a Mott gap, but q₊ ≠ q₋ generates a chemical potential μ ∝ λ²/Δ that dopes it first, leaving a Luttinger liquid. Reversing the sign of the λ² terms gaps the spin sector (c → 1 numerically), so the estimators do see a marginally relevant coupling when one is present.

## Status

Numerical for q = 2. For q = 3, numerical at the accessible scales (ξ ≲ XI3 cells) plus a one-loop argument for longer scales. A spectral check on rings, or the full q = 3 chain at larger χ, would upgrade it.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), sub-session of the coordinating session, 2026-10-27. No Codex/GPT was used.
- Code: new here (`cells.py`, `umps_md.py`, `sparse_mpo.py`, `fdmrg_md.py`, `effective.py`, `marginal.py`, `run_fes.py`, `analyze.py`), built on `../programmable_cell/majorana_bank_check.py` (cell) and `../d2_followup_220_1282/vumps.py` (VUMPS core).
- Checks: the cell reproduces α and exact 2ⁿ-fold degeneracy (residual 10⁻¹⁵). The sparse big-site MPO, the dense big-site MPO and the three-sub-site factorization give identical energies (12 digits). The second-order effective Hamiltonian matches two-cell exact diagonalization with residual ∝ λ³. The one-loop projection respects SU(2) × U(1) to machine precision and is calibrated on the Hubbard chain. All estimators were run on Ising and on the first-order chains (exact c = 1/2, 1, 2).
- Not rerun by anyone else.

## Method

Each cell is one site. For VUMPS (Zauner-Stauber et al. 2018, multi-site, with a low-rank null-space expansion so the D = 324 site never needs two-site tensors) the MPO is sparse. Unit cells: 4 cells for q = 2 (the finite-χ state of the first-order chain is a period-4 spiral; 1- and 2-cell runs fail or converge badly) and 2 cells for q = 3. Two cells per unit cell detect one-cell translation breaking through the cell-resolved bond energies, entropies and charges.

Estimators. (i) Finite-entanglement scaling: S at the cell boundary against ln ξ (ξ from the transfer matrix, in cells), c = 6 dS/d ln ξ over the last three χ. (ii) Open chains of L cells by two-site DMRG on the factorization [pL | sL ⊗ logical ⊗ sR | pR], whose inner bonds (q times the cell-boundary χ) make it a lossless big-site MPS: Calabrese–Cardy fits of S(ℓ) at cell boundaries with a staggered term, and the lowest gap from a penalty-projected second state (z from L·gap). (iii) In VUMPS, e(χ) − e∞ ∝ ξ^−p with p = 1 + z.

The second-order Schrieffer–Wolff chain (`effective.py`) is exact through O(λ²) and contains a two-cell term T2 and a three-cell term T3 (middle cell excited). `marginal.py` expands its quartic terms in continuum fields c_{a,j} → ψ_{aR} + (−1)^j ψ_{aL} and projects onto so(4) = su(2)₊ ⊕ su(2)₋ currents (su(2)₊ contains the U(1) charge).

## Results

### Calibration

| model (exact c) | χ | ξ (cells) | c (FES) |
|---|---|---|---|
| Ising (1/2) | 8–64 | 34–2190 | 0.498 |
| first-order chain r = 2 (1) | 8–128 | 11–491 | 0.988 |
| first-order chain r = 4 (2) | 8–256 | 2–44 | CAL4 |

The first-order r = 2 chain reproduces the exact open-chain gap 4t sin(π/2(L+1)) to six digits, which checks the excited-state DMRG.

### q = 2

TABLE-Q2

Every run has Q = 0 per cell, cell-to-cell differences in S and bond energy below 10⁻⁸, and p ≈ 1.8–2.1 (z = 1). The leakage out of the Bell pair grows from 0.07 (g = 0.2, λ = 0.05) to 0.37 (g = 0.5, λ = 0.8), so the strongest points are far from perturbative, yet c does not move. OBC-Q2

Second order (`effective.py`): T2 contains only the constant and −0.120 λ² Z_jZ_{j+1} (g = 0.5; −0.026 λ² at g = 0.2), and T3 a correlated hopping. With first-order hopping λα² this is XXZ with Δ_xxz = −0.46 λ/Δ (g = 0.5) or −0.24 λ/Δ (g = 0.2), so K > 1 and umklapp (dimension 4K) is irrelevant. Charge per cell is half-odd, so a trivially gapped phase is excluded and dimerization would need K < 1/2, the opposite side.

### q = 3

TABLE-Q3

Second-order couplings per λ² (g = 0.5; g = 0.2 in brackets), in a normalization where the Hubbard chain with U = 1 has +0.5 in the charge algebra and −0.5 in the spin algebra (positive means marginally relevant):

- spin su(2)₋, isotropic: −0.495 (−0.170), marginally irrelevant;
- charge su(2)₊: g_∥ = −0.048 (−0.044) along the U(1) axis, g_⊥ = −0.796 (−0.253) for the charge-2 (umklapp) components. Since |g_⊥| ≫ |g_∥| this flows to a Mott gap at commensurate filling;
- chemical potential: T2 contains 0.0423 λ² (0.0297 λ²) i(c₀c₁ + c₂c₃) per cell per bond, absent for q₊ = q₋. Its doping, δn = 2μ/(πv) with v = 4λα², is 0.061 per cell at λ = 0.4, which is the charge per cell found by VUMPS on the effective chain (−0.0617). Against a Mott gap of order v exp(−π/2g̃_⊥), μ wins for every λ in the range, so umklapp is detuned.

The same projection on the truncated cell of 68 states (everything below 2Δ) changes these numbers by less than 2%. The effective chain gives c = 1.97 at λ = 0.4 and 0.8 (ξ up to 23), and EFFLOW. With the λ² terms multiplied by −4 the spin coupling becomes relevant, and c falls to 1.01 at ξ ≈ 15 while the charge sector stays doped (Q = 0.144). With them multiplied by +4 c stays near 2.

## Translation to the fixed-d₀ construction

In `programmable_cell` §4 the cell is a segment of b sites of the C7 selector chain, with d₀ = 4k_c q. The big site here is its logical content: the port and partner live on the register factor ℂ^q of the two end sites and the next ones inward, and the n logical qubits on ancilla registers, so D = q⁴2ⁿ costs no extra local dimension. Δ and g are energy terms of the port-idle history cell, reproduced up to O(1/J_c), and λ is the bond term of item 5. The numerics fix the window that the construction must hit. λ must sit below the cell gap (here up to λ = 0.8Δ, beyond first-order perturbation theory), and the relevant ratio of second to first order is λ/(α²Δ), between 0.2 and 3 here. In the microscopic chain the cell gap is the clock-tower gap ∝ J_c/b², so λ, and with it the logical bandwidth ∝ λα², scale down with b, but nothing in the infrared analysis depends on the absolute scale. Two lessons carry over. First, q₊ ≠ q₋ (needed for c > d₀ − 1 in route B) generically produces a chemical potential at O(λ²), which removes umklapp, a useful effect. Second, the non-abelian current couplings of SU(q₊) × SU(q₋) at first order are marginal, and their signs must be computed for the target (q₊, q₋). For q = 3 only SU(2) is present and its sign is benign. For q₊ = 8k_c these signs are unknown.

## Reproducing

numpy and scipy only. `python3 cells.py` (cell checks, seconds), `python3 effective.py` (second order, a minute), `python3 marginal.py` (couplings, a minute). VUMPS: `run_fes.py {sparse|trunc|effective|logical|ising} --q --qp --n --g --lam --ncells --chis ... --out`, with `--lanczos --eig_factor 0.1` for q = 3. Open chains: `fdmrg_md.py encoded --q --g --lam --L --chis --nstates 2`. `analyze.py fes|obc data/*.json` makes the tables. Task lists: `fes_tasks.txt`, `fes_q3_tasks.txt`, `obc_tasks.txt`, `eff_tasks.txt`, with the slurm scripts beside them.

Typhon jobs (partition typhon): 17512493 (controls and q = 2 VUMPS), 17512848 (q = 3 VUMPS, `*_arr.json`), 17512851 (q = 2, g = 0.2), 17512979 and 17512980 (open chains). Several runs (effective chains, truncated cells, some open chains) ran as `srun --overlap` steps inside the 17512493 allocations because of queue limits. Times: q = 2 VUMPS to χ = 192 in hours on 8 cores; q = 3 at D = 324 takes about an hour per χ step at χ ≈ 32 and grows like χ³. Data in `data/`.

## Relation to the rest of the repository

This tests §5 item 2 of `contrib/programmable_cell/` and the infrared step of REPORT §9.6 numerically, and complements the analytic `contrib/ir_control/`. Its §3 claim that the SO(2n) point sits on a separatrix and that incommensurate detuning removes the non-Gaussian channels is borne out for q = 3 in one respect (the chemical potential removes umklapp) and qualified in another (the SU(2) current coupling within degenerate Fermi points survives detuning, and its sign decides). "Where things stand" could say that the encoded logical chain flows to c = r/2 for q = 2 and, at accessible scales, for q = 3.
