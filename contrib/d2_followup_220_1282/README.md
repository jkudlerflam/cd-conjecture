# Points 220 and 1282 have c = 1

## Claim

The chains are H = Σ_j [σ_j^T J σ_{j+1} + f·σ_j] on qubits, nearest neighbour, one-site translation invariant and Hermitian, with (J, f) of points 220 and 1282 in `candidates_params.json` (copied from `../d2_classification`, bond convention of `ed.py`). The quantity is the Virasoro c.

Both points have c = 1. Each is an incommensurate "floating" spiral, a single compact boson, with dynamical exponent consistent with z = 1.

| point | c, finite-entanglement scaling | spiral q (×2π) | spiral dimension x_q | ring c_S, L = 21–26 |
|---|---|---|---|---|
| 220 | 1.00 ± 0.01 (ξ ≤ 730) | 0.186 | 0.16 | 1.08–1.12 |
| 1282 | 1.01 ± 0.02 (ξ ≤ 155) | 0.300 | 0.34 | 1.18–1.23 |

Neither point is a counterexample to c ≤ d − 1. The doubling values 6ΔS/ln 2 = 1.93 and 2.16 are open-boundary artifacts (section 4).

## Status

Numerical. A quasiparticle-ansatz velocity would add a third estimate, the Casimir energy of the rings.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), a subagent of the coordinating session, 2026-10-27. No Codex/GPT was used.
- New code: `vumps.py` (multi-site VUMPS, any finite-range MPO), `fdmrg.py` (open-chain DMRG) and the analysis scripts. Ring ED uses `../d2_classification/ed.py` and `run_point.py` unchanged.
- Checks. VUMPS gives −4/π for Ising to 6·10⁻¹¹, is unchanged by a complex rotation of the Ising chain, and matches the free-fermion control below to 2·10⁻⁴ at χ = 64. `fdmrg.py --edgefield` reproduces the L = 96 energy of `d2_classification/dmrg.py` to 10⁻⁹. Every estimator was run on controls with known c. Nobody else has rerun anything.

## Method and results

### 1. Finite-entanglement scaling (boundary-free)

For an infinite MPS at bond dimension χ a CFT gives S = (c/6) ln ξ + const, with ξ = −n/ln|λ₁/λ₀| from the n-site cell transfer matrix. I run VUMPS (Zauner-Stauber et al., PRB 97, 045145) for χ = 8 to 128, and c_fit is 6 × the slope of the cell-averaged S against ln ξ over the last four χ.

A one-site cell cannot hold the spiral order that the finite-χ state wants, and VUMPS there does not converge (as for XX). A scan of cells n = 1 to 37 (`data/cellscan/`) favors n = 11 and 16 for 220 and n = 10 for 1282, commensurate with the spiral. The soft sliding mode makes plain VUMPS stall, so I take damped steps A_C ← e^(−τH_AC)A_C, C ← e^(−τH_C)C, with τ chosen so the variational energy never rises. A run stops when the gradient is below 10⁻⁸, or below 10⁻³ with S stationary to 10⁻⁷ over 20 steps.

| model (exact c) | cell | χ in fit | ξ in fit | c_fit | κ |
|---|---|---|---|---|---|
| Ising (1/2) | 1 | 48–128 | 1240–7810 | 0.498 | 1.86 |
| XX (1) | 2 | 48–128 | 140–490 | 1.000 | 1.28 |
| XXZ, Δ = 0.5 (1) | 2 | 64–128 | 170–430 | 1.008 | 1.27 |
| XXZ, Δ = −0.5 (1) | 2 | 48–128 | 180–630 | 1.014 | 1.28 |
| XY + DM + field, free fermions (1) | 5 | 24–64 | 78–290 | 0.992 | 1.34 |
| cluster-Ising, three-site (3/2) | 2 | 48–128 | 41–109 | 1.500 | 1.02 |
| two decoupled XX chains, next-nearest neighbour (2) | 4 | 48–128 | 25–60 | 1.94 | 0.89 |
| **220** | 11 | 48–128 | 208–730 | 1.002 | 1.28 |
| 220 | 16 | 24–64 | 48–164 | 1.03 | 1.25 |
| **1282** | 10 | 32–96 | 38–155 | 1.008 | 1.28 |
| 1282 | 23 | 24–64 | 25–90 | 1.014 | 1.30 |

The free-fermion control (J = diag(1, 0.6, 0), D_z = 0.5, f = (0, 0, 0.6)) is the closest analogue, with a large DM term, no conserved magnetization and Fermi points at 0.44π and 0.77π (`freefermion.py`).

For the candidates, the slopes between consecutive χ are as follows.
- 220 (n = 11), χ = 24 to 128: 1.027, 1.019, 1.005, 1.005, 1.001, 1.002.
- 1282 (n = 10), χ = 24 to 96: 0.98, 0.99, 1.04, 1.00, 0.99. A continuation of the n = 23 run to χ = 96 (ξ = 153, gradient 3·10⁻⁴) adds 0.99.

Single bonds of the cell give 0.94 to 1.10. The exponent κ of ξ ∝ χ^κ also matches the c = 1 controls (theory 1.34, 1.05, 0.87 for c = 1, 3/2, 2). A c = 3/2 or 2 chain would show slopes near 1.5 or 2 already at ξ ~ 50, as the controls do. One run (1282, n = 13) became non-injective at χ = 64 and is discarded.

### 2. Correlations in the infinite state

The cell magnetization splits into a uniform part m₀ and a spiral amplitude |M_q| ∝ ξ^(−x_q). The same fit gives 0.250 for XX and 0.98 for the free-fermion density wave (exact 1/4 and 1).

- 220 has m₀ = 0.33. In the n = 11 cell the spiral "unlocks" above χ = 16 (every bond then carries the same S), and the leading transfer-matrix pair has phase 11q ≡ ±0.041, so q = 0.1855. The n = 16 cell (q locked at 3/16) gives x_q = 0.16.
- 1282 has m₀ = 0.17, and q = 3/10 beats 7/23, 8/27 and 5/17 in energy. x_q = 0.34 in both the n = 10 and n = 23 cells.

So spin correlations decay as r^(−2x_q) cos(2πqr), quasi-long-range spiral order with one sliding mode, and K = 1/(4x_q) ≈ 1.5 (220) and 0.75 (1282) in a normalization where XX has K = 1.

### 3. Ring ED with momentum resolution

Six levels in every momentum sector, to L = 26 for the candidates (27 and 28 still running) and to L = 24 for two incommensurate c = 1 controls, the free-fermion chain and the twisted XXZ chain `xxzm07_dm_h`.

- Towers. The lowest excitation is a "rotor" pair at relative momentum ±N₀ (units 2π/L), with N₀/L = 0.16 to 0.21 (220) and 0.29 to 0.32 (1282) for most L, so N₀ ≈ qL. States at ±2N₀ mod L follow (cleanly for 220, less so for 1282). These are charge-one and charge-two states of the sliding mode, with L·ΔE/2π = 0.003 to 0.4 depending on qL mod 1. Descendants cannot be isolated in this dense spectrum, so I quote no velocity or Casimir c.
- z. Over L = 12 to 26 the median of L·ε_n, n = 5 to 20 (all momenta), scales as L^(−0.02) for 1282 and L^(−0.16) for 220, so z ≈ 1.0 and 1.2, where z = 2 would give L^(−1). Fits of ln ε_n against ln L (n = 10, 20) give z = 0.94 to 1.11 (1282), 0.99 to 1.26 (220) and 1.08 to 1.20 (free-fermion control).
- c_S. The ring entanglement fit (l ≥ L/4) gives 1.08 to 1.12 for 220 and 1.18 to 1.23 for 1282 at L = 21 to 26, without a trend, against 1.01 to 1.03 for both controls. Rings with a nearly degenerate ground state (220 at L = 14, 19, 1282 at L = 15, 25) are excluded. The overshoot matches the near-ferromagnetic controls of `d2_classification`, while c = 3/2 would need an 18 to 28% undershoot.

### 4. Open chains, and why the doubling estimate fails

`analyze_obc.py` fits S(l) over l ∈ [L/4, 3L/4] to the Calabrese–Cardy form alone (plain), with a bond-energy term (E_b), with a fitted A cos(2πkl + φ)/[(2L/π) sin(πl/L)]^p (osc, k in brackets), and after averaging over one period (smooth). "Full" means the field f on the end sites, as in `dmrg.py`, and "half" means f/2 there.

| point | ends | L | χ | plain | E_b | osc (k) | smooth |
|---|---|---|---|---|---|---|---|
| 220 | full | 96 | 128 | 0.96 | 0.99 | 0.95 (0.183) | 0.98 |
| 220 | half | 192 | 160 | 0.96 | 0.99 | 1.00 (0.186) | 0.99 |
| 1282 | full | 48 | 192 | 0.76 | 1.08 | 1.15 (0.302) | 1.02 |
| 1282 | full | 96 | 128 | 1.67 | 1.69 | 1.67 (0.429) | 1.64 |
| 1282 | half | 96 | 128 | 0.98 | 0.97 | – | 1.02 |
| 1282 | full | 192 | 96 | 0.96 | 0.95 | 0.96 (0.299) | 0.97 |
| 1282 | half | 192 | 160 | 1.10 | 1.09 | – | 1.14 |

For 220, S(l) oscillates at k = q ≈ 0.186 (a period of 5.4 sites) with a peak-to-trough amplitude of about 0.25. The central bond lies near a trough at L = 48 and near a crest at L = 96, so the central-cut difference picks up roughly 6 × 0.25/ln 2 ≈ 2. Every fit that removes the oscillation gives 0.90 to 1.02 at L = 48 and 0.95 to 1.00 at L = 96 and 192.

For 1282, the full-field L = 96 profile gives 1.67 under every fit. Halving the field on the two end sites brings the same window to 0.98, and doubling L with the full field brings it to 0.96 (χ = 96, still rising by about 0.04 per χ step). A bulk central charge cannot depend on a boundary field, so this is a boundary crossover reaching at least L/4 ≈ 25 sites into the chain. The half-field profile at L = 192 is not a clean logarithm either (local slopes 0.97, 0.5 and −0.5 over l = 24–48, 48–72, 72–96), and its window value of 1.1 comes from this shape. Open chains of 1282 at L ≤ 192 cannot fix c, and the boundary-free estimators give 1.

## Reproducing

numpy and scipy, plus numba for the ED (on Typhon, `PYTHONNOUSERSITE=1` with the default anaconda). Outputs are in `data/`, and the χ = 128 states (`*.npz`) only in `/scratch/lustre/jkudlerflam/d2_followup_220_1282/data/` on Typhon.

- `python3 run_fes.py MODEL --n N --mode damped --chis 8 12 16 24 32 48 64 96 128 --save_state` (MODEL is 220, 1282 or a control in `models.py`, 10 min to 2 h on 4 cores), then `analyze_fes.py` and `analyze_order.py` for sections 1 and 2.
- `python3 fdmrg.py MODEL --L 192 --chis 32 64 96 128 160 [--edgefield]` (about 1 h on 8 cores), then `analyze_obc.py 220 1282` for section 4.
- `../d2_classification/run_point.py NAME --k 6 --procs 16 --sizes ...` (`extra_points.json` defines the free-fermion control), then `ed_analysis.py 220 1282 xydm xxzm07_dm_h` for section 3.
- Slurm scripts are `fes_array.slurm` (`fes_tasks.txt`), `fes_round2.slurm`, `obc_packed.slurm` (`obc_tasks.txt`), `ed_big.slurm` and `ed_ctrl.slurm`.
- Typhon jobs (partition typhon, about 170 core-hours so far): 17511443 and 17511561 (VUMPS, with stalled tasks resubmitted and Ising run on the laptop), 17511439 (ED to L = 28, 16 cpus, 370 GB), 17511570 (ED controls), 17511574 and 17511763 (open chains, the latter at L = 384 and full end fields), 17511877 (VUMPS to χ = 256). ED at L = 27, 28, VUMPS at χ = 192, 256 and the L = 384 chains were still running when this was written, and their partial outputs are not used.

## Relation to the rest of the repository

This settles the two points flagged in the coordinator's note in `../d2_classification/README.md`. Their provisional "(b) c = 1" verdicts stand on boundary-free evidence, and "Where things stand" is unchanged. Damped VUMPS in commensurate cells could also settle the class (c) points 269, 1289 and 871.
