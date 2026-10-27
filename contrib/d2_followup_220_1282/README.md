# Points 220 and 1282 have c = 1

## Claim

Chains: H = Σ_j [σ_j^T J σ_{j+1} + f·σ_j] on qubits (Pauli matrices), nearest neighbour, one-site translation invariant, Hermitian, with (J, f) of points 220 and 1282 in `candidates_params.json` (copied from `../d2_classification`). The quantity is the Virasoro c of the critical theory.

Verdict, for both points: c = 1, a z = 1 incommensurate "floating" spiral (one compact boson).

| point | c from finite-entanglement scaling | spiral q (units of 2π) | spiral dimension x_q | ring ED c_S, L = 21 to 26 |
|---|---|---|---|---|
| 220 | 1.00 ± 0.01 (ξ up to 730) | 0.186 ± 0.002 | 0.16 | 1.08 to 1.12 |
| 1282 | 1.01 ± 0.02 (ξ up to 90) | 0.300 ± 0.004 | 0.33 | 1.18 to 1.23 |

Neither point is a counterexample to c ≤ d − 1. The doubling estimates 6ΔS/ln 2 of 1.93 (220) and 2.16 (1282) are artifacts of the open boundaries, explained in section 4.

## Status

Numerical. It would be upgraded by a quasiparticle-ansatz velocity, which would turn the ring energies into an independent Casimir estimate of c.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), a subagent of the coordinating session, 2026-10-27. No Codex/GPT was used.
- New code: `vumps.py` (multi-site VUMPS with a general MPO), `fdmrg.py` (open-chain DMRG with a general MPO), the analysis scripts. Ring ED reuses `../d2_classification/ed.py` and `run_point.py` unchanged.
- Checks: VUMPS reproduces −4/π for Ising and XX to 10⁻⁹ and the free-fermion energy of the XY + DM control (below) to 2·10⁻⁴ at χ = 64. A complex, rotated Ising chain gives the same numbers as the real one. At L = 96, `fdmrg.py --edgefield` reproduces the energy of `d2_classification/dmrg.py` to 3·10⁻⁹. Every estimator was run on controls with known c.
- Nothing has been rerun by anyone else.

## Method and results

### 1. Finite-entanglement scaling (boundary-free)

For an infinite MPS at bond dimension χ, a CFT gives S = (c/6) ln ξ + const, with ξ = −n/ln|λ₁/λ₀| from the n-site cell transfer matrix (Pollmann et al. 2009). I use VUMPS (Zauner-Stauber et al. 2018) with two-site null-space expansion, χ = 8 to 128, and report c_fit = 6 × slope of S against ln ξ over the last four χ.

A one-site cell fails here. The finite-χ state wants spiral order, which a one-site cell cannot hold, and VUMPS converges to non-injective states or not at all (the same happens for XX). I therefore scanned cells n = 1 to 37 at χ = 8 to 24 (`data/cellscan/cellscan.json`). The lowest energies come from cells commensurate with the spiral, n = 11 and 16 for 220 and n = 10, 20, 30 for 1282, with 23 and 13 next. Plain VUMPS also oscillates in these cells because of the soft sliding mode, so I use damped steps (A_C ← e^(−τH_AC) A_C, C ← e^(−τH_C) C, with τ adjusted so the variational energy never rises). A run stops when the gradient is below 10⁻⁸, or when it is below 10⁻³ and S has changed by less than 10⁻⁷ in 20 steps. For the candidates the gradient stalls at 10⁻⁷ to 10⁻⁶ (the sliding mode) while S and ξ have converged.

| model (exact c) | cell | χ (fit) | ξ (fit) | c_fit |
|---|---|---|---|---|
| Ising (1/2) | 1 | 32–96 | 590–4550 | 0.500 |
| XX (1) | 2 | 48–128 | 140–491 | 1.000 |
| XXZ Δ = 0.5 (1) | 2 | 64–128 | 174–427 | 1.008 |
| XXZ Δ = −0.5 (1) | 2 | 48–128 | 178–628 | 1.014 |
| XY + D_z + h_z, free fermions, no U(1), incommensurate (1) | 5 | 24–64 | 78–290 | 0.992 |
| cluster-Ising, three-site (3/2) | 2 | 48–128 | 41–109 | 1.500 |
| two decoupled XX chains, next-nearest neighbour (2) | 4 | 48–128 | 25–? | 1.94 |
| **220** | 11 | 48–128 | 208–730 | 1.002 |
| 220 | 16 | 24–64 | 48–164 | 1.03 |
| **1282** | 23 | 24–64 | 25–90 | 1.014 |
| 1282 | 10 | 24–48 | 26–64 | 1.00 |

The XY + D_z + h_z control (J = diag(1, 0.6, 0) plus D_z = 0.5, f = (0, 0, 0.6)) is the closest analogue. It has a large DM term, no conserved magnetization and Fermi points at 0.44π and 0.77π (`freefermion.py`). For 220 the pairwise slopes between consecutive χ from 24 to 128 are 1.027, 1.019, 1.005, 1.005, 1.001, 1.002. For 1282 (n = 23) from 24 to 64 they are 0.98, 1.00, 1.03, 1.00. The ξ–χ exponent κ is 1.28 (220) and 1.30 (1282), against 1.27 to 1.30 for the c = 1 controls and 1.02 for c = 3/2. Other cells (220 with n = 16, 21, 1282 with n = 10) agree, except one run (1282, n = 13) that jumped to a non-injective state at χ = 64 and is discarded. A c = 3/2 or 2 chain would show slopes near 1.5 or 2 already at ξ ~ 50, as the controls do.

### 2. Correlations in the infinite state

The cell magnetization splits into a uniform part m₀ and a spiral component M_q.

- 220: m₀ = 0.33. In the n = 11 cell the spiral "unlocks" above χ = 16 (all bonds carry the same S), and the leading transfer-matrix pair has phase 11q ≡ ±0.041, so q = 0.1855. In the n = 16 cell (q locked to 3/16) |M_q| ∝ ξ^(−0.16), so x_q ≈ 0.16.
- 1282: m₀ = 0.17, q = 0.30. The n = 10 cell (q = 3/10) is lower in energy than 7/23 and 8/27, and |M_q| ∝ ξ^(−0.33) in both the n = 10 and n = 23 cells.
- The same fit gives x = 0.250 for XX and 0.98 for the 2k_F density wave of the free-fermion control (exact 1/4 and 1).

So both chains have quasi-long-range spiral order, ⟨S₀·S_r⟩ ~ r^(−2x_q) cos(2πqr), with a single soft U(1) sliding mode. That is a c = 1 Luttinger liquid with K ≈ 1.5 (220) and 0.75 (1282), normalized so that XX has K = 1.

### 3. Ring ED with momentum resolution

`run_point.py` with six levels in each momentum sector, to L = 26 (220) and 25 (1282), plus the two c = 1 incommensurate controls to L = 24. L = 27 and 28 are still running (see Reproducing).

- Towers. The lowest excitation of every ring is a "rotor" pair at relative momentum ±N₀ (units 2π/L), with N₀/L = 0.17 to 0.20 (220) and 0.29 to 0.32 (1282). These are the charge ±1 states of the sliding mode at momentum ±2πq. Higher rotor states sit at ±2N₀ and ±3N₀ mod L. Their scaled energies L·ΔE/2π are 0.003 to 0.17, oscillating with qL mod 1 as a twisted boson does. There is one family, not the two independent families c = 2 would give. The descendants could not be isolated from this dense spectrum, so no velocity or Casimir c is quoted.
- z. The scaled levels L·ε_n (n ≤ 16, all momenta) are flat from L = 12 to 25 for 1282 and drift by about 30% for 220. A z = 2 spectrum would halve. Fits of ln ε_n against ln L give z = 0.96 to 1.12 (1282) and 0.99 to 1.26 (220) for n = 10 and 20, and 1.08 to 1.20 for the free-fermion control.
- c_S. The ring entanglement fit gives 1.08 to 1.12 (220, L = 21 to 26, no drift) and 1.18 to 1.23 (1282, L = 21 to 24). The free-fermion and twisted-XXZ controls give 1.01 to 1.03. Rings with a near-degenerate ground state (220 at L = 14 and 19, 1282 at L = 25) give c_S ≈ 2 and are excluded. The 10 to 20% overshoot matches the near-ferromagnetic controls of `d2_classification`, and it would have to be a 25% undershoot for c = 3/2.

### 4. Open chains, and why the doubling estimate fails

Windowed fits over l ∈ [L/4, 3L/4] (`analyze_obc.py`), plain, with a bond-energy term, with a fitted oscillation A cos(2πkl + φ)/[(2L/π) sin(πl/L)]^p, and after averaging over one period:

| point | end field | L | χ | plain | bond energy | oscillation (k) | smoothed |
|---|---|---|---|---|---|---|---|
| 220 | full | 96 | 128 | 0.96 | 0.99 | 0.95 (0.183) | 0.98 |
| 220 | half | 192 | 160 | 0.96 | 0.99 | 1.00 (0.186) | 0.99 |
| 1282 | full | 96 | 128 | 1.67 | 1.69 | 1.67 | 1.64 |
| 1282 | half | 96 | 128 | 0.97 | 0.97 | – | 1.01 |
| 1282 | half | 192 | 96 | 1.08 | 1.08 | – | 1.13 |

220: S(l) oscillates with k = q ≈ 0.186, a period of 5.4 sites and peak-to-trough amplitude about 0.25. The central bond happens to sit near a trough at L = 48 and near a crest at L = 96, so 6ΔS/ln 2 picks up 6 × 0.25/ln 2 ≈ 2. Every fit that removes the oscillation gives 0.95 to 1.00 at L = 48, 96 and 192.

1282: the L = 96 profile in `d2_classification` (full field on the end sites) gives c ≈ 1.67 under every windowed fit. Moving half the end field off the end sites (the sum-of-bond-terms convention) brings the same L = 96 window to 0.97. A bulk central charge cannot depend on a boundary field, so this is a boundary crossover over about L/4 ≈ 25 sites, of the kind the coordinator anticipated. At L = 192 the window gives 1.0 to 1.1, and the window-averaged doubling 96 → 192 gives 1.06.

## Reproducing

numpy and scipy, plus numba for the ED. On Typhon use `PYTHONNOUSERSITE=1` with the default anaconda. Data are in `data/`, and the χ = 128 states are only on Typhon (`/scratch/lustre/jkudlerflam/d2_followup_220_1282/data/*.npz`).

- `python3 run_fes.py MODEL --n N --mode damped --chis 8 12 16 24 32 48 64 96 128 --save_state` runs VUMPS. MODEL is 220, 1282 or a control in `models.py`. It takes 10 min (controls) to 1 h (n = 23, χ = 128) on 4 cores. `analyze_fes.py` and `analyze_order.py` then give sections 1 and 2.
- `python3 fdmrg.py MODEL --L 192 --chis 32 64 96 128 160 [--edgefield]` runs the open chain, 1 to 3 h on 8 cores. `analyze_obc.py 220 1282` gives section 4.
- ED: `../d2_classification/run_point.py NAME --k 6 --procs 16 --sizes ...`, with `extra_points.json` for the free-fermion control. `ed_analysis.py 220 1282 xydm xxzm07_dm_h` gives section 3.
- Slurm scripts: `fes_array.slurm` (tasks in `fes_tasks.txt`), `fes_round2.slurm`, `obc_packed.slurm` (`obc_tasks.txt`), `ed_big.slurm`, `ed_ctrl.slurm`.
- Typhon jobs, partition typhon:
  - 17511443 and 17511561: VUMPS. Some tasks were cancelled and resubmitted, and Ising was run on the laptop.
  - 17511439: ED to L = 28, 16 cpus and 370 GB.
  - 17511570: ED controls.
  - 17511574 and 17511763: open chains, including L = 384.
  - 17511877: VUMPS continued to χ = 256.

  Outputs of the jobs still running at the time of writing are not in this README.

## Relation to the rest of the repository

This resolves the two points flagged in the coordinator's note of `../d2_classification/README.md`. Their provisional "(b) c = 1" verdicts stand, now on boundary-free evidence. "Where things stand" does not change. The method, damped VUMPS in commensurate cells calibrated on the free-fermion XY + DM chain, could also settle the class (c) points 269, 1289 and 871.
