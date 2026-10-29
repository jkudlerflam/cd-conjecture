# Microscopic DMRG of the d = 13 mobile bond at its smallest members

Mechanism check on the full thirteen-state chain of `codex/work_C88/bond13.py` (eleven material letters plus V_0, V_1), at the smallest members of the same bond family, with couplings chosen so that every scale is representable in double precision. Target c = n + 1 for the mobile member with n payload flavours (p = n + 1 data sites, cell length b = p + 1), and c = n for the control with the vacancies frozen out (full occupancy).

## Claim

1. [numerical, exact diagonalization of frozen rings] At the smallest members the logical (payload) bandwidth of the literal bond is far below every other scale unless the token cost Δ is pushed close to resonance. At the default diagnostic point of `bond13.params` (D = 4, g = 0.07, λ = 0.1, t_d = 0.003) the n = 1 bandwidth on a three-cell ring is 6 × 10⁻¹⁴, below double-precision resolution of the chain's O(1) energies. Resolving it needs D = 0, Δ ≲ 1 and g, λ ≳ 0.2. At D = 0, Δ = 0.5, g = λ = 0.3, t_d = 0.2 (the point used below) the n = 1 bandwidth is 9.2 × 10⁻⁵ against a cell gap of 0.28 and a next-level gap of 0.19. The n = 2 member at the same point has bandwidth 3.7 × 10⁻⁶ on a two-cell ring.
2. [numerical, not yet available] Calibrated central charges of the mobile member and its control. The DMRG runs are launched (see Status). No c value is claimed in this file yet.
3. [heuristic] Since the mobile material bond is multiplied by J and by the contact probability (about 0.3 at half occupancy), the n = 1 payload velocity at J/t = 0.5 is of order 10⁻⁵ t, against a vacancy velocity of order t. Finite-χ DMRG drops the slowest mode first (this was the main systematic in `../vacancy_n2_dmrg/`), so I expect the open-chain estimators to need large χ before the payload mode shows up, and the n = 2 member to be out of reach at any χ I can afford.

## Status

Partial. The code is written and the frozen-ring scales are measured. The open-chain DMRG runs had completed only their first χ = 64 sweep when this file was written, so no calibrated c is available for any member or control. Not done: VUMPS or iDMRG finite-entanglement scaling (no complex symmetric uMPS code was adapted in the time available), and the cell-length selection beyond what the fixed (N, Mk) sector enforces.

Two points to check when the runs finish. The control `ctl1_L45` reported E = +117.23 after its first sweep while `ctl1_L30` reported -5.64, so the L = 45 start may be stuck in a high-energy configuration. And the redteam colour field t²/(JK) is 2 × 10⁻⁴ at K = 10⁴, J = 0.5, which is comparable to the n = 1 logical bandwidth (about J α_c × 10⁻⁴). The redteam conditions F1 and F2 (J/t ≪ 0.0056, D J/t small) are met here in the second sense (D = 0) but not the first. At b = 3 the selector window is much wider than at b = 15, and in a fixed (N, Mk) sector μ_sel drops out, so F1 only enters through wrong-length pair defects, which the runs do not measure directly.

## How it was established and checked

- Author: Claude (Opus 5.5) subagent for Jonah Kudler-Flam, 2026-10-29. No Codex call (repository-bound numerics). Nobody else has rerun anything.
- `bond13.py` is imported read-only. Nothing outside this folder was modified.
- `cqdmrg.py` is `../vacancy_n2_dmrg/qdmrg.py` with complex dtypes and Hermitian conjugates in the environments, Lanczos inner products and densities (the d = 13 bond is complex through φ_A, ψ_A and the port phase θ). Conserved charges are Q (the bond13 charge, twice the logical fermion number on the code), the occupied number N and the marker number Mk.
- `check_ed.py` compares `cqdmrg.py` with exact diagonalization of open chains of 6 to 8 sites in the seed's full connected sector (finite K = 30, illegal words included). It was running on the Typhon login node when this file was written. Its output (`data/check_ed.json`) should be confirmed before any DMRG number is used.

## Method and results

**Scales** (`scales.py`, `data/scales.jsonl`). Frozen material ring of M cells (full occupancy, so markers and cell lengths are fixed), projected on legal words (K = ∞), from the all-zero vacuum seed and from a seed with one bit set. The M lowest one-fermion levels form the logical band. Rows have μ_sel = 0 and f = 0, and the cell gap is the first excitation in the charge-0 sector.

| p | M | D | Δ | g | λ | t_d | cell gap | logical bandwidth | next gap |
|---|---|---|---|---|---|---|---|---|---|
| 2 | 3 | 4 | 2.59 (default) | 0.07 | 0.1 | 0.003 | 1.414 | 5.8e-14 | 1.414 |
| 2 | 3 | 0 | 2 | 0.2 | 0.2 | 0.2 | 1.409 | 6.4e-10 | 1.409 |
| 2 | 3 | 0 | 1 | 0.2 | 0.2 | 0.2 | 1.298 | 6.0e-8 | 1.246 |
| 2 | 3 | 0 | 0.5 | 0.2 | 0.2 | 0.2 | 0.303 | 1.5e-5 | 0.249 |
| 2 | 3 | 0 | 0.5 | 0.3 | 0.3 | 0.2 | 0.279 | 9.2e-5 | 0.193 |
| 2 | 3 | 0 | 0.5 | 0.5 | 0.5 | 0.2 | 0.230 | 1.8e-3 | 0.075 |
| 3 | 2 | 0 | 0.5 | 0.3 | 0.3 | 0.2 | 0.444 | 3.7e-6 | (band continues) |

[numerical] The bandwidth grows by about four orders of magnitude between Δ = 2 and Δ = 0.5, which is the expected behaviour of a high-order virtual process (two ports, two conversions and the token transit) with denominators set by the token cost. The hierarchy cell gap ≫ bandwidth survives at Δ = 0.5 with g = λ = 0.3 (ratio about 2000), but only about 40 at g = λ = 0.5.

**Open-chain DMRG** (`run_mobile.py`, `tasks.txt`). Couplings D = 0, Δ = 0.5, g = λ = 0.3, t_d = 0.2, K = 10⁴, J = 0.5 (one run at J = 1), t = 1, flux and θ at the `bond13.params` defaults. Product start: cells [R_0, R_a, …, E] with the logical bits spread evenly and V_0 vacancies after markers (the unique legal colouring). Occupancy N/L = 5/8, payload filling Nlog/(Mk(p - 1)) = 0.4 for n = 1 and 0.35 for n = 2, neither half filling nor a low-order commensurate density.

| run | p | L | N | Mk | Nlog | target c |
|---|---|---|---|---|---|---|
| n1_L48, n1_L72, n1_L96 | 2 | 48, 72, 96 | 30, 45, 60 | 10, 15, 20 | 4, 6, 8 | 2 |
| n1J1_L72 | 2 | 72 | 45 | 15 | 6 | 2 (J = 1) |
| ctl1_L30, L45, L60 | 2 | 30, 45, 60 | = L | 10, 15, 20 | 4, 6, 8 | 1 (frozen control) |
| n2_L64 | 3 | 64 | 40 | 10 | 7 | 3 |
| ctl2_L40 | 3 | 40 | 40 | 10 | 7 | 2 (frozen control) |

χ grows through 64, 128, 192, 256, 384, 512. First-sweep timings on the loaded login node were 3 to 5 minutes per sweep at χ = 64, so χ = 512 will take days.

**Estimators** (`analyze.py`). (A) Difference estimator of `../vacancy_n2_dmrg`: c = m + slope of S_model(ℓ) - S_ref(ℓ) in the chord variable, with S_ref the exact profile of m free open chains carrying the sector's numbers (m = 2 for the mobile n = 1 member: the vacancy fluid with N particles and the payload with Nlog on L sites), and Friedel terms at 2πN/(L + 1), 2πNlog/(L + 1) and 2πMk/(L + 1). (B) Ratio m c_fit(model)/c_fit(ref) with the same fit. (C) Half-chain pairs between sizes, calibrated by the same pair on the reference. For the frozen controls the entropy is sampled at cell boundaries and compared with a free chain of Mk sites. For n = 2 the payload Fermi seas are not known in closed form, so `analyze.py` reports only the raw fit for the n = 2 control and does not build an n = 2 mobile reference. The vacancy density and markers per occupied site are printed for every run (both are fixed by the sector, so they serve as sanity checks rather than measurements).

## Reproducing

numpy and scipy only, `python3 -B` from this folder. `bond13.py` is found at `../../codex/work_C88/`.

- `python3 -B scales.py --p 2 --M 3 --D 0 --Delta 0.5 --g 0.3 --lam 0.3 --td 0.2` prints one row of the scales table (seconds for M = 3, p = 2).
- `python3 -B check_ed.py` writes `data/check_ed.json` (minutes to tens of minutes, dominated by the closure at L = 8).
- `python3 -B run_mobile.py NAME <line of tasks.txt>` writes `data/obc/NAME.json` after every χ. `job.slurm` is the array over `tasks.txt`.
- `python3 -B analyze.py` prints every estimator and writes `data/summary.json`.

Runs in flight on Typhon (2026-10-29, about 09:10 local), all under `/scratch/lustre/jkudlerflam/d13_microscopic_dmrg/contrib/d13_microscopic_dmrg/`:
- login node at nice 19: tasks 1, 2, 4, 5, 8 (`n1_L48`, `n1_L72`, `ctl1_L30`, `ctl1_L45`, `n2_L64`) and `check_ed.py`, logs in `logs/login_*.out` and `logs/check_ed.out`;
- SLURM array 17526426, tasks 3, 6, 7, 9 (`n1_L96`, `ctl1_L60`, `n1J1_L72`, `ctl2_L40`), pending on QOSMaxJobsPerUserLimit, logs in `logs/d13dmrg_*.out`.

To collect: `scp 'typhon:/scratch/lustre/jkudlerflam/d13_microscopic_dmrg/contrib/d13_microscopic_dmrg/data/obc/*.json' data/obc/` and `scp typhon:/scratch/lustre/jkudlerflam/d13_microscopic_dmrg/contrib/d13_microscopic_dmrg/data/check_ed.json data/`, then `python3 -B analyze.py`. The JSON files hold energies, bond entropies and charge densities only (no MPS), well under 1 MB each.

## What remains

(i) Confirm `check_ed.json`, then read the calibrated c of the n = 1 member and its control at converged χ. (ii) A VUMPS or iDMRG FES run on the period-b frozen control and on the mobile chain, which needs a complex, symmetry-resolved uMPS code. (iii) The n = 2 member seems out of reach with this bond at any resolvable coupling, because its bandwidth is a further factor of about 25 below n = 1. A resolvable n = 2 check would need either Δ even closer to the token threshold (where the code picture weakens) or a reduced logical model of the kind used in `../vacancy_n2_dmrg`.
