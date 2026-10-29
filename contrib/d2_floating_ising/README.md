# Can an Ising line be tuned inside a floating phase? A direct test at d = 2

## Claim

Class: one-site translation-invariant nearest-neighbour qubit chains H = Σ_j [σ_j^T J σ_{j+1} + f·σ_j] (Pauli matrices, the bond convention of `contrib/d2_classification/ed.py`), in the single-Z₂ stratum (C67) and the stratum with no uniform internal symmetry (C68). The question is the mechanism singled out by C68 Section 7 and left open by `contrib/d2_grid/` (direction (iii)): an independent neutral Ising sector going critical inside an incommensurate one-boson ("floating") liquid, which would give c = 3/2.

1. **[numerical]** On six two-parameter cuts through floating regions (1,728 ring points, L = 12 to 20, momentum resolved), no neutral level outside the boson tower softens anywhere in the liquid. In the Z₂ cuts (A, B, C, J; 1,156 points, 586 of them liquid-like) the lowest neutral, parity-preserving level at zero relative momentum sits at x ≥ 1.58 in units of the current level, and near x = 2 (the boson's J J̄) almost everywhere. An extra Ising energy field would put a level at x = v_I/v, for instance x = 1/2 in the Ising × XX control.
2. **[numerical]** Every level that a velocity-free screen flagged as an extra low q = 0 level (80 points, 42 of them stable in L) has the opposite Z₂ parity and a magnetization differing by ΔM ≈ 2 from the ground state. It is the adjacent charge sector of the same boson (a "moving zero mode"), and the exact U(1) line g = 0 of cut A produces the same flag with ΔM = 2 exactly.
3. **[numerical]** Finite-entanglement scaling (FES) at all 18 followed points (the flagged points, liquid points bordering Ising-like points, and three extra NS points) gives c between 0.93 and 1.09. The C68 anchor gives 0.986 (χ ≤ 64, ξ = 236). The calibration gives 1.000 for XX, 0.46 to 0.49 for Ising, and 1.54 for Ising × XX with the same code. Open chains at L = 64 and 128 give c between 0.92 and 1.10 at the 13 followed points, except two values below 0.85 that come from unconverged or end-dominated runs (listed below).
4. **[numerical]** So no coexisting Ising sector and no c > 1 was found. Where an Ising-like line (ring c_S ≈ 1/2) meets the liquid, as on the C67 junction family, the liquid side is c = 1 by both estimators.

This is a targeted search on six planes, not a scan of the strata. It does not exclude a narrowly tuned Ising line between grid points (spacing 1/16 of each cut range), or one whose velocity is so small that it stays above x ≈ 1.5 at L ≤ 20.

## Status

Numerical. The ring screen is calibrated on XX, XXZ and an exact U(1) liquid. The two thermodynamic estimators (FES and open-chain entanglement) are calibrated on Ising, XX and Ising × XX. Many FES runs had reached only χ = 16 to 48 when this was written (100-minute budget), so their c values carry the pair-slope scatter of about ±0.1 listed below. They are still 0.4 or more below the Ising × XX control at the same χ. The Ising × XX open-chain control at L = 128 had not converged (χ = 32, truncation 2 × 10⁻⁵), so the open-chain estimator is calibrated on Ising × XX only at L = 64, and through the U(1) DMRG of `contrib/d2_staggered/` (1.486 at L = 128).

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), a subagent of the c–d conjecture campaign, 2026-10-29. No Codex/GPT was used. Nobody else has rerun anything.
- Reused code, unchanged and imported: the momentum-resolved ring ED `contrib/d2_classification/ed.py`, the damped VUMPS `vumps.py` and the open-chain `fdmrg.py` of `contrib/d2_followup_220_1282/`. The cell-scan logic follows `contrib/d2_grid/fes_auto.py`, and the Ising × XX d = 4 term follows `contrib/d2_staggered/fes_point.py`. New code here only defines the families and cuts, drives these, and fits.
- Checks.
  - Ring detector controls (`ring_controls.py`, `data/ring_controls.json`). XX gives sJ = 3.97 to 3.98 (exact v = 4), x0 = 1.000 (the charge pair at q = 0) and neutral x0n = 2.000 at L = 16 to 20. XXZ at Δ = 0.5 gives x0 = 1.32 (exact 1/K = 4/3). The U(1) anchor in a field, (g, u, D, b) = (0, 0.4, 0.5, 1.05), exactly one boson, gives x0 = 0.67, 0.54, 0.41 at L = 16, 18, 20 (a charge level approaching a crossing) and neutral x0n = 2.04, 2.06, 2.07. This last control is the reason the raw q = 0 level cannot be used without the parity and magnetization test.
  - Level assignment (`levelchar.py`, `data/levelchar.json`). For eight flagged points the q = 0 low level has P opposite to the ground state and ΔM = 1.99 to 2.12 (Z₂ cuts), and the first level with the ground state's P and M is at s ≈ 2 sJ.
  - The anchor's ring spectrum reproduces the C68 continuum data. With sJ ≈ v the charge levels sit at s = 1.68, against v/(4K) = 5.444 × 0.308 = 1.68, and c_S(20) = 1.005.
  - FES calibration with the same code and χ list: XX 1.000 (χ ≤ 128, ξ = 491, κ = 1.28), Ising 0.46 to 0.49 depending on the window (χ ≤ 128, ξ = 5,568), Ising × XX 1.536 (χ ≤ 128, ξ = 97, κ = 0.97). Open chains: Ising 0.509 (L = 128), XX 0.980 (L = 128) and c_dbl = 0.995, Ising × XX 1.36 at L = 64 (χ = 64, truncation 5 × 10⁻⁶).

## Method and results

### Families and cuts (`models_fi.py`)

| cut | family | fixed | axes (17 × 17 grid) |
|---|---|---|---|
| A | Z₂: (Jx, Jy, Jz, D_z, f_z) = (1+g, 1−g, u, D, b) | u = 0.4, D = 0.5 | g ∈ [0, 0.8], b ∈ [0, 2.4] |
| B | same | g = 0.3, D = 0.5 | u ∈ [−0.8, 1.2], b ∈ [0, 2.4] |
| C | same | u = 0.4, b = 1.2 | g ∈ [0, 0.8], D ∈ [0, 0.8] |
| J | C67 eq. (13): (Jx, Jy, u, D, b) = (−1, −0.5, 0.25, D, −g_f) | | D ∈ [0, 0.4], g_f ∈ [1.6, 2.4] |
| N | C68 N family: S = diag(1+η, 1−η, 0.4), D = (η, 2η, 0.5), f = b (3, 4, 5)/√50 | | η ∈ [0, 0.4], b ∈ [0, 2.4] |
| P | C68 Section 7 plane through point 1282: J0 + s e_z e_zᵀ, f0 + t e_z | | s, t ∈ [−1, 1] |

The C68 anchor is (g, u, D, b) = (0, 0.4, 0.5, 0), a corner of cut A. Cuts A to C contain the pocket-to-Ising crossover of the free XY + DM + field chain (δ = 2g against 2D, and b near Jx + Jy), where an Ising line approaches the liquid. Cut J is the interacting Ising/chiral-liquid junction of C67. Cuts N and P are generic NS planes. In cut N the field strength is decoupled from η, so that η = 0 is not a symmetric point only when b ≠ 0.

### Stage 1. Ring spectra (`ringscan.py`, `analyze_ring.py`)

For each point and L = 12, 14, 16, 18, 20: the four lowest levels in every momentum sector, and the ground-state entropies at L = 16 and 20 with the chord fit c_S. With s = L(E − E0)/2π and q the momentum relative to the ground state, the diagnostics are

- sJ, the lowest level at q = ±2π/L. In one boson this is the U(1) current, x = 1, so sJ ≈ v.
- x0 = s0/sJ, the first excited level at q = 0. An Ising energy field is neutral, Z₂-even and always at q = 0, so an extra Ising sector puts a level at x0 = v_I/v (1/2 for Ising × XX, where v_I = 2 and v = 4). One boson gives neutral levels at x = 2 there, plus charged levels whose momentum happens to be commensurate on the ring.
- x0n (`levelscan.py`), the same with the extra requirement that the level has the ground state's Z₂ parity (Z₂ cuts) and |ΔM| < 1, M = Σ_j ⟨Z_j⟩.

A point is liquid-like if c_S(20) > 0.7 and sJ(20)/sJ(16) < 1.15, Ising-like if 0.3 < c_S(20) < 0.7 with the same sJ test, and flagged if liquid-like with x0 < 1.3 at L = 16, 18, 20 and not increasing.

| cut | points | liquid-like | Ising-like | flagged | flag stable in L | min over liquid of max(x0n(16), x0n(18)) |
|---|---:|---:|---:|---:|---:|---:|
| A | 289 | 203 | 0 | 16 | 3 | 1.97 |
| B | 289 | 209 | 3 | 21 | 14 | 1.58 |
| C | 289 | 134 | 6 | 32 | 15 | 1.99 |
| J | 289 | 40 | 21 | 0 | 0 | 10.4 |
| N | 289 | 201 | 0 | 11 | 10 | 0.66 (see below) |
| P | 283 | 92 | 1 | 0 | 0 | 1.90 |

ASCII maps of each cut are in `data/ring_maps.txt`. Cut A, for example, is liquid for g ≤ 0.6 at every b except a gapped band near b ≈ 0.4 to 0.7 for g ≥ 0.4, and gapped for g ≥ 0.65. The flagged points form a horizontal band at b ≈ 1.05 to 1.2 that starts on the exact U(1) column g = 0.

[numerical] In the Z₂ cuts the neutral channel is clean. No liquid point has a neutral q = 0 level below x = 1.58 at both sizes, and the flagged low levels are all the ΔM ≈ 2, opposite-parity charge state.

[numerical, weaker] In the NS cuts M is not conserved, so the |ΔM| < 1 test cannot separate neutral from charged levels. The N points with small x0n have ΔM between 0.3 and 0.9 and the same level pattern as the Z₂ charge states (a level at s ≈ 3, then the next at s ≈ 2 sJ), and they jump between L = 16 and 18 (for example 2.43 then 0.13 at η = 0.1, b = 1.05). For that reason four of them went to FES, which is the deciding test there.

### Stage 2. The two thermodynamic estimators (`fes_fi.py`, `obc_fi.py`, `analyze_mps.py`)

FES: damped VUMPS, a χ = 8, 16 cell scan over n = 1 to 10 (1 to 8 for the extra N points) choosing the lowest-energy injective cell, then χ = 8, 12, 16, ... with c = 6 dS/d ln ξ over the largest four χ. Open chains: two-site DMRG at L = 64 (χ ≤ 128) and 128 (χ ≤ 64 reached), full field on the end sites, S(l) = (c/6) ln[(2L/π) sin(πl/L)] fitted separately on even and odd bonds over L/8 ≤ l ≤ 7L/8, and c_dbl from the window means at L = 64 and 128.

| point | why | cell | χ reached | ξ | c (FES) | last pair slopes | c_obc L = 64 / 128 | c_dbl |
|---|---|---:|---:|---:|---:|---|---|---:|
| anchor (0, 0.4, 0.5, 0) | reference | 4 | 64 | 236 | 0.986 | 1.04, 0.95, 0.99, 1.01 | | |
| point 1282 (P, s = t = 0) | reference | 10 | 24 | 26 | 1.085 | 1.12, 1.18, 0.98 | | |
| B_14_10 (g, u, D, b) = (0.3, 0.95, 0.5, 1.5) | flag | 7 | 48 | 121 | 0.980 | 1.14, 0.95, 1.02, 0.98 | 0.93 / 0.95 | 1.00 |
| B_15_11 (0.3, 1.075, 0.5, 1.65) | flag | 7 | 32 | 69 | 1.007 | 0.88, 1.20, 0.92, 1.02 | 0.93 / 0.96 | 1.00 |
| B_15_12 (0.3, 1.075, 0.5, 1.8) | flag | 7 | 32 | 70 | 1.011 | 0.87, 1.17, 0.94, 1.03 | 1.01 / 0.98 | 0.97 |
| B_14_11 (0.3, 0.95, 0.5, 1.65) | flag | 7 | 48 | 122 | 0.980 | 1.14, 0.95, 1.03, 0.98 | 1.01 / 0.98 | 0.97 |
| B_13_10 (0.3, 0.825, 0.5, 1.5) | flag | 7 | 48 | 125 | 0.981 | 1.13, 0.95, 1.03, 0.97 | 1.00 / 0.98 | 0.97 |
| B_15_13 (0.3, 1.075, 0.5, 1.95) | flag | 7 | 32 | 71 | 1.012 | 0.85, 1.17, 0.94, 1.03 | 1.06 / 1.00 | 1.02 |
| C_08_12 (0.4, 0.4, 0.6, 1.2) | flag | 7 | 32 | 78 | 1.010 | 0.89, 1.05, 0.97, 1.04 | 1.02 / 0.99 | 1.02 |
| C_10_13 (0.5, 0.4, 0.65, 1.2) | flag | 8 | 32 | 101 | 0.992 | 0.98, 0.97, 0.99, 1.02 | 1.03 / 0.73 (χ = 32, unconverged) | |
| B_03_06 (0.3, −0.425, 0.5, 0.9) | junction | 9 | 32 | 171 | 0.982 | 1.07, 0.91, 1.00, 1.02 | 0.88 / 0.96 | 0.99 |
| B_02_06 (0.3, −0.55, 0.5, 0.9) | junction | 8 | 24 | 99 | 0.951 | 0.85, 1.10, 0.96 | 1.37 / 0.96 | 1.00 |
| B_03_07 (0.3, −0.425, 0.5, 1.05) | junction | 6 | 24 | 118 | 0.995 | 1.06, 0.91, 1.01 | 0.90 / 0.95 | 1.01 |
| J_10_06 (D, g_f) = (0.25, 1.9) | junction | | | | running | | 1.35 / 1.05 | 1.07 |
| N_01_07 (η, b) = (0.025, 1.05) | flag | 5 | 24 | 45 | 0.996 | 0.95, 1.04, 1.01 | 1.27 / 1.10 | 0.99 |
| N_01_08 (0.025, 1.2) | flag | 5 | 24 | 49 | 0.992 | 0.94, 1.04, 1.01 | 1.18 / 0.84 | 0.74 |
| N_00_07 (0, 1.05) | low x0n | 5 | 32 | 63 | 1.000 | 0.96, 1.04, 1.00, 0.96 | | |
| N_04_08 (0.1, 1.2) | low x0n | 5 | 16 | 33 | 0.977 | 0.96, 1.01 | | |
| N_09_07 (0.225, 1.05) | low x0n | 7 | 16 | 27 | 0.968 | 0.94, 1.03 | | |
| N_09_06 (0.225, 0.9) | low x0n | 7 | 12 | 20 | 0.933 | 0.93 | | |
| Ising × XX control | c = 3/2 | 2 | 128 | 97 | 1.536 | 1.61, 1.52, 1.52, 1.59 | 1.36 / unconverged | |

κ = d ln ξ/d ln χ is 1.36 to 1.44 at the followed points with at least four χ, against the c = 1 value 1.34 and 1.05 for c = 3/2 (0.97 measured on the control). The low-χ pair slopes scatter by about ±0.15 on every point, including the anchor. At χ = 16 to 48 the Ising × XX control already gives slopes 1.5 to 1.6, so a c = 3/2 point would separate at the χ reached.

[numerical] Every followed point is c = 1 on both estimators, within the calibrated scatter. The open-chain values at L = 64 that exceed 1.2 (B_02_06, J_10_06, N_01_07) fall to 0.96, 1.05, 1.10 at L = 128. The two values below 0.85 are C_10_13 (stopped at χ = 32 with truncation 4 × 10⁻⁷) and N_01_08 (a pinned spiral, the end-sensitivity documented for point 1282 in `contrib/d2_grid/`), and FES gives 0.99 at both.

[heuristic] The mechanism that would have produced c = 3/2 here needs a Z₂-even neutral mode whose gap closes at fixed q = 0 while the boson stays finite. In these planes the only soft q = 0 levels belong to the boson's own charge tower, and the Ising lines that do appear (cut J, the edges of cuts B and C) border the liquid without entering it. This matches the C67 picture, where the junction is a velocity-collapse endpoint rather than a crossing of two finite-speed sectors.

### What this does not cover

- Six planes, 17 × 17 each. A line narrower than the grid spacing, or one present only in other directions of the 4-dimensional (Z₂) or 8-dimensional (NS) strata, is not tested. The C68 response-matrix step (retuning mixed relevant couplings once a neutral branch is seen) was never triggered, because no neutral branch was seen.
- The neutral-level test is clean in the Z₂ cuts but only indirect in the NS cuts, where the decision rests on FES at four points per cut and on the absence of flags in cut P.
- FES runs at the followed points reached χ = 12 to 64 (ξ = 20 to 236) at the time of writing. Several were still running and write into `data/fes/` as each χ finishes.
- A c > 1 sector whose velocity is much smaller than the boson's would appear at x0n = v_I/v ≪ 1 and would have been seen. One with v_I ≳ 0.8 v would sit at x0n ≳ 0.8 and would still have been seen in the Z₂ cuts (minimum 1.58). Velocities v_I > 1.5 v are not excluded by the ring test alone, only by the entanglement estimators at the followed points.

## Reproducing

numpy, scipy and numba (on Typhon, `PYTHONNOUSERSITE=1` with the default anaconda). Run from this folder. `../d2_classification/ed.py` and `../d2_followup_220_1282/{vumps,fdmrg}.py` are imported, not copied.

- Ring scans: the 24 lines of `tasks_ring.txt` (`sbatch --array=0-23 cmd_array.slurm tasks_ring.txt`), each `python3 ringscan.py CUT --sizes 12 14 16 18 20 --ent 16 20 --procs 8 --start S --stop S+73`; `tasks_ring2.txt` ran the second half of each chunk in parallel (duplicates are removed by tag). About 30 s per point per core at L ≤ 20. Then `python3 analyze_ring.py --maps` (writes `data/ring_summary.json`, stored here gzipped, gunzip before use), `python3 levelscan.py --L 16 18 --procs 9` (`data/levelscan.json`), `python3 levelchar.py TAG ... --L 16 18`, `python3 ring_controls.py`.
- Selection: `python3 select_fes.py --nflag 10 --nj 3` writes `tasks_fes.txt` and `tasks_obc.txt`; the four extra NS points are `tasks_fes_n.txt`.
- Calibration: `tasks_calib.txt`. Production: `sbatch --cpus-per-task=4 --mem=8G --time=01:30:00 --array=0-15 cmd_array.slurm tasks_fes.txt`, the same with `tasks_obc.txt` (0-27) and `tasks_fes_n.txt` (0-3). Then `python3 analyze_mps.py` (writes `data/mps_summary.json`).
- Data here: `data/ring/` (every ring point, all levels and entropies), `data/fes/`, `data/obc/` (S on every bond), and the summaries above. No MPS states were saved.
- Typhon usage, partition typhon: jobs 17526129 and 17526167 (ring scans), 17526153 and 17526154 (calibration), 17526197, 17526198, 17526242 (FES and open chains), about 60 core-hours in total.

## Relation to the rest of the repository

This is the direct test listed as not attempted in `contrib/d2_grid/README.md` (direction (iii)) and named as the next useful computation in `codex/work_C68/RESULT.md` Section 7. It complements `contrib/d2_staggered/` (no neutral Ising sector on the staggered branch) with the same negative answer in the single-Z₂ and symmetry-free strata, on the planes tested. Natural next steps would be (i) the C68 K-stratum tilted-field continuation, (ii) larger χ at the followed points, and (iii) an NS neutral-level test using the magnetization along the local spiral axis instead of Σ Z.
