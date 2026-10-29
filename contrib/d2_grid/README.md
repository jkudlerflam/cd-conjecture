# A calibrated survey of the generic d = 2 residual

## Claim

Class: periodic, one-site translation-invariant, nearest-neighbour Hermitian qubit chains,
H = Σ_j [σ_j^T J σ_{j+1} + f·σ_j], with J = diag(Jx, Jy, Jz) + A(D), A_ab = ε_abc D_c, in the principal frame of the symmetric part (the C6 convention, `codex/work_C6_d2/`). The quantity is the total Virasoro central charge c of the infrared CFT.

The survey covers the four residual strata of C6 that are not controlled analytically (no uniform internal symmetry, exactly one unitary Z₂, time reversal only, and the plane-reflection antiunitary), and the staggered-field box a(XX − YY) + uZZ + bZ, a = 1, −1 < u < √(1 + b²/4), of `contrib/d2_gap_reduction/`.

1. **[numerical]** No point with c > 1 was found. Every point that small-size estimators flagged as possibly above 1 (173 points, out of about 103,000 sampled) was followed to the infinite chain with finite-entanglement scaling (FES). Each one is either a c = 1 point, a gapped point with a long but finite correlation length, or one of two points whose FES slope stays below 1 and falls slowly as χ grows (to 0.80 and 0.50 at the largest χ).
2. **[numerical]** The largest FES central charge among converged runs is c = 1.04 ± 0.05 (point Z2_21_295, χ ≤ 64, a single-Z₂ point). The 132 flagged points that pass the c = 1 criterion have c_fit = 1.001 with a standard deviation of 0.012 across points (range 0.965 to 1.044), against 1.000 for XX, 0.992 for the free-fermion spiral control, 1.500 for cluster-Ising and 2.01 for two decoupled XX chains with the same code.
3. **[numerical]** This includes the neighbourhoods of 880 "junctions", places on the sampled 2-spheres where an Ising-like line (ring c_S ≈ 1/2) touches both a c ≈ 1 region and a gapped region. None of them carries c > 1. Twenty-one junction-patch points survived both screens, and FES gives c = 1 for nine of them and a gap for the other twelve.
4. **[numerical]** In the staggered-field box the ring maximum is c_S(16) = 1.27, on the nonconformal boundary u = r. FES at the C63 overshoot point (u, b) = (1.4, 2.0) gives c = 1.02, and at three other points near the boundary 0.99 to 1.04. Interior points tested are gapped.
5. **[heuristic]** The large ring values (c_S up to 3.6 at L = 16, 20% of the symmetry-free points above 1.1) come from near-classical spirals. Their finite-ring ground states are momentum eigenstates built from a rotor tower of the sliding mode, which adds a "zero-mode" entanglement that grows like ln l and has nothing to do with c.

This is not a proof. It supports "no counterexample at d = 2" in the generic residual at the resolution described below.

## Status

Numerical, with the sampling limits stated in "What the survey does not cover". The estimator that decides every verdict is infinite-chain FES with damped VUMPS, calibrated on five exactly known chains. The ring and open-chain estimators are used only as screens, and the survey shows quantitatively why they cannot decide on their own in this family.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), a subagent of the coordinating session, 2026-10-28. No Codex/GPT was used. Nobody else has rerun anything.
- Reused code, unchanged: the momentum-resolved ring ED `contrib/d2_classification/ed.py` (its self-test matches dense diagonalization in every sector to 10⁻¹³), and `vumps.py`, `fdmrg.py`, `analyze_obc.py`, `models.py` from `contrib/d2_followup_220_1282/`. New code here only builds parameters, drives these, and fits.
- Checks.
  - The stratum construction was checked by brute force. `strata.symmetry_report` tests every principal-axis π rotation and every improper involution that preserves diag S against the C6 conditions (G_u, G_a). Random points of each stratum have exactly the intended symmetry (NS none, Z₂ only R_z, TR only T, RA only K_z), and the staggered family has R_z, K_x, K_y.
  - Every estimator was run on chains with known c (tables below). Ring entanglement on Ising, XX and the free-fermion spiral, open-chain fits on Ising, XX and the free-fermion spiral, FES on Ising, XX, the free-fermion spiral, cluster-Ising (3/2) and two decoupled XX chains (2).
  - The screening thresholds were tested on a control group of 300 random points with ring c_S ≈ 1 (4 flagged, against 142 of the 2,625 ring-flagged points).
  - FES verdicts were cross-checked against changing the cell. Four of the points flagged as "anomalous" in the first FES pass gave c = 1.00 to 1.01 once the cell scan was extended to n ≤ 24 (table 5).
  - The C63 run (`codex/work_C63/`) names no candidate in the generic residual. Its strongest staggered-family overshoot, (u, b) = (1.4, 2), is included and gives c = 1.

## Method and results

### Stage 0. Sampling

Because positive rescaling is irrelevant, each stratum is a sphere in its coordinate subspace of v = (Jx, Jy, Jz, Dx, Dy, Dz, fx, fy, fz). I take random great 2-spheres through each stratum (three random orthonormal vectors in its subspace, `strata.sphere_basis`) and put N Fibonacci points on each. A generic 2-sphere meets a codimension-one critical surface in curves and a codimension-two junction of two such surfaces in isolated points, so random 2-spheres sample junctions as well as phases.

| stratum | coordinates | spheres × points | ED sizes |
|---|---|---|---|
| NS, no uniform internal symmetry | all nine | 40 × 600 | 12, 16 |
| Z2, one unitary Z₂ about z | Jx, Jy, Jz, Dz, fz | 25 × 400 | 12, 16 |
| TR, time reversal only (f = 0) | Jx, Jy, Jz, Dx, Dy, Dz | 25 × 500 | 12, 16 |
| RA, plane-reflection antiunitary | Jx, Jy, Jz, Dz, fx, fy | 25 × 500 | 12, 16 |
| staggered box, a = 1 | u ∈ [−1, 1.6], b ∈ [0.05, 2.5] | 27 × 25 grid | 12, 16 |
| junction patches | 7 × 7 tangent grid, half-width 0.2 rad | 880 × 49 | 12, 16 |

Per point (`scan.py`): the two lowest levels in every momentum sector, the ground state's block entropies, and the ring chord fit c_S(L) = 3 × slope of S(l) against ln[(L/π) sin(πl/L)] over l ≥ L/4.

Fibonacci spacing is about 0.15 rad on a 600-point sphere, and the patches refine a 0.4 rad square around each junction candidate to 0.07 rad.

### Stage 0 results. What the rings show

| stratum | points | c_S(16) < 0.25 | max c_S(16) | c_S(16) > 1.1 | > 1.3 | min(c_S(12), c_S(16)) > 1.2 |
|---|---:|---:|---:|---:|---:|---:|
| NS | 24,000 | 14,206 (59%) | 2.66 | 4,749 | 958 | 1,461 |
| Z2 | 10,000 | 6,195 (62%) | 2.53 | 924 | 192 | 180 |
| TR | 12,500 | 3,317 (27%) | 3.06 | 2,295 | 460 | 572 |
| RA | 12,500 | 8,606 (69%) | 3.60 | 1,200 | 287 | 412 |

The c_S(16) < 0.25 column is an approximate gapped fraction. The TR stratum is mostly gapless, as expected from the Lieb–Schultz–Mattis constraint for a Kramers spin-1/2 with translation, which rules out a unique gapped symmetric state. The maps show extended c ≈ 1 regions (floating spirals) and gapped regions, separated by lines.

Ring values above 1.3 are common and do not mean c > 1. A typical example is TR_20_263 with c_S = 3.10 at L = 12 and 3.06 at L = 16. Its single-site entropy is ln 2 (zero local magnetization in a chain with no symmetry that forbids it) and its low spectrum is a tower of rotor states at momenta that follow qL. This is the momentum-projected spiral. A translation-averaged product spiral already has S(l) growing like ln l on a ring, which the chord fit reads as c ≈ 3. The same mechanism makes the ring values jump with L (NS_11_249 gives 1.19, 1.96, 3.50 at L = 18, 20, 22; TR_19_122 gives 1.90, 2.12, 3.40). The ring calibration is fine away from spirals.

| control (exact c) | c_S(12) | c_S(16) | c_S(20) |
|---|---:|---:|---:|
| Ising (1/2) | 0.501 | 0.501 | 0.500 |
| XX (1) | 1.016 | 1.007 | 1.004 |
| XY + D_z + h_z, free-fermion spiral (1) | 1.048 | 1.053 | 1.012 |
| XXZ Δ = −0.95 (1) | 1.168 | 1.104 | 1.068 |

### Stage 1. Open-chain screen

All 2,625 points with min(c_S(12), c_S(16)) > 1.2, a control group of 300 random points with 0.8 < c_S < 1.2 at both sizes, and seven known chains went through two-site DMRG on open chains of L = 32 and 64 (χ = 24 then 48, full uniform field on the end sites, `obcscreen.py`). Open ends pin the spiral and remove the rotor zero mode. With x(l) = (1/6) ln[(2L/π) sin(πl/L)] and the window l ∈ [L/4, 3L/4], it records c_eb(64), the window fit with the bond-energy correction of Laflorencie et al., and c_dbl = 6[S̄(64) − S̄(32)]/ln 2, S̄ the window mean. A point is flagged when min(c_dbl, c_eb(64)) > 1.2.

| chain (exact c) | c_dbl | c_eb(64) |
|---|---:|---:|
| Ising (1/2) | 0.528 | 0.508 |
| XX (1) | 0.995 | 0.966 |
| free-fermion spiral (1) | 0.957 | 0.967 |
| XXZ Δ = −0.95 (1) | 1.035 | 0.995 |
| XXZ Δ = 0.5 (1) | 0.972 | 0.932 |
| point 220 of d2_classification, c = 1 floating spiral | 0.585 | 0.910 |
| point 1282, c = 1 floating spiral | 1.132 | 0.685 |

| group | points | median c_dbl | median c_eb(64) | flagged |
|---|---:|---:|---:|---:|
| ring-flagged | 2,625 | 1.00 | 0.98 | 142 |
| ring c ≈ 1 control | 300 | 0.99 | 0.95 | 4 |

The open chain removes the ring overshoot for 95% of the ring-flagged points. The last two control rows show the remaining weakness: on spirals pinned by the ends the open-chain values scatter by ±0.4 around 1 at these sizes.

The 880 junction candidates (an I point, 0.3 < c_S < 0.75 at both sizes, whose eight nearest neighbours include a B point with 0.8 < c_S < 1.3 at both sizes and a G point with c_S(16) < 0.25) were refined by a 7 × 7 patch at L = 12 and 16. The best point of each patch (largest min(c_S(12), c_S(16))) was kept if above 1.1, which gave 431 points, and these went through the same open-chain screen. Twenty-one were flagged.

### Stage 2. Finite-entanglement scaling, the deciding estimator

For each flagged point (146 from the spheres, 21 from the junction patches, 6 from the staggered box), `fes_auto.py` runs damped VUMPS at χ = 8 and 16 in every cell n = 1 to 12 (later 16, and 24 for the anomalous points), takes the lowest-energy converged injective cell, and grows χ to 64 (128 or 256 for selected points). For a CFT, S = (c/6) ln ξ + const, with ξ = −n/ln|λ₁/λ₀| of the cell transfer matrix. I report c_fit = 6 × slope over the last four χ, the pair slopes, and κ = d ln ξ/d ln χ (theory for a z = 1 CFT 1.34, 1.05, 0.87 for c = 1, 3/2, 2).

Calibration, same code and settings.

| chain (exact c) | cell | χ | ξ at largest χ | c_fit | κ |
|---|---|---|---:|---:|---:|
| Ising (1/2) | 1 | 8 to 128 | 3,460 | 0.48 (pairs 0.46 to 0.52) | 1.62 |
| XX (1) | 2 | 8 to 128 | 491 | 1.000 | 1.28 |
| free-fermion spiral, XY + D_z + h_z (1) | 5 | 8 to 64 | 290 | 0.992 | 1.34 |
| cluster-Ising, three-site (3/2) | 2 | 8 to 128 | 109 | 1.500 | 1.02 |
| two decoupled XX chains, NNN (2) | 4 | 8 to 128 | 56 | 2.01 (spread 0.17) | 0.83 |

Verdicts (`summary.py`). "c = 1" means |c_fit − 1| < 0.06 over at least five converged χ and 1.15 ≤ κ ≤ 1.45. "c = 1, one jump" means the median of the last five pair slopes is within 0.06 of 1 with only one slope off by more than 0.3, which happens once when a spiral unlocks inside the cell (as for point 220 in `d2_followup_220_1282`). "Gapped" means the last pair slopes fall below 0.5 and keep falling, with S saturating while ξ grows.

| source | points | c = 1 | c = 1, one jump | gapped | open |
|---|---:|---:|---:|---:|---:|
| NS spheres | 94 | 89 | 1 | 4 | 0 |
| Z2 spheres | 2 | 2 | 0 | 0 | 0 |
| TR spheres | 16 | 8 | 1 | 7 | 0 |
| RA spheres | 34 | 21 | 3 | 8 | 2 |
| junction patches (all strata) | 21 | 9 | 0 | 12 | 0 |
| staggered box | 6 | 3 | 1 | 2 | 0 |
| total | 173 | 132 | 6 | 33 | 2 |

For the 132 c = 1 points, c_fit = 1.001 ± 0.012 (mean and standard deviation across points), κ = 1.31 ± 0.04, cells n = 1 to 23, and ξ = 51 to 863 at the largest χ. The largest single value is Z2_21_295, c_fit = 1.044 with a spread of 0.05 over shifted windows (pairs 1.31, 1.36, 0.92, 1.05, 1.05, 1.02, ξ ≤ 51). The gapped points have ξ = 5 to 65 at χ = 64 to 128.

### Stage 3. The points that needed more

These are every point whose first FES pass was not clean, or whose screen values were extreme.

| point | ring c_S, L = 12/16 [18/20/22] | open chain L = 96 and 192, c_eb | first FES | final FES | verdict |
|---|---|---|---|---|---|
| NS_20_482 | 1.40, 1.37 [1.35, 1.22, 1.22] | 0.78, 0.81 | n = 12, 1.55 (χ ≤ 24) | n = 15, 1.001 (χ ≤ 96, ξ = 278) | c = 1 |
| NS_30_144 | 1.48, 1.23 [1.42, 1.16, 1.38] | 2.00, 1.10 | n = 5, pairs 1.5 to 0.9 | n = 14, 1.002 (χ ≤ 96, ξ = 312) | c = 1 |
| TR_19_122 | 1.67, 1.75 [1.90, 2.12, 3.40] | 0.99, 0.14 | n = 2, pairs 1.4 to 1.0 | n = 23, 1.011 (χ ≤ 96) | c = 1 |
| NS_34_153 | 1.41, 1.36 [0.00, 1.30, 1.29] | 0.39, 1.00 | n = 4, pairs 1.8 to 1.0 | n = 23, 0.965 (χ ≤ 96) | c = 1 |
| NS_20_461 | [1.21, 1.10, 1.09] | 1.27, 1.13 | n = 11, χ ≤ 24 | n = 13, 1.005 (χ ≤ 96) | c = 1 |
| NS_36_224 | [1.36, 1.25, 1.26] | 0.79, 0.96 | n = 11, χ ≤ 24 | n = 23, 0.972 (χ ≤ 64) | c = 1 |
| NS_23_470 | [1.63, 1.51, 1.42] | 0.46, 0.74 | n = 9, one jump | n = 20, 1.014 (χ ≤ 96) | c = 1 |
| NS_11_249 | 2.47, 2.22 [1.19, 1.96, 3.50] | 1.80, 2.32 (unconverged) | n = 2, pairs 2.2 to 3.0 (χ ≤ 64) | n = 2 and n = 4 (same state), pairs 1.01, 0.58, 0.45, 0.30 at χ = 64 to 192, ξ = 50 | gapped |
| RA_24_106 | 2.11, 1.22 | 0.87, 1.41 | n = 3, pairs up to 5.0 | pairs 0.60, 0.42, 0.27, 0.13, 0.06 at χ = 48 to 192, ξ = 59 | gapped |
| RA_23_218 | 1.59, 1.56 [1.57, 1.50, 1.46] | −1.6, 1.42 | n = 2 | pairs 0.24, 0.14, 0.07, 0.02 (χ ≤ 128, ξ = 65) | gapped |
| RA_23_433 | 1.39, 1.52 [1.57, 1.53, 1.52] | −4.5, 0.77 | n = 2 | pairs 0.74, 0.82, 0.73, 0.50 at χ = 32 to 128, ξ = 240, κ = 0.86 | open, below 1 |
| RA_5_88 | [1.39, 1.27, 1.20] | 1.02, 1.03 | n = 2 | pairs 0.92, 0.89, 0.90, 0.84, 0.80 at χ = 32 to 192, ξ = 198 | open, below 1 |

NS_11_249 deserves a remark, since every small-size estimator put it near c = 2 and its first FES pass did too (pair slopes 2.2, 1.7, 2.2, 3.0, 2.4 for χ = 8 to 48, with ξ growing only from 10 to 23, κ = 0.36). From χ = 64 on the pair slopes fall (1.01, 0.58, 0.45 at χ = 64, 96, 128) while ξ grows slowly to 41. A CFT of any c gives a constant slope with κ near the theory value, so this is the crossover of a gapped chain with a long correlation length, and the small-χ slopes measure the growth of short-distance entanglement before saturation. At χ = 192 the slope is 0.30 and ξ = 50. Cells n = 2 to 24 were scanned and n = 4 converges to the same two-site-periodic state, so this is not a cell artifact.

The open-chain values of the right column scatter between −4.5 and 2.3 for points that the infinite chain shows to be c = 1 or gapped. At L ≤ 192 a pinned spiral or a gapped chain with ξ ≈ 50 to 300 makes the window fit meaningless, as found for point 1282 in `d2_followup_220_1282`.

As a second estimator on c = 1 points, open chains at L = 96 and 192 (χ up to 128, full end fields) were run for a random sample of 16 of the 132 c = 1 points (`data/obc_c1_sample.json`). Keeping runs with χ ≥ 64 and truncation below 10⁻⁷, the one-period-smoothed window fit gives a median of 0.96 (range 0.50 to 1.44, 15 points) at L = 96 and 1.00 (range 0.84 to 1.45, 12 points) at L = 192, and the bond-energy fit gives medians 1.04 and 1.03. The median agrees with FES, but single points scatter by up to ±0.45 even at L = 192, which is why this estimator cannot decide a single point here. The point with the largest L = 192 value, RA_18_269 (1.45), gives 0.53 at L = 96 and c = 1.00 from FES.

### Staggered-field box

On the 27 × 25 grid the ring c_S(16) is near 1 for small b at every u (the XXZ end), decreases into the gapped region as b grows, and reaches its maximum 1.27 at u = 1.0, b = 0.05, on the boundary u ≈ r = √(1 + b²/4), where the factorization line has z = 2 and C63 proved a false-positive entanglement mechanism. FES:

| (u, b) | why | cell | χ | c_fit | verdict |
|---|---|---|---|---:|---|
| (1.4, 2.0) | C63's strongest overshoot | 1 | to 128 | 1.018 | c = 1 |
| (1.3, 1.7) | near u = r | 1 | to 128 | 1.037 | c = 1 |
| (1.0, 0.05) | largest ring value, essentially on u = r | 6 | to 128 | 1.06 (median of pairs 1.01, κ = 0.43) | c = 1, slow ξ growth near the z = 2 line |
| (0.9, 0.5) | C63's clean point | 1 | to 96 | 0.985 | c = 1 |
| (0.5, 0.3) | interior | 1 | to 128 | 0.00 (ξ = 25) | gapped |
| (−0.5, 0.1) | small b, u < 0 | 1 | to 128 | 0.00 (ξ = 9) | gapped |

The two interior points are gapped where the L = 16 ring gives c_S = 0.35 and 0.03. So the gapless part of the box inside the C63/d2_gap_reduction bounds looks like a strip along u = r plus the small-b XXZ edge. I did not map the strip's inner boundary.

### What the survey does not cover

- The sampling is random. A feature narrower than about 0.15 rad on a sphere, or a junction whose two lines are not both visible as ring features at L = 12 and 16, can be missed. The patch refinement only follows junctions of the I/B/G type.
- A point reaches FES only if the ring gives c_S > 1.2 at both L = 12 and 16 (or wins a junction patch with c_S > 1.1) and the open chain then gives c > 1.2 on both measures. A c = 3/2 point with ring c_S between 1.1 and 1.2 at L ≤ 16 would be missed. For c = 3/2 the calibration argues against such a miss (cluster-Ising and Ising ⊗ XXZ give ring values 1.46 to 1.56 at these sizes in `d2_multicritical`), but an interacting c slightly above 1 could be missed.
- An open-chain value can also be too low on a spiral (point 220 gives c_dbl = 0.59). A c > 1 point whose open-chain estimate happened to fall below 1.2 would not have reached FES. The 300-point control group bounds how often this happens for c = 1 points, not for c > 1.
- FES with χ ≤ 64 resolves ξ up to a few hundred sites. A gapped point with a larger ξ, or a c = 1 point with very slow FES convergence, is classified by its trend only.
- Only one estimator decides. For a random sample of 16 of the c = 1 points, open chains at L = 96 and 192 were run as a second estimator (results above). A ring Casimir estimate with the velocity from the dispersion was not attempted, because the rotor towers of the floating phases make the velocity ambiguous (as found in `d2_classification`).

## Reproducing

numpy, scipy and numba (on Typhon, `PYTHONNOUSERSITE=1` with the default anaconda). Run from this folder. `../d2_classification/ed.py` and `../d2_followup_220_1282/{vumps,fdmrg,analyze_obc,models}.py` are imported, not copied.

- Stage 0. `python3 scan.py sphere NS 1 --N 600 --procs 16` and the 116 lines of `scan_tasks.txt` (slurm `scan_array.slurm`), about 1.4 s per point per core at L = 12 and 16. `python3 analyze_scan.py` gives the stratum table. `python3 make_patches.py` writes `patch_points.json` (880 candidates, stored here gzipped, gunzip before use), scanned by `tasks_patch.txt` (`cmd16_array.slurm`).
- Stage 1. `screen_points.json` (built from the scans as described above, with the controls first) and `obcscreen.py --points screen_points.json --start i --stop i+25 --edgefield` (`tasks_screen.txt`, `cmd1_array.slurm`, about 45 s per point). `select_patch.py` writes `patch_screen_points.json`. It was run when 42,929 of the 43,120 patch points had finished, and the one point added by the complete data is appended at the end (`tasks_pscreen.txt`, `tasks_pscreen2.txt`).
- Stage 2. `candidates.json` holds every flagged point (9-vector and origin). `python3 fes_auto.py NAME --cells ... --chis ...` (`tasks_fes1.txt` to `tasks_fes4.txt`, `tasks_stag.txt`, `tasks_anom.txt`; `cmd_array.slurm`, 4 cores, 10 min to 2 h per point). Calibration `tasks_calib.txt` with `fes.py` and `obc.py`. Long runs `tasks_long.txt` (`cmd8_array.slurm`). Then `python3 analyze.py fes` and `python3 summary.py` (writes `data/summary.json`).
- Stage 3. `tasks_anom.txt`, `tasks_anom_ed.txt` (ring ED at L = 18, 20, 22), `tasks_obc_c1.txt`, then `python3 analyze.py obc`.
- Data here. `data/fes/` (every FES run, per χ), `data/obc/` (open chains at L = 96 and 192), `data/obcscreen/` (gzipped screening output), `data/scan_compact.jsonl.gz` and `data/patch_compact.jsonl.gz` (per point c_S, ground-state momentum, S(L/2) and the four lowest scaled gaps, from `compact.py`), `data/scan_summary.json`, `data/summary.json`. The raw ED files with every sector's energies and entropy profile (175 MB) are in `/scratch/lustre/jkudlerflam/d2_grid/contrib/d2_grid/data/scan/` on Typhon. No MPS states were saved.
- Typhon usage, partition typhon (a few apollo tasks): about 1,000 core-hours, jobs 17520291 (sphere scans), 17520316 and 17521652 (calibration), 17520648 and 17521602, 17522228 (open-chain screens), 17520791 (patches), 17521353, 17521385, 17521386, 17521491, 17521516, 17521649, 17521933 (FES), 17522009 (staggered FES), 17521650, 17521651 (anomalies), 17522061 (open-chain sample), 17522197 (long FES). A few of the L = 192 sample chains and the χ = 256 long runs were still running when this was written and are not used.

## Relation to the rest of the repository

This is the numerical half of open item 5 of STATUS.md for the generic qubit residual. It extends `d2_classification`, `d2_multicritical` and `d2_two_fluid`, which hunted in structured families, to random sampling of every residual stratum, and it moves the decision from ring and open-chain entropies to infinite-chain FES for every flagged point. The picture that emerges is simple. The generic residual consists of gapped phases and c = 1 floating spirals, the time-reversal stratum is mostly gapless, and the places where Ising-like lines meet c = 1 regions are either c = 1 or gapped at the FES level. Three directions would sharpen it: (i) a denser junction search with a criterion that does not rely on ring entropies, for instance FES at small χ along the lines themselves, (ii) larger χ on the two "open" RA points and on NS_11_249, and (iii) a direct test, at a few c = 1 floating points, of whether an Ising line can be tuned into the floating phase with the remaining parameters, which is the one mechanism for c = 3/2 that the d2_multicritical heuristic does not exclude in the symmetry-free strata.
