# Staggered-charge branch A_st^o: one boson, no neutral Ising sector

## Claim

Consider the one-site translation-invariant nearest-neighbour qubit chain

    H = sum_j [ X_j X_(j+1) - Y_j Y_(j+1) + u Z_j Z_(j+1) + b Z_j ],   b > 0,  -1 < u < r = sqrt(1 + b^2/4),

the residual staggered-field branch A_st^o of C63 (codex/work_C63/RESULT.md, Section 1). It conserves the staggered charge Q = sum_j (-1)^j Z_j. C63 Section 3 shows that, with a lattice U(1) acting nontrivially on the gapless fields, 1 < c < 3/2 is impossible and c = 3/2 needs a neutral Ising stress tensor. C63 Section 7 asks whether such a neutral Ising sector can become critical while the charged oscillator keeps a finite speed.

[numerical] At every tested point of A_st^o the answer is no. Where the chain is critical, the low-energy spectrum and both entanglement estimators are those of a single compact boson with c = 1. No neutral level outside the boson tower was found in any charge sector, while the same diagnostics pick out the extra Ising sector of an Ising x XX control immediately. The largest asymptotic central-charge estimate on the branch is

    c = 1.05 +/- 0.05   (finite-entanglement scaling at (u, b) = (1.364, 2), near the BKT edge of the liquid),

against 1.50 for the Ising x XX control with the same estimator and chi range. Every other critical point gives an FES value between 0.95 and 1.04, and open-chain values at the largest size (L = 192 or 256) between 0.80 and 1.17 that are still falling wherever they exceed 1.02 (the values below 0.9 are at b = 4, where the open-chain estimator is biased low by the polarized ends). Open-chain entropy fits close to the factorization boundary u = r reach values up to 5 at L = 32. They decrease monotonically with L at every such point (for example 4.44, 3.00, 2.05, 1.61, 1.28, 1.17 for L = 32 to 256 at r - u = 0.005, b = 0.5), and the thermodynamic-limit estimator gives 1 at the same points. C63 Section 6 supplies the mechanism. These are crossovers, not central charges.

Points tested: b = 0.5, 1, 2, 4, with u from the gapped interior up to r - 0.005. There are 22 production points with open chains L = 32, 64, 96, 128 (L = 192 at 13 of them and L = 256 at 10), six ring points with L = 20 to 40, 60 cut points at L = 64, and thermodynamic-limit FES at all 22 production points.

## Status

Numerical, on a finite set of points of a two-parameter family, with the estimators calibrated on XX, XXZ, the critical transverse-field Ising chain and an additive Ising x XX chain. It is not an exhaustive scan of A_st^o. It does not exclude a narrowly tuned point between the sampled ones, or a crossover length beyond the reach of L = 256 and chi = 256. It would be upgraded by a controlled continuum argument for the near-boundary regime (where the charge stiffness grows without bound and the boson velocity vanishes), or by an independent DMRG code reproducing the tables.

## How it was established and checked

- Author: Claude Opus 5.5 subagent of the c-d conjecture campaign (Jonah Kudler-Flam), 2026-10-28. All code in this folder was written for this task. The FES runs use the damped VUMPS of `contrib/d2_followup_220_1282/vumps.py`, imported unchanged.
- No independent referee yet.
- Checks performed:
  - The U(1) DMRG reproduces exact diagonalization of the same MPO at L = 10 (open and periodic, two charge sectors, ground and first excited state, three (u, b) points) to 1e-12, and the d = 4 Ising x XX and TFI MPOs at L = 5 and 8 to 1e-12.
  - Open XX at L = 64 matches the free-fermion energy to 6e-10 (relative 7e-12).
  - Independent L = 192 runs (local, ground state only, and Typhon, with excited states) at three near-boundary points agree to 2e-11 in E0 and 2e-7 in c_obc.
  - Ring calibrations reproduce known Luttinger data. At L = 40, XX gives c = 1.002, K = 1.00, v = 3.99 (exact 1, 1, 4). XXZ at u = 0.5 (Delta = -1/2 in the rotated frame) gives c = 1.002, K = 1.504, v_c = 2.607, v_n = 2.606 (exact 1, 3/2, 3 sqrt(3)/2 = 2.598).
  - Excited states are orthogonal to the lower states of the same sector at convergence (the overlaps are stored as `overlaps` in each data file).
  - Every claim below is tagged [numerical] (with sizes and estimator) or [heuristic].

## Method

### Model and charge

On an even chain, alternate-site X rotations (a product unitary, so all entanglement and energies are unchanged) give C63 eq. (4),

    H' = sum_j [ X_j X_(j+1) + Y_j Y_(j+1) - u Z_j Z_(j+1) ] + b sum_j (-1)^j Z_j,

with n_j = (1 - Z_j)/2 conserved and Q = L - 2N. The ground state lies at N0 = L/2, and N0 + q and N0 - q are degenerate by particle-hole combined with reflection (open) or translation (ring). Open chains carry the full field b on every site, including the two ends.

### Symmetry-resolved DMRG (`udmrg.py`)

Finite two-site DMRG written in numpy/scipy for this task. Every bond index carries a U(1) charge, bond bases are sorted by charge, the SVD is done block by block, and the effective Hamiltonian is applied as a precomputed list of dense block products, so each run lives exactly in one sector N. Periodic chains use the same code with an MPO that carries the wrap-around bond across the chain. Excited states in a sector are found by projecting out the already converged MPS (penalty 20 |phi><phi|, exact overlap environments). The bond dimension is grown through chi = 64, 128, 256 (up to 384 for rings and the d = 4 control). The largest discarded weight at the final chi is 3e-11 on open chains and 1e-8 on rings.

For each point and size the run stores, in the sector q = 0, the ground state and three excited states, the entanglement S(l) and the block charge variance F(l) = <Q_A^2> - <Q_A>^2 on every bond (F is exact from the charge labels of the Schmidt states), and the charge-resolved Schmidt spectrum at the central and quarter bonds. In q = 1 it stores the ground state and one excited state, and in q = 2 the ground state. The L = 192, 256 runs keep only the q = 0 ground and first excited states and the q = 1 ground state (or only the q = 0 ground state, files `gs_*`).

### Estimators

[numerical, calibrated] Open chains, with chord length lambda = (2L/pi) sin(pi l/L), fitted separately on even and odd bonds over L/8 <= l <= 7L/8 (to remove the parity oscillation) and averaged,

    S(l) = (c/6) ln lambda + const,   F(l) = (K/(2 pi^2)) ln lambda + const.

Rings use lambda = (L/pi) sin(pi l/L), S = (c/3) ln lambda, F = (K/pi^2) ln lambda (an MPS bond of a ring is a two-endpoint interval). For one boson the charge tower is E(q) - E0 = (pi v/L) q^2/(2K) on open chains and (2 pi v/L) q^2/(4K) on rings, so

    v_c = 2 K L D1/pi,   x_n = (E_n - E0)/(pi v_c/L) on open chains, (E_n - E0)/(2 pi v_c/L) on rings,

with D1 = E0(N0 + 1) - E0(N0). A single boson predicts D2/D1 -> 4, open-chain neutral levels x_n = 1, 2, 2, 3, ..., ring neutral levels x_n = 1, 1 (momentum +/-2pi/L), then 2 and K, and the same oscillator gaps inside every charge sector (x_c). An additional neutral Ising sector with speed v_I adds levels at x = (v_I/v_c) h with h in {1/2, 3/2, 2, ...} (free boundary) in every charge sector. Since x_n1 = v_n/v_c, the same number compares the neutral and charged speeds directly.

[numerical, calibrated] Finite-entanglement scaling (`fes_point.py`) uses dense damped VUMPS on the original, one-site-invariant Hamiltonian with a two-site cell, chi = 8 to 256, and fits S = (c/6) ln xi, with xi from the transfer-matrix gap. This is a thermodynamic-limit estimator that involves neither open boundaries nor the finite-size crossover near u = r.

## Results

### Calibration

| Control | true c | OBC c_obc at L = 128 / 256 | c_inf from L = 64, 96, 128 | ring c at L = 40 | FES c (joint fit, largest four chi) | neutral spectrum x_n |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| XX (u = b = 0) | 1 | 0.980 / 0.990 | 1.002 | 1.002 | 0.998 | 1.020, 2.040, 2.040 (L = 128) |
| XXZ, u = 0.5, b = 0 | 1 | 0.986 / 0.993 | 1.002 | 1.002 | 0.989 | 1.015, 2.031, 2.031 (L = 128) |
| critical TFI | 1/2 | 0.509 / - | 0.500 | - | 0.493 | L(E1 - E0) = 3.093 at L = 32, against pi v h = pi for h = 1/2 |
| Ising x XX (d = 4) | 3/2 | 1.486 / - | 1.489 | - | 1.501 | 0.262, 0.786, 1.040 (L = 64) |

The last row is the signature this test looks for. With v_I = 2 and v_c = 4, the Ising free-boundary levels h = 1/2 and 3/2 appear at x = 0.25 and 0.75 (0.262 and 0.786 at L = 64), below the boson oscillator at x = 1.04, and they appear identically in the q = 1 sector (x_c = 0.262). Both entanglement estimators read 3/2 within 0.015.

### Production points, open chains

[numerical, open chains, L = 64 to 256, chi <= 256] The entanglement coefficient c_obc, charge stiffness K_F, and the neutral (x_n, sector q = 0) and charged-sector (x_c, sector q = 1) spectra. Complete rows, including L = 32, the charge-gap ratios and velocities, are in `summary.json`.

| (u, b) | regime | c_obc at L = 64 / 96 / 128 / 192 / 256 | K_F at Lmax | x_n at L = 128 | x_c at 128 |
| --- | --- | --- | ---: | --- | ---: |
| (0.80, 0.5) | gapped side of BKT, xi > L | 0.875 / 0.857 / 0.839 / - / - | 1.51 | 1.170, 1.684, 1.977 | 0.779 |
| (0.83, 0.5) | near BKT | 0.991 / 0.990 / 0.988 / 0.986 / 0.984 | 2.01 | 1.000, 1.637, 1.824 | 0.813 |
| (0.88, 0.5) | liquid | 1.077 / 1.062 / 1.054 / - / - | 2.70 | 0.945, 1.800, 1.857 | 0.909 |
| (0.93, 0.5) | liquid | 1.071 / 1.043 / 1.031 / 1.020 / 1.015 | 3.28 | 0.973, 1.944, 1.951 | 0.972 |
| (0.98, 0.5) | liquid, r - u = 0.05 | 1.101 / 1.054 / 1.036 / 1.022 / 1.015 | 4.72 | 0.974, 1.948, 1.961 | 0.974 |
| (1.0258, 0.5) | r - u = 0.005 | 2.995 / 2.050 / 1.612 / 1.278 / 1.169 | 17.3 | 0.739, 1.478, 1.593 | 0.739 |
| (0.95, 1) | gapped | 0.542 / 0.434 / 0.352 / - / - | 0.50 | 2.649, 3.147, 3.694 | 0.922 |
| (1.0, 1) | liquid, K near 2 | 1.021 / 1.002 / 0.996 / 0.993 / 0.993 | 2.20 | 0.977, 1.803, 1.871 | 0.890 |
| (1.068, 1) | liquid, r - u = 0.05 | 1.078 / 1.041 / 1.028 / - / - | 4.15 | 0.980, 1.959, 1.971 | 0.980 |
| (1.098, 1) | r - u = 0.02 | 1.381 / 1.168 / 1.103 / 1.059 / 1.041 | 6.89 | 0.928, 1.855, 1.886 | 0.928 |
| (1.113, 1) | r - u = 0.005 | 3.291 / 2.128 / 1.631 / 1.278 / 1.169 | 15.4 | 0.720, 1.439, 1.543 | 0.720 |
| (1.35, 2) | gapped | 0.400 / 0.253 / 0.172 / - / - | 0.19 | 6.364, 10.362, 11.425 | 3.319 |
| (1.364, 2) | near BKT | 0.918 / 0.903 / 0.902 / 0.906 / 0.910 | 1.81 | 1.126, 2.647, 3.047 | 1.473 |
| (1.38, 2) | liquid | 0.893 / 0.926 / 0.951 / - / - | 2.96 | 1.057, 2.054, 2.107 | 1.038 |
| (1.394, 2) | liquid, r - u = 0.02 | 1.112 / 1.023 / 1.008 / 1.002 / 1.000 | 4.47 | 1.008, 2.017, 2.040 | 1.008 |
| (1.4042, 2) | r - u = 0.01 | 1.588 / 1.212 / 1.106 / 1.051 / 1.034 | 6.89 | 0.937, 1.873, 1.923 | 0.937 |
| (1.4092, 2) | r - u = 0.005 | 2.625 / 1.693 / 1.353 / 1.145 / - | 10.8 | 0.719, 1.436, 1.576 (L = 96) | 0.718 (96) |
| (2.2, 4), (2.21, 4), (2.216, 4) | gapped (FES xi <= 23) | open-chain fits not meaningful, see below | - | - | - |
| (2.2261, 4) | near BKT | 0.804 / 0.741 / 0.755 / 0.796 / - | 1.91 | 1.477, 3.255, 3.367 | 1.775 (96) |
| (2.2311, 4) | r - u = 0.005 | 1.217 / 0.934 / 0.882 / 0.896 / - | 3.98 | 1.176, 2.341, 2.537 (L = 96) | 1.174 (96) |

A dash means that size was not run for that point.

[numerical] Deep in the liquid (K_F >= 3) the spectrum is the boson tower to within about 3 percent at L = 128. x_n2 and x_n3 are twice x_n1 and nearly degenerate, x_c equals x_n1 (the same oscillator in the q = 1 sector), D2/D1 = 4.00, and x_n1 = v_n/v_c rises toward one with L (0.949, 0.965, 0.973, 0.982, 0.986 at (0.93, 0.5) for L = 64 to 256). There is no level below the oscillator in either sector, and the charged and neutral speeds agree to the same accuracy.

[numerical] Near the BKT edge (K_F close to 2, at (0.83, 0.5), (1.0, 1), (1.364, 2), (2.2261, 4)) the boson pattern is distorted by amounts that barely change with size. D2/D1 lies between 3.6 and 4.6, x_c between 0.81 and 1.78, and the level-two doublet is split. The distortion tracks K_F, not u or b. It appears wherever K_F is near 2 and is gone once K_F exceeds 3. [heuristic] This is the expected effect of the marginal staggered-field vertex cos(2 phi) (dimension K) on an open chain. Its first-order matrix elements scale like L^(1-K), the same order as the level spacing at K = 2, and they decay only logarithmically. An independent Ising sector would not track K, would appear identically in the q = 0 and q = 1 sectors, and would raise both entanglement estimators by 1/2, as the control shows. None of this is seen. At the b <= 2 points c_obc sits at 0.90 to 0.99 and changes by at most 0.01 between L = 96 and 256, and FES gives 1.00 to 1.05.

[numerical] Close to u = r (r - u <= 0.02), c_obc is large at small L and falls monotonically at every point, to 1.17 at L = 256 for r - u = 0.005 (b = 0.5 and 1), to 1.15 at L = 192 for r - u = 0.005 (b = 2), and to 1.00 to 1.04 at L = 256 for r - u = 0.01 to 0.02. K_F reaches 12 to 25 and falls slowly with L, v_c drops to 0.4 or less, and x_n1 = x_c sits below one (0.6 to 0.9) while L is below the crossover length, rising with L. The boundary line is the factorization line of C63 Section 6, with L + 1 degenerate ground states, a z = 2 defect branch, and an exact half-charge entropy growing as (1/2) ln L, which an open-chain fit reads as c = 3. The near-boundary values are therefore the crossover from that line, with a crossover length that grows as u -> r. FES at the same points gives 0.95 to 1.04. The trend at every one of them is downward, and none has a c > 1 asymptote at the sizes reached.

[numerical] At b = 4 the open-chain profiles of the gapped points (u = 2.19 to 2.22) have S(l) and F(l) decreasing toward the centre, so the chord fits return small negative numbers. This is charge bound to the strongly polarized ends, not a bulk feature, and FES confirms a gap (xi = 7.6, 12.6, 22.8 at chi = 256 for u = 2.2, 2.21, 2.216). The two b = 4 critical points give open-chain values below one (0.80 and 0.90 at L = 192) and FES values 1.00 and 1.04, so the same boundary effect seems to depress the open-chain estimate there. It does not affect the c > 1 question.

### Rings

[numerical, rings, L = 20, 24, 32, 40, chi <= 384]

| (u, b) | c_ring at L = 24 / 32 / 40 | K_F at 40 | lowest ring levels x_n at L = 40 (q = 0) | x_c at 40 | v_c, v_n at 40 |
| --- | --- | ---: | --- | ---: | --- |
| (0, 0) XX | 1.007 / 1.003 / 1.002 | 1.00 | 1.002 (four-fold) | 0.999 | 3.989, 3.996 |
| (0.5, 0) XXZ | 1.009 / 1.004 / 1.002 | 1.50 | 1.000, 1.000, 1.493, 1.493 | 0.997 | 2.607, 2.606 |
| (0.93, 0.5) | 1.046 / 1.023 / 1.014 | 3.43 | 0.980, 0.980, 1.934, 1.934 | 0.978 | 1.082, 1.060 |
| (1.098, 1) | 1.149 / 1.098 / 1.067 | 7.18 | 1.107, 1.107, 2.169, 2.169 | 1.105 | 0.453, 0.501 |
| (1.38, 2) | 1.076 / 1.048 / 1.033 | 3.40 | 0.982, 0.982, 1.798, 1.937 | 0.968 | 0.596, 0.586 |
| (2.216, 4) | 0.837 / 0.789 / 0.732 | 0.71 | gapped (1.48, 1.53, 1.53, 2.47) | 0.915 | - |

[numerical] At every critical ring point the lowest neutral level is the momentum +/-2pi/L oscillator doublet at x = 1, followed by level-two states, with the same oscillator in the q = 1 sector. Charged and neutral velocities agree to 2 percent at L = 40, except at (1.098, 1), where they are converging toward each other from 0.35 and 0.69 at L = 20 to 0.45 and 0.50 at L = 40. The ring entanglement coefficient decreases toward one with L at every point.

### Finite-entanglement scaling

[numerical, VUMPS, two-site cell, chi = 8 to 256, c from S = (c/6) ln xi]

| (u, b) | chi range | xi at largest chi | c, joint fit of the largest four chi | c from consecutive chi pairs (last three) |
| --- | --- | ---: | ---: | --- |
| XX | 8 to 256 | 1164.6 | 0.998 | 1.057, 0.971, 0.994 |
| XXZ u = 0.5 | 8 to 256 | 1481.4 | 0.989 | 0.974, 1.005, 0.972 |
| Ising | 8 to 192 | 9473.3 | 0.493 | 0.471, 0.533, 0.461 |
| Ising x XX | 8 to 256 | 192.6 | 1.501 | 1.589, 1.493, 1.440 |
| (0.83, 0.5) | 8 to 256 | 1371.8 | 1.030 | 1.007, 1.048, 1.020 |
| (0.88, 0.5) | 8 to 256 | 1484.7 | 0.992 | 0.963, 1.015, 0.977 |
| (0.93, 0.5) | 8 to 256 | 1624.7 | 0.990 | 0.965, 1.017, 0.964 |
| (0.98, 0.5) | 8 to 256 | 1972.1 | 0.990 | 0.984, 0.995, 0.987 |
| (1.0258, 0.5) | 8 to 96 | 1398.2 | 0.982 | 1.026, 0.937, 0.975 |
| (1.0, 1) | 8 to 256 | 1210.8 | 1.003 | 0.997, 1.007, 1.000 |
| (1.068, 1) | 8 to 256 | 1533.5 | 0.976 | 0.988, 0.961, 0.989 |
| (1.098, 1) | 8 to 192 | 1413.2 | 0.959 | 0.961, 0.996, 0.917 |
| (1.113, 1) | 8 to 96 | 1060.2 | 0.952 | 0.984, 0.942, 0.917 |
| (1.364, 2) | 8 to 256 | 1005.2 | 1.051 | 1.036, 1.061, 1.048 |
| (1.38, 2) | 8 to 256 | 1168.8 | 0.981 | 0.996, 0.979, 0.967 |
| (1.394, 2) | 8 to 256 | 1308.3 | 0.961 | 0.999, 0.951, 0.937 |
| (1.4042, 2) | 8 to 128 | 789.1 | 0.972 | 1.089, 0.920, 0.961 |
| (1.4092, 2) | 8 to 64 | 465.7 | 1.022 | 0.971, 1.037, 1.060 |
| (2.2261, 4) | 8 to 128 | 448.1 | 1.004 | 1.032, 1.014, 0.954 |
| (2.2311, 4) | 8 to 48 | 199.6 | 1.042 | 1.124, 1.113, 0.917 |
| (0.8, 0.5), (0.95, 1), (1.35, 2), (2.2, 4), (2.21, 4), (2.216, 4) | 8 to 256 | 670, 73, 82, 7.6, 12.6, 23 | gapped, xi saturates and S stops growing | - |

The near-boundary FES runs stop at smaller chi because VUMPS converges slowly there (the runs to chi = 256 were still going when this file was written; they write into `data/fes/` as each chi finishes). The pair-to-pair scatter of the calibrated XX and Ising x XX runs, about +/-0.05 and +/-0.07, sets the error bar. All critical points of A_st^o lie between 0.95 and 1.05 and the Ising x XX control at 1.50. The maximum on the branch, 1.051 at (1.364, 2), is within that scatter of one and 0.45 below the control.

### Cuts at fixed b, L = 64

[numerical, open chains, L = 64, chi = 128, 60 points, rows in `summary_cuts.json`] At b = 0.5 (u = 0.74 to 1.02), b = 1 (0.90 to 1.11), b = 2 (1.32 to 1.40) and b = 4 (2.19 to 2.23), and at u = -0.9, -0.5, -0.2, 0.2 for b = 0.5, 1, 2:

- u <= 0.2 is gapped at every b (c_obc below 0.01, charge gaps D1 between 0.7 and 6 at L = 64).
- Along each cut c_obc, K_F and x_n1 vary smoothly. The sequence is gapped, a BKT-like onset with K_F passing through 2, a liquid with x_n1 between 0.92 and 1.02 and x_n2 close to 2 x_n1, and finally the boundary crossover with K_F > 6. There is no interior point where a second neutral level comes down, where x_n1 dips and recovers, or where c_obc rises above its neighbours by more than the smooth crossover trend. The one non-monotone feature, c_obc = 0.83 at (1.3736, 2) between 0.90 and 0.91, is a dip, not a peak.

An exploratory L = 32 and 64 scan (`data/scan/`, `summary_scan.json`, 32 points at r - u = 0.005 to 0.7) located the liquid strips used above and shows the same pattern.

## Conclusion

[numerical] On the staggered-charge branch, at b = 0.5 to 4 and from the gapped interior to within 0.005 of the factorization line, the charged boson is the only critical sector. Its velocity stays finite away from u = r and vanishes only on approach to the line, which is also where the open-chain entropy coefficients become large, while the neutral spectrum shows no independent level. The neutral Ising scenario of C63 Section 7 is not realized at any tested point. The largest asymptotic c is 1.05 +/- 0.05 (FES), and every finite-size value above 1.02 falls monotonically with size.

[heuristic] Two questions remain open: (i) the asymptotic form of the near-boundary crossover, which a continuum treatment of the dilute z = 2 defect gas near u = r could settle; (ii) the strip at b > 4, where the liquid becomes very narrow and the open-chain estimator is affected by the strongly polarized ends.

## Reproducing

Python 3 with numpy and scipy only. The FES driver also imports `contrib/d2_followup_220_1282/vumps.py`. Run from this folder.

```sh
# one open-chain point (sectors q = 0, 1, 2; three neutral excited states; one in q = 1)
OPENBLAS_NUM_THREADS=4 python3 run_point.py stag --u 0.93 --b 0.5 --L 128 --chis 64 128 256
# ring
OPENBLAS_NUM_THREADS=4 python3 run_point.py stag --u 0.93 --b 0.5 --L 40 --pbc --chis 64 128 256 384 --nexc 4
# controls
python3 run_point.py tfi --L 128 --chis 32 64 128 --nexc 3
python3 run_point.py isxx --L 128 --chis 64 128 256 384
# thermodynamic-limit FES
python3 fes_point.py stag --u 0.93 --b 0.5 --chis 8 12 16 24 32 48 64 96 128 192 256
python3 fes_point.py isxx --chis 8 12 16 24 32 48 64 96 128 192 256
# all production runs (SLURM): one line of tasks_main.txt or tasks_big.txt per array task
sbatch --array=0-153 job_array.slurm
sbatch --export=ALL,TASKS=tasks_big.txt --array=0-31 job_array.slurm
# fits and tables
python3 analyze.py                                   # summary.json (open chains and FES)
python3 analyze.py --sub pbc --json summary_rings.json
python3 analyze.py --pattern 'cut*' --json summary_cuts.json
python3 analyze.py --sub scan --json summary_scan.json
python3 tables.py                                    # markdown tables, including the L-extrapolation
```

The L = 64 cuts used `run_point.py stag --u U --b B --L 64 --chis 64 128 --nexc 2 --nexc1 1 --out cut_uU_bB_L64`. The ground-state-only L = 192, 256 runs used `--chis 64 128 256 --qs 0 --nexc 0 --out gs_uU_bB_LL`. The exploratory scan used `--chis 48 96 --nexc 1 --nexc1 0`. Recorded wall times on four cores: 20 to 105 minutes for an L = 128 open-chain point with all seven states, 20 to 75 minutes for an L = 40 ring point, and up to two hours for an FES point.

| File | Content |
| --- | --- |
| `udmrg.py` | U(1)-resolved finite DMRG, MPOs for the rotated staggered chain, TFI and Ising x XX |
| `run_point.py` | charge-resolved ground and excited states, S, F and Schmidt spectra for one point |
| `fes_point.py` | VUMPS finite-entanglement scaling in the original frame |
| `analyze.py`, `tables.py` | fits, extrapolation and the tables of this file |
| `data/obc`, `data/pbc`, `data/fes`, `data/scan` | raw results (JSON, no MPS tensors) |
| `summary*.json` | fitted quantities per point and size |
| `tasks_main.txt`, `tasks_big.txt`, `job_array.slurm` | Typhon job definitions |

## Relation to the rest of the repository

This is the focused test named in C63 Section 7 (codex/work_C63/RESULT.md) for the A_st^o row of the C63 covering table. It moves that row from "small ED, globally unresolved" to "DMRG to L = 256 and FES at 22 points, one boson everywhere tested", and does not touch the unitary-Z2, antiunitary-only or symmetry-free rows. It complements, and does not duplicate, the grid survey in `contrib/d2_grid/`.
