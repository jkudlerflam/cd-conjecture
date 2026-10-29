# Dressed-defect costs in the vacuum-clock vacancy model by DMRG

## Claim

[Numerical] In the open-chain version of the legally projected vacuum-clock model of C89 and C94 (cells of selected length b = 15, vacancies hopping with t = 1, clock and marker energies entering as J times the contact projector), a pinned cell costs a positive energy of about 0.86 J at material density 0.6 and L = 50, and the second-order residual for a pin has the favourable sign (it raises the cost, as C89's operator-valued calculation found). This is converged to about one percent in bond dimension at J/t = 0.04 and 0.01.

[Numerical, unresolved] For wrong-length cells the test asked for did not converge in the time available. The first-order cost of a (14, 16) pair is very small in this geometry (about 1.9e-5 J for W = 30), and the DMRG energy differences at bond dimension 64 have errors of 2e-7 to 5e-7, comparable to the cost itself. The sign of the pair cost and the value of ε(J/t) for wrong-length cells are therefore not determined here. The negative entries in the raw tables are within that error and should not be read as negative costs.

## Status

Partial. The code is validated against exact diagonalization and runs at L = 100, but at the snapshot taken for this report only the pin comparison at L = 50 is converged in bond dimension. Wrong-length-cell comparisons at L = 50 and every comparison at L = 100 exist only at bond dimension 32 or 64, where the DMRG error (up to 1e-3 at L = 100, χ = 32) is larger than the residual ΔR being sought. The runs were left going on Typhon (directory `defect_cost_dmrg` in scratch) with χ up to 128, and `analyze.py` will regenerate the tables from their output.

## How established and checked

- `check_ed.py`: DMRG against an independent sparse exact diagonalization built in C89's construction (gap compositions times clock product, open chain) for nine configurations including wrong-length cells, a merged cell and pins, J/t = 0.04, 0.1, 0.3, dimensions up to 14,784. Energy agreement within 2.2e-11, free (J = 0) energies agree with the closed form below within 1e-13, and the J = 0 contact profiles agree within 1e-13. All pass (`check_ed.json`).
- The J = 0 energy is exact: with material at both ends, the free problem is V free fermions on the L − 2 interior sites. At L = 50 the χ = 128 DMRG free energy differs from it by 1.7e-6, which is the gap cutoff G = 5 (below).
- `check_gapcut.py`: the cutoff g_j ≤ G on the gap occupations. At N = 12, V = 8 the free energy error is 1e-6 at G = 5 and 1e-9 at G = 6 (a Gaussian tail, as expected for emptiness formation). Its effect on defect costs (b = 6 analogue, J = 0.01) is 4.6e-9 for a wrong-length pair and 1.5e-8 for a pin, so G = 5 is not the limiting error.
- Bond-dimension convergence is judged by the change in the defect cost between the last two χ, reported as `chi_chg`.

## Method and results

Model (`model.py`). Material sites j = 0..N−1 in order, site j carrying the gap g_j (vacancies between material j and j+1, unbounded with norm-one shifts, truncated at G) and a clock bit. The last site carries no gap, so the open chain has material at both ends. Terms, with C89/C94's signs and coefficients: vacancy hopping −t between adjacent internal gaps, clock hopping −J A_j between neighbouring sites of one non-pinned cell, and the marker energy J μ/2 (A_{m−1} + A_m) at each cell's last site, with μ = selector(15) and A_j = [g_j = 0]. A pinned cell has no clock. The one difference from `clock_model` is that the ring is replaced by an open block (A_{N−1} is absent).

DMRG (`qdmrg_sd.py`). The block-sparse two-site code of `contrib/vacancy_n2_dmrg/qdmrg.py`, extended to one bond matrix per bond and a filter on the bond charges kept in each SVD. Charges are (vacancies, clocks). The filter imposes exactly one clock per non-pinned cell and none in a pinned cell (the cumulative clock count at each cell boundary is fixed), and puts no vacancy on the last site.

Residual. Following C94, R = E − E_0 − E_1, with E_1 the first-order (degenerate perturbation theory) clock energy computed from the exact J = 0 contact profile ⟨A_j⟩_0. On the ring this reduces to C94's J α Σ cell_energy, since there ⟨A_j⟩_0 = α. Then ΔR = R_σ − R_good, W = wrong-cell length plus pin count, and ε = max(0, −ΔR)/(J W).

Runs: N = 30 (two cells, V = 20, L = 50) and N = 60 (four cells, V = 40, L = 100), J/t = 0.04, 0.02, 0.01, 0.005, 0.0025, defects: a (14, 16) pair, a merged 30 cell, one pin, and at N = 60 two pairs and two pins. Vacancy scan at N = 30, J = 0.01, V = 18 to 22.

Converged results (L = 50, χ = 64 against 32 for the change):

| defect | J/t | W | cost | χ change | cost/(J W) | ΔR | ΔR/(J² W) |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| pin | 0.04 | 1 | 3.450e-2 | 8e-5 | 0.862 | +1.6e-4 | +0.10 |
| pin | 0.01 | 1 | 8.592e-3 | 4e-5 | 0.859 | +8.9e-6 | +0.09 |
| pin, V = 21 | 0.01 | 1 | 8.266e-3 | 6e-5 | 0.827 | +2.6e-5 | +0.26 |

The pin cost per J is about 0.86 at material density 0.6 and drifts down slowly with V (0.93 J at V = 18 to 0.79 J at V = 22, χ = 32). That is larger than the 0.583 J ring reference at half density, consistent with the higher contact probability at higher material density. The positive ΔR for pins at χ = 64 agrees with C89. At χ = 32 the pin ΔR is negative (about −3e-5), so the χ = 32 numbers for ΔR are not trustworthy even for pins.

Unresolved (L = 50, (14, 16) pair, W = 30): cost −5.8e-7 at J = 0.04 (χ change 5.3e-7) and +1.2e-7 at J = 0.01 (χ change 1.9e-7). The first-order cost is J × 1.85e-5, of the same size as the DMRG error at χ = 64. In this geometry both cells touch the open ends, and the Friedel profile of ⟨A_j⟩_0 reduces the first-order pair cost well below the ring value α × 9.0e-4 J, which makes this a hard test at this size. The merged 30 cell has a positive cost of order 1e-3 at χ = 32 at every J, but it is not χ-converged (it is not monotone in J).

Vacancy crossings (N = 30, J = 0.01, χ = 32). At zero temperature the grand-canonical cost min_V[E_σ(V) + sV] − min_V[E_good(V) + sV] is at least min_V of the canonical costs, so positivity at crossings follows from canonical positivity at every V. The pin costs are positive at all five V (7.9e-3 to 9.3e-3) and so at every crossing tilt. For the pair the canonical costs are again below resolution. `crossing.json` records the scan.

L = 100. All runs were still at χ ≤ 64 at the snapshot. The χ = 32 data have DMRG errors near 1e-3, far above ΔR, and are kept only as raw output.

What this does and does not show. The pin defect cost stays positive and of order J at L = 50, and its residual has the favourable sign, at the two converged J. Nothing here bears on ε(J/t) → 0 for wrong-length cells. Reaching it needs energy differences good to about 1e-9 at L ≥ 100, which means χ of 128 or more with many more sweeps than this budget allowed, or a geometry with the defect far from the open ends (heuristic). A cheaper route at these sizes would be to compute ΔR at second order directly from the DMRG free state (a contact-contact response contracted with the clock sources), which avoids subtracting two large energies.

## Reproducing

From this folder: `python3 -B check_ed.py` (a few minutes), `python3 -B check_gapcut.py`, then `python3 -B task.py <line>` for each line of `tasks.txt` and `tasks_cross.txt` (the runs used 46 parallel single-thread processes on the Typhon login node at nice 19, `xargs -P 46 -L 1 python3 -B task.py < tasks.txt`). Each task writes `out/<tag>.json` after every bond dimension. `python3 -B analyze.py` writes `results.json`, `crossing.json` and the tables (`analysis_stdout.txt` is the snapshot used here). Files: `model.py` (Hamiltonian, exact free energy, open-chain ED), `qdmrg_sd.py` (DMRG), `dmrg_run.py` (driver).
