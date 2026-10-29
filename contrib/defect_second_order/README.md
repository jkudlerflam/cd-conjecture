# Second-order defect residual on rings from the exact free tagged kernel

## Claim

[Numerical] In the zero-port vacuum-clock model of C89/C94 on a ring (selected b = 15, t = 1, isolated cell clock ground states), the second-order coefficient ΔB = B_defect − B_good, with ΔR = −J² ΔB + O(J³/t²), has been computed exactly (to about 1e-9) for rings up to N = 90 material sites (L = 180 at material density 1/2, L = 150 at density 0.6). ΔB is independent of ring size to four or five digits from N = 30 to N = 90, and for pairs it is independent of the separation of the two wrong cells.

[Numerical] Wrong-length pairs (14,16), (13,17), (12,18), (10,20) and merged long cells (30, 45) have ΔB < 0. Their second-order residual has the favourable sign, so ε(J/t) = max(0, J ΔB)/W vanishes at second order and the cost is raised. For the (14,16) pair, cost = J a ΔE_material [1 + 0.847 J/t] at density 1/2 and [1 + 0.837 J/t] at density 0.6, with t ΔB/W = −7.55e-6 and −1.12e-5.

[Numerical] Splitting a good cell into two short cells (7,8), (5,10), (3,12), (2,13) gives ΔB > 0, the unfavourable sign. Here ε(J/t) = (J/t) t ΔB/W, with t ΔB/W = 0.0217 for (7,8) at density 1/2 and 0.0276 at density 0.6. Relative to the first-order cost, cost = J a ΔE_material [1 − 8.85 J/t] (density 1/2) and [1 − 7.49 J/t] (density 0.6). The second-order cost therefore stays positive for J/t below about 0.11 at both densities. The relative coefficient is smaller for more unequal splits (−0.31 for (2,13)).

[Heuristic] The truncation at second order is not controlled uniformly in cell length (C94 §6). At the one size where it can be compared with ED (N = 30, V = 1) the measured residual for long cells exceeds the second-order value by 20 to 40 percent at J = 0.01 to 0.005, so the coefficients above are the J → 0 slopes, not a bound at finite J.

## Status

The exact kernel and the second-order coefficients are done for N = 30, 60, 90 at densities 1/2 and 0.6, plus the C89 check point N = 30, V = 1. Runs at N = 120 (L = 240 and 200) were still going when this was written. `python3 -B analyze.py` folds them into `results.json` once `out/N120_*.json` exists. No statement about orders beyond J² is made.

## How established and checked

- `check_kernel_ed.py`: the determinant kernel against C94's sparse-resolvent kernel (`codex/work_C94/inherited_second_order.py`, conjugate gradient in the full gap basis) at eight sizes up to (N, V) = (30, 1) and (10, 6). Largest difference 2e-10, which is the conjugate-gradient tolerance (`check_kernel_ed.json`).
- `check_B.py`: the full operator-valued coefficients B_good and B_defect against C89's `second_order_results.json` for all eight of its patterns, including (14,16), (10,20) and the merged 30 cell at N = 30, V = 1. Agreement within 2.6e-9. C89 in turn compared −ΔB with the ED residual ΔR/J² at J = 0.04 to 0.005 (for example (14,16): predicted 3.11e-4, measured 3.84e-4, 3.49e-4, 3.31e-4 at J = 0.02, 0.01, 0.005), so the chain back to exact diagonalization is closed at small size.
- Sum rule: Σ_m C_m(0) = ⟨A_0⟩ = a holds to 2e-13 at every size, and imaginary parts are below 1.3e-13.
- Quadrature: trapezoid in u = log τ, step 0.25, τ from e^-20 to 12 L/t. Changing the step and range moves the kernel by 2e-10. Pair ΔB agrees between N = 30, 60, 90 to about 1e-7 absolute, which bounds the combined finite-size and numerical error.
- The first-order ring value a ΔE_material for (14,16) at density 0.6 is 4.02e-4, that is a × 9.0e-4 with a = 0.447, matching the ring reference quoted in `contrib/defect_cost_dmrg`.

## Method and results

Second order (C94 eq. (13)). B_σ = Σ_ij χ_N(j − i) ⟨h_i h_j⟩_φ, with χ_N(r) the tagged contact resolvent kernel of the free bath and h_i the clock source at bond i (marker halves μ/2 and the clock hop). For a product of cell ground states, Σ_i ⟨h_i h_{i+r}⟩ splits into the mean sources and the within-cell covariance (`sources.py`).

Kernel (`kernel.py`). The material label shift commutes with H. In its sector m the bath is N free fermions on the L-ring with momenta 2π(n + ξ + m/N)/L and total momentum zero (C96). In the frame anchored at material 0, A_0 = n_0 n_1, and n_0 n_1 |FS⟩ is again a Slater determinant (δ_0, δ_1 and the Fermi-sea orbitals vanishing at sites 0 and 1), with norm² ⟨n_0 n_1⟩. So the sector correlator C_m(τ) = ⟨A_0 P_m e^{−τ(H−E_0)} P_m A_0⟩ is (L/N)² (1/L) Σ_s of single N × N determinants det(Ψ^† g_s Ψ), with g_s the twisted single-particle propagator translated by s. The s sum (the zero-momentum projection) is done by FFT, and each determinant is stabilized by factoring out the N lowest modes of grid m. Then χ_N(r) = Σ_m e^{iθ_m r} ∫_0^∞ dτ [C_m(τ) − δ_{m0} a²]. Cost is about N^3 L per τ point and sector, polynomial rather than the exponential sum over C96 form factors.

Results (t = 1, a ΔE = first-order cost per J, rel = −t ΔB/(a ΔE), so cost = J a ΔE [1 + rel J/t + O(J²)]). N = 90 rings; the N = 30 and 60 values agree to the digits shown except the pin and the (2,13) split (t ΔB = 0.0842 at N = 30). Full table in `analysis_stdout.txt`.

| density | defect | W | a ΔE | t ΔB | t ΔB/W | rel |
| --: | :-- | --: | --: | --: | --: | --: |
| 0.5 | (14,16), any separation | 30 | 2.675e-4 | −2.266e-4 | −7.55e-6 | +0.847 |
| 0.5 | (13,17) | 30 | 1.091e-3 | −9.31e-4 | −3.10e-5 | +0.854 |
| 0.5 | (10,20) | 30 | 7.877e-3 | −7.11e-3 | −2.37e-4 | +0.903 |
| 0.5 | long 30 | 30 | 1.315e-2 | −0.3735 | −1.245e-2 | +28.4 |
| 0.5 | long 45 | 45 | 3.300e-2 | −0.751 | −1.669e-2 | +22.8 |
| 0.5 | split (7,8) | 15 | 3.676e-2 | +0.3252 | +2.168e-2 | −8.85 |
| 0.5 | split (2,13) | 15 | 0.2679 | +0.0852 | +5.68e-3 | −0.318 |
| 0.5 | pin | 1 | 0.5833 | −0.199 | −0.199 | +0.341 |
| 0.6 | (14,16), any separation | 30 | 4.024e-4 | −3.367e-4 | −1.122e-5 | +0.837 |
| 0.6 | long 30 | 30 | 1.977e-2 | −0.4850 | −1.617e-2 | +24.5 |
| 0.6 | split (7,8) | 15 | 5.529e-2 | +0.4143 | +2.762e-2 | −7.49 |
| 0.6 | pin | 1 | 0.8773 | −0.259 | −0.259 | +0.295 |

The pin ΔB drifts slowly with N (−0.206, −0.201, −0.199 at N = 30, 60, 90, density 1/2), the only entry not converged to three digits. Its sign agrees with the open-chain DMRG and with C89.

Dependence on L and W. At these sizes ΔB is already at its large-L value, so ε at second order is a function of the defect only. For pairs ΔB/W grows roughly as (ℓ − 15)², tracking the first-order cost, so rel is nearly constant (0.84 to 0.90). For split cells ΔB is of order 0.1 to 0.4 per extra cell boundary, and ε/(J/t) ranges from 0.006 to 0.028.

Contrast with high occupancy. At N = 30, V = 1 (density 0.97, the C94 test) the signs flip: the long 30 cell is unfavourable (t ΔB/W = 3.0e-4, as in C94) and the (7,8) split is favourable. The sign of the second-order residual for length defects depends on density, not only on the defect.

## Reproducing

From this folder, numpy and scipy only.

- `python3 -B check_kernel_ed.py` and `python3 -B check_B.py` (seconds; they import C94's resolvent code read-only from `codex/work_C94`, and C89's JSON).
- `OMP_NUM_THREADS=1 python3 -B run.py N V 10` for each line of `sizes.txt` (N = 90 takes about 70 s on 10 processes, N = 120 a few minutes). Writes `out/N<N>_V<V>.json` with the kernel χ_N(r) and the defect table.
- `python3 -B analyze.py > analysis_stdout.txt` writes `results.json` and the table.
