# C116 numerical checks

2026-10-29. All reported computations completed. Larger cluster checks were not run.

## Method and normalization

[Numerical] `check_rings.py` diagonalizes every conserved-charge block and concatenates the complete spectra. It checks the full Hilbert dimension, Hermiticity, the trace and the second spectral moment. No low-energy truncation is used. The local run includes XX and Heisenberg lengths 4, 6, 8 and SU(3) lengths 3, 6. The largest charge block has dimension 90, and the largest full Hilbert space has dimension 729. The run took about 0.07 seconds on this machine.

[Assumed calibration data] The standard exact infrared and bulk values in these conventions are

| Model | Bond | d | M | e_inf | A |
| :--- | :--- | ---: | ---: | :--- | :--- |
| XX | Sx Sx + Sy Sy | 2 | 1/2 | -1/pi | 1 |
| Heisenberg | S dot S | 2 | 1/2 | 1/4 - log 2 | pi/2 |
| SU(3) | permutation P | 3 | 1 | 1 - log 3 - pi/(3 sqrt(3)) | 4 pi/3 |

These infrared normalizations are inherited from the campaign. The ground-energy bulk constants are calibration inputs, not inferred from the small rings. The principal ratio R does not use e_inf.

[Numerical] The tested full traces include beta/L = 1/4, 1/2, 1, 2, 4, 8. Every check of |R| <= 4 beta M passed. The following values at beta/L = 4 violate R <= log(d^2):

| Model | Compared lengths | R | log(d^2) |
| :--- | :--- | ---: | ---: |
| XX | 4 and 8 | 3.370382150 | 1.386294361 |
| Heisenberg | 4 and 8 | 5.581805924 | 1.386294361 |
| SU(3) | 3 and 6 | 16.733383896 | 2.197224577 |

[Proved interpretation] These violations reject a temperature-independent positive spatial transfer rank d^2. They are not violations of a bound on A. The positivity-free analytic contradiction is in PROOF.md equation (17).

## Finite-size and thermal effects

[Numerical] At beta/L = 8, the estimator A_R = 4 R/(pi tau) gives

| Model | Smaller length L | A_R | Exact infrared A |
| :--- | ---: | ---: | ---: |
| XX | 4 | 1.095974415 | 1 |
| Heisenberg | 4 | 1.776966651 | 1.570796327 |
| SU(3) | 3 | 5.326401778 | 4.188790205 |

These small rings have substantial finite-size corrections. No extrapolation for Heisenberg or SU(3) is claimed.

[Numerical] As a separate ground-energy diagnostic, 6L(e_inf L - E_L)/pi gives

| Model | L | Diagnostic |
| :--- | ---: | ---: |
| XX | 4 | 1.076962164 |
| XX | 6 | 1.032936138 |
| XX | 8 | 1.018288711 |
| Heisenberg | 4 | 1.737294191 |
| Heisenberg | 6 | 1.648887214 |
| Heisenberg | 8 | 1.618276732 |
| SU(3) | 3 | 5.101408624 |
| SU(3) | 6 | 4.426429163 |

[Numerical] The full XX trace is also evaluated in O(L) work by summing all fermion occupations through even/odd parity recursions on both Jordan-Wigner boundary grids. This agrees with direct diagonalization to maximum absolute log-trace error 5.6844e-14 on the tested small rings. For the comparison of lengths 256 and 512 at tau = 8,

    A_R = 0.9994308758,
    A_Q = 2 Q/(pi tau) = 1.0000236273,
    target A = 1.

The different thermal corrections in R and Q explain their different approach. No approximation of the trace by its ground state is made.

## Multiple velocities, massive spectator and nonsymmetric bond

[Numerical] Two rescaled XX sectors, as in PROOF.md equation (9), have d = 4, M = 1/2, and A = 1 at every epsilon. Full-trace product evaluations at L = 2048 and tau = 1 give

| epsilon | B | R | Q |
| ---: | ---: | ---: | ---: |
| 0.5 | 4.5 | 0.365729852 | 7.644147104 |
| 0.1 | 12.1 | 0.557769428 | 19.935838722 |
| 0.02 | 52.02 | 0.616902087 | 82.827207025 |
| 0.005 | 202.005 | 0.628662231 | 323.332826298 |

The fast growth of Q is the reciprocal-velocity effect proved in PROOF.md. These are elementary O(L) sums, including the paired length 4096, not large matrix computations. No fit in epsilon is used.

[Numerical] The exact massive Ising expression (21) was checked for J = 0.5, 0.05, 0.005 and L = 4, 16, 64, 256, 1024 at tau = 4. Every value lies between zero and log 2. The zero contribution to the extracted A follows analytically for every J > 0, not just from these samples.

[Numerical] `check_comparisons.py` also tests an oriented complex qubit bond with no imposed symmetry on rings of lengths 3 and 6. It verifies the gluing and Golden-Thompson bounds, scalar-shift invariance of R, and exact cancellation of a coboundary. Here the test uses the representative's half-spread, an upper bound on M; it does not numerically optimize the gauge. The proof of the optimized estimate is independent of this test. All checks passed.

[Numerical] For the twisted XX bond theta = pi/4, full diagonalization gives E_4 = -1, E_8 = -2.613125929753 and R_4(16) = -8.466897057. This verifies the analytic even-ring sign counterexample in PROOF.md. The largest matrix in this additional check has dimension 256.

[Numerical] Direct quadrature of the XX infinite-chain pressure,

    p(beta) = integral_0^(2 pi) log(1+exp(-beta cos k)) dk/(2 pi),

gives Xi_4(16) = 2.126715160 and Xi_4(32) = 4.445647934, both above log 4. The analytic bound at beta = 32 already proves the contradiction, so no conclusion depends on quadrature error estimates.

## Reproduction and unavailable cluster run

[Work record] Run the small checks inside this folder:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 check_rings.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 check_comparisons.py

The results are `checks_small.json` and `comparisons.json`. Larger ED sizes are enabled only by `--cluster`, and the script refuses that option without SLURM_JOB_ID. The prepared `run_typhon.slurm` requests one CPU, 4 GB and 20 minutes, for XX/Heisenberg through L = 12 and SU(3) through L = 9. No GPU or computation on the login node is needed.

[Blocked computation] Two Typhon connection attempts, separated by analytic work, failed with `Connection timed out during banner exchange`. No SLURM job was submitted, no remote file was written, and no larger diagonalization was substituted locally. The cluster-size checks remain unperformed. Both completed local scripts run in less than a second each.
