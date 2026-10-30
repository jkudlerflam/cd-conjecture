# C123 checks

2026-10-29.

[Numerical reproduction] Run these commands from work_C123. NumPy and SciPy are the only numerical dependencies.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 ed_checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 transfer_checks.py

[Work record] The completed runs took approximately 10.0 seconds, 57.1 seconds, and less than one second. Output is xx_results.json, ed_results.json and transfer_results.json. Progress output from the first two runs is in xx_run.txt and ed_run.txt. No computation was submitted remotely.

## 1. XX finite-grid determinants

[Numerical] Fifty-one calibrations compare every individual positive history against dense full-Hilbert-space multiplication at n=4,6,8, beta=0.2,2,8, N=1,2,3, with either site 0 or both sites 0,n/2 measured. The n=8,N=3 two-site cases are omitted to keep the direct dense calibration small. The maximum relative error, normalized by the largest history weight in its case, is 1.95763e-14. Both n modulo four classes are tested. No parity sector is dropped.

[Numerical] There are 147 survey cases at L=4,8,32,128,512,2048,8192, a=0.125,0.25,0.5,1,2,4,16, and N=1,2,4. Nine further cases use L=32,512,8192, a=0.25,1,4, N=6. Thus 156 cases compute the actual numerator and denominator. Doubled rings reach length 16384. The largest history determinant is order 24, or order 26 with the zero-mode border.

[Proved implementation] Every history uses 2N projector appearances on the reflected contour. Doubled-ring histories use 4N. The physical parity projection is performed before squaring a one-ring history, retaining cross terms between the two physical sectors. The zero-mode contribution is evaluated by a bordered determinant, with no regulator. Bulk ground energies are removed analytically, and the stable exact doubling defect is restored at the end.

[Numerical] All cases satisfy the one-seam and two-seam pinching inequalities and the sandwich R_L-3N log 2<=G<=R_L+2N log 2. The large-ring G/N values are listed in RESULT.md. The scaled residual L^2(G-R_L+N log 2) at a=1,N=4 is

| L | Scaled residual |
| ---: | ---: |
| 4 | -0.753286120 |
| 8 | -0.688680976 |
| 32 | -0.667311726 |
| 128 | -0.666235502 |
| 512 | -0.666168925 |
| 2048 | -0.666164769 |
| 8192 | -0.666164488 |

[Proved and numerical distinction] This sequence checks the L^(-2) theorem at one spacing and insertion count. It does not establish uniformity as N grows with L. The uniform XX cost theorem instead follows from the inherited partition-ratio estimate.

## 2. Full-spectrum interacting ED

[Assumed conventions] XXZ uses h=Sx Sx+Sy Sy+Delta Sz Sz, d=2, M=1/2 for all listed anisotropies. Heisenberg is Delta=1. SU(3) uses the swap bond h=P, d=3, M=1. Each pinching is in the standard spin or color basis.

[Numerical] Models are XX, XXZ at Delta=-0.9,-0.5,0.5, Heisenberg, and SU(3). Qubit ring pairs are L,2L=4,8 and 6,12. SU(3) pairs are 3,6 and 4,8. For each pair use a=0.25,1,2,4 and N=1,2,4,8, for a total of 192 cases. At each case beta=N a L/M. Every one-ring history is enumerated, including all 3^8 SU(3) histories. The doubled-ring denominator is propagated directly as a density operator in each conserved block.

[Proved implementation] Qubit blocks conserve magnetization. SU(3) blocks conserve all color counts. Spin flip and color permutations identify equivalent blocks for the summed channel. For a fixed recorded history the omitted blocks are restored with their actual transformed color labels before summing. Individual sector contributions are added before squaring a history weight. Multiplying squared sector weights by a degeneracy would be wrong and is not done. The sum of retained block dimensions with multiplicity is asserted to equal d^n. Every eigenvalue in each block is retained.

[Numerical] The largest block dimensions are 924 for 12 qubit sites and 560 for 8 SU(3) sites. The largest direct dense XX calibration has dimension 256. No low-energy truncation, selected-sector approximation, or eigenvalue floor is used. These modest calculations do not implement the larger cluster-scale ED requested as a possible extension.

[Numerical] The sum of explicit one-ring histories agrees with independent channel propagation with maximum relative error 3.14485e-14. All pinching and R/G inequalities pass. The 32 XX cases agree with the independent Gaussian computation to maximum absolute error 6.75016e-13 in G. The starting C121 values at a=1,N=1 are reproduced as 0.762550822 for XX L=4, 1.888215629 for Heisenberg L=4, and 2.292662706 for SU(3) L=3.

[Numerical] At N=8, the spacing dependence on the largest vacuum-compatible rings computed is

| Model and L | a=0.25 | a=1 | a=2 | a=4 |
| :--- | ---: | ---: | ---: | ---: |
| XXZ, Delta=-0.9, L=6 | -0.015324 | -0.303410 | -0.084516 | 0.581886 |
| XXZ, Delta=-0.5, L=6 | -0.295413 | 0.341737 | 1.389574 | 3.473920 |
| XXZ, Delta=0.5, L=6 | -0.108114 | 1.446196 | 3.595692 | 7.887054 |
| Heisenberg, L=6 | 0.008222 | 1.889414 | 4.520688 | 9.758379 |
| SU(3), L=3 | 0.024716 | 2.979480 | 7.223387 | 15.596937 |

[Numerical interpretation] These are G/N, not G. Near the ferromagnetic XXZ point the modest sizes and slower velocity produce stronger temperature dependence. No monotonicity in a is assumed. SU(3) L=4 data are saved but are outside the vacuum-compatible subsequence. Only one vacuum-compatible SU(3) pair was computed, so no SU(3) vacuum-sequence size trend is established.

## 3. Temporal-transfer identities

[Numerical] transfer_checks.py builds both Hilbert-Schmidt transfer matrices explicitly for a two-site XX numerator and four-site XX denominator. The short numerator ring has its two periodic bonds as defined by the literal periodic sum. This test concerns the finite-dimensional identity, not the L>=4 XX theorem. Three values tau=0.2,1,3 and N=1,2,3,4 give 12 checks. The spectral moment formulas agree with independently enumerated histories to maximum relative error 3.24275e-14.

[Numerical] Both compressed transfer matrices have nonnegative spectra to numerical precision. The top-eigenspace overlap is positive, and every case satisfies the spectral upper bound in PROOF (17). The test retains the complete top eigenspace when degeneracies are present.

## 4. Limitations and work record

[Work record] The first ED survey used N<=4. It passed, then was extended to N=8 to examine the temporal growth rate. The final ED file contains the complete 192-case run. No failed mathematical assertion was discarded. This is floating-point evaluation of exact finite-dimensional formulas, not interval arithmetic.

[Unperformed] No larger interacting ED, SU(3) length-6 versus length-12 comparison, paired-fermion multi-pinch implementation, cluster submission, or independent agent review was performed. The no-network-command instruction precludes a Typhon submission in this session. No general interacting spectral or overlap estimate was verified beyond these finite examples.
