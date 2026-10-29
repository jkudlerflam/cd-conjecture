# C118 numerical checks

2026-10-29.

[Numerical] All completed checks use full traces. The interacting spectra contain every conserved-color or magnetization sector, with Hermiticity, Hilbert dimension, trace and second-moment checks inherited from the inspected C116 spectrum routine. Local rings are XX and Heisenberg L = 4,6,8, and SU(3) L = 3,6. The largest charge block is 90. Direct density-matrix checks use at most dimension 729. The whole suite takes about two seconds on one BLAS thread.

[Assumed normalizations] XX uses h = Sx Sx + Sy Sy, M = 1/2. Heisenberg uses h = S dot S, M = 1/2. SU(3) uses h = P, M = 1. These are the settled campaign normalizations.

## Candidate inequalities

[Numerical] The following table distinguishes proved statements from sampled extensions. A failed candidate is reported as a failure, even though the script itself successfully detects that failure.

| Candidate | XX | Heisenberg | SU(3) | Status |
| :--- | :--- | :--- | :--- | :--- |
| R <= 2 log d + (10 pi/3) M log_2(d) beta/L | Pass | Pass | Pass | Proved for XX here, only sampled for the interacting chains |
| R <= beta Delta_L | Pass | Pass | Pass | Sampled only, no general theorem |
| R convex in beta | Fails | Fails | Fails | Both signs of R'' occur |
| R concave in beta | Fails | Fails | Fails | Both signs of R'' occur |
| R <= I(two halves) | Fails | Fails | Fails | Fails by beta = 16 on the listed small rings |
| Averaged seam norm <= 4M/L | Fails | Fails | Fails | Product-state witnesses also prove failure at large L |
| Exact relative-entropy identities (20)-(23) | Pass | Pass | Pass | Proved algebraically |
| Jensen seam comparison (15) | Pass | Pass | Pass | Proved algebraically |
| Constant-angle partial swap preserves the two-copy bond | Fails | Fails | Fails | Nonzero commutator and rotated-bond mismatch |

[Numerical] The thermal and curvature sweep uses 121 logarithmically spaced aspect ratios from 0.025 to 32. The smallest sampled slack in the first row is

    XX: 1.517192495
    Heisenberg: 1.517191638
    SU3: 2.612405210

[Numerical] The extrema of R'' with respect to beta, evaluated directly as 2 Var(H_L)-Var(H_(2L)), are

| Model | L | Minimum | Maximum |
| :--- | ---: | ---: | ---: |
| XX | 4 | -0.024378482 | 0.157450102 |
| Heisenberg | 4 | -0.026822553 | 0.283208992 |
| SU3 | 3 | -0.514274221 | 2.806739373 |

## XX large rings and quadrature

[Numerical] The full spin trace is computed by the two Jordan-Wigner grids, with an even/odd occupation recursion and long-double arithmetic. Its maximum discrepancy against small-ring ED is 2.14e-14. The large-ring tests cover L = 4,6,8,16,64,256,1024,4096, with doubled length up to 8192, at tau = 0.01,0.05,0.1,0.25,0.5,1,2,4,8,32. This is linear-cost free-fermion work, not a large diagonalization. All eighty tests of the proved inequality pass. Independent integration tests the pressure-correction bound on sixteen L,beta pairs.

[Numerical] At tau = 8 the ratio and the proved ceiling are

| L | R | 2 log 2 + 5 pi tau/3 |
| ---: | ---: | ---: |
| 4 | 6.886210344 | 43.274196409 |
| 6 | 6.538655873 | 43.274196409 |
| 8 | 6.423300256 | 43.274196409 |
| 16 | 6.314967076 | 43.274196409 |
| 64 | 6.281680941 | 43.274196409 |
| 256 | 6.279609395 | 43.274196409 |
| 1024 | 6.279479964 | 43.274196409 |
| 4096 | 6.279471979 | 43.274196409 |

[Proved interpretation] The slack is substantial because the bound uses absolute Fourier-coefficient estimates. The asymptotic XX coefficient A = 1 is already known from C101 and is not inferred from this loose bound.

## Seam and information checks

[Numerical] The maximum residual in the two Gibbs relative-entropy identities and their mutual-information decomposition is 2.60e-12. Jensen's comparison holds in every test. At beta = 16,

| Model | L | R | Mutual information | Sum of marginal relative entropies |
| :--- | ---: | ---: | ---: | ---: |
| XX | 4 | 3.370382150 | 1.832322103 | 5.657854099 |
| Heisenberg | 4 | 5.581805924 | 2.096886450 | 12.043550829 |
| SU3 | 3 | 22.311179590 | 3.322852827 | 40.864071625 |

[Numerical] Constant-angle partial-swap diagnostics, on the four on-site factors of a replicated two-site bond, are

| Model | Commutator norm | Bond mismatch at theta = pi/4 |
| :--- | ---: | ---: |
| XX | 2.828427125 | 1.414213562 |
| Heisenberg | 3.464101615 | 1.732050808 |
| SU3 | 6.928203230 | 3.464101615 |

[Numerical] Sixteen tests of the classical comparison R <= log 3 pass for a nonsymmetric three-state transfer matrix. The proof needs entrywise nonnegativity and does not assume matrix positive semidefiniteness.

## Reproduction and scope

[Work record] Run inside work_C118:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 checks.py

The complete output is checks_local.json. run_summary.json contains the saved console summary. ring_tools.py contains the full-spectrum routine adapted from the read C116 source. The script's positive assertions test proved identities and the explicitly labeled sampled extension. Failed candidates are diagnosed from the saved values; their failure is not silently treated as success.

[Unperformed] run_typhon.slurm is prepared for the larger full-spectrum run. It requests one CPU, 64 GB and 24 hours, with single-threaded BLAS. The cluster option includes XX and Heisenberg through length 12 and SU(3) through length 12, enabling the vacuum-compatible SU(3) comparison L = 6 versus 12. The largest SU(3) charge block has dimension 34650, with one dense matrix occupying about 9.6 GB before solver workspace and copies. The requested runtime has not been benchmarked. The Python script refuses --cluster unless SLURM_JOB_ID is present. No Typhon connection was attempted and no job was submitted, as required by the task. These larger results remain unavailable.

[Work record] Final validation confirmed Python syntax and saved-result counts. The SLURM guard was exercised with SLURM_JOB_ID absent and refused the run with "Large ED requires a SLURM allocation." This expected refusal is not a failed numerical check.
