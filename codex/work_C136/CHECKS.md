# C136 numerical checks

2026-11-07. [Numerical] Every numerical statement here is floating point. No spectral enclosure or independent referee certification is claimed.

## Bethe and ED checks

The local solver is copied from C134, which adds return_roots to C129's solver. check_c129.py and half_scan.py are local copies as well. No Python imports or bytecode writes go into the predecessor folders.

The commands run were

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 check_c136.py --mode small --output small.jsonl
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 check_c136.py --mode moderate --output moderate.jsonl
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 check_c136.py --mode spot --output spot.jsonl
```

[Numerical] All three exited successfully, including reruns after adding the all-filling constant. They contain 216, 672, and 144 paired cases, respectively. The anisotropies are the six requested values. The first mode checks every n=1,...,floor(L/2) at each L=2,...,12. The second does the same at L=16,17,32,33,64,65. The third uses n in {1,2,round(sqrt(L)),floor(L/4),floor(L/2)-1,floor(L/2)} at L=128,129,256,257. Every pair uses (2L,2n), also at odd L.

[Numerical] ED cross-checks every small-ring energy with L<=12 and every doubled-ring energy with 2L<=12. Maximum absolute energy discrepancy is 1.24345e-14. The largest residual is 3.975e-14, the largest scaled last-Newton energy correction is 1.763e-11, and the largest discrepancy against the cancellation-free free defect formula is 1.798e-11. These diagnostics are not error certificates.

[Numerical] All cases check the C129 free trial bound, dilute bound, ferromagnetic band, and the new all-filling compact constant. Below half filling they also check the smaller finite-boundary constant, and at exact half filling they check C134's constant. The minimum variational slack is -5.685e-14, at floating-point scale. The largest Q=L Lambda/M_0 is 2.3431457505076203 at Delta=0,L=2,n=1. It is the known free value 8-4sqrt(2), not a proposed universal ceiling.

## Independent continuum diagnostics and failures

The commands were

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 check_c136.py --mode intermediate --output intermediate.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 check_c136.py --mode intermediate --index 1 --max-order 1296 --output intermediate_refined_1.json
```

[Numerical, incomplete] The first run had 14 unresolved cases out of 54. The refinement at Delta=-0.999 resolved three additional cases. The combined result is 43 successful and 11 unresolved cases, recorded in intermediate_combined.summary.json. Initial outputs are retained in intermediate_initial.json and intermediate_initial.summary.json. The code was subsequently extended to retain each adaptive-order attempt and to support a larger cluster refinement ceiling.

The unresolved output consists of nine Delta=-0.9999 cases with

```text
f(a) and f(b) must have different signs
```

and two Delta=-0.999 cases with

```text
L=17, n=8:  Quadrature unresolved: 9.313208249039917e-08
L=32, n=15: Quadrature unresolved: 5.762225274708044e-08
```

These are failures of the independently discretized Fredholm integral equation and continuum-boundary search, not failures of the finite Bethe solve or observed violations of a proved inequality. At small alpha=pi-gamma, the scattering kernel has width of order alpha. Underresolved Gauss discretization can distort its norm and prevent the mass root from being bracketed. Increasing the quadrature order repairs some cases, but all eleven are kept as unresolved. A small residual of a discretized equation would not by itself resolve this issue.

[Numerical] The successful cases test the exact dressed energy and mass identities, cumulative-root bound, density comparison, density mismatch, curvature, quadratic tangent gap, velocity variation, exponentially localized boundary corrections through two derivatives, direct midpoint error, and the finite-boundary constant. Their maxima are

| Diagnostic | Largest observed value |
| :--- | ---: |
| Dressed quadrature identity error | 1.777e-15 |
| Mass identity error | 1.361e-15 |
| Change under continuum refinement | 2.577e-9 |
| Cumulative-root error divided by its bound | 1.000000000000002 |
| Density discrepancy divided by its bound | 0.084902 |
| Mass mismatch divided by its bound | 0.055780 |
| Midpoint error divided by its bound | 0.041430 |
| Boundary correction divided by its exponential bound | 0.037751 |
| Sampled velocity variation divided by its bound | 8.253e-7 |

Sampled variation is a lower estimate of full variation. The proof supplies the bound, not this sample. The smallest sampled r_Q/r_infinity is 0.139513, above the conservative kappa in every successful case.

## Fields, kernels, constants, and clipping

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 audit_c136.py
```

[Numerical] This exited successfully. It checks 450 field/DM/length combinations through L=6 with doubled ED through 12, and 216 sector twist cases. The smallest field-reduction slack is -4.441e-16; the smallest combined field/twist slack is 0.2056167584. The largest twist spread divided by J pi^2/(4L) is 0.8105694691. Endpoint sector energies at Delta=-1 agree with -JL/4 to 3.553e-15.

[Numerical] Ninety further cases at L=8,...,12 check the optional explicit high-field sector-placement estimate in PROOF Section 9. The selected density never exceeds 1/8 in that sample. The convex-combination sector lower bound has minimum slack -3.553e-15. This extra estimate is not needed for the all-filling theorem.

[Numerical] Independent integrals check the kernel L1 norm and its beta=1/2 exponential moment to 1.34e-15. A 301-point logarithmic delta grid in [1e-6,1] checks the closed majorant K_all<=10^15 delta^(-25/3); the largest ratio is 0.001512. The C134 crossover constants are evaluated down to epsilon=1e-6 and reproduce the powers derived in PROOF Section 8. For example, at epsilon=1e-4, the band length for A=1 is 100, while C134's cutoff is approximately 62556.

[Numerical, scope limitation] The new all-filling constants make its sufficient cutoff length much larger than the local interacting Bethe grid. That interacting large-size clipping branch was not tested directly. C134's sharper half-filled clipping checks are used where applicable. Four additional free cases at L=2^20 and 2^22, with n=L/2 and n=L/2-1, check the new clipping branch using geometric sums and analytic root counts. These are constant-cost scalar evaluations, not large local diagonalizations or root arrays. They test the midpoint and tail split and the clipping jumps. They do not test interacting Fermi-boundary corrections at those sizes.

## Prepared for Typhon, not run

[Assumed, pending] run_typhon.slurm is a six-job CPU array, one job per anisotropy. It runs every L=2,...,256 and every sector n=1,...,floor(L/2), totaling 98,304 pairs and paired rings through 512. run_intermediate_typhon.slurm refines the 54 continuum cases with adaptive Gauss orders through 4374. Both request one CPU, 8 GB, and 24 hours; supply a site-required account or CPU partition when submitting. Both retain failures and exit nonzero when any case remains unresolved.

No job was submitted and no network command was used. The full requested all-sector grid through a few hundred sites remains pending. There is no proved constant uniform in Delta to check as a final global ceiling.
