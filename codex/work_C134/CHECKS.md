# C134 numerical checks

[Numerical] All reported values are floating point. No numerical result substitutes for a proof or an interval-certified spectral bound.

## Commands executed

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 half_scan.py --mode small --output half_small.jsonl
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 half_scan.py --mode moderate --output half_moderate.jsonl
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 check_analytic.py
```

The third command was rerun after adding the forced iterative-solver checks. All runs passed.

## Scope and measured errors

[Numerical] The anisotropy grid is {-0.999,-0.99,-0.9,0,0.5,0.9,0.99,1}. The small run has 88 cases, one for every L=2,...,12 and grid point. Each length-L energy is compared with sparse or dense sector ED. Spin-flipped sectors are also diagonalized. Doubled Bethe rings reach 24 sites; ED does not cover those beyond 12. Maximum Bethe/ED discrepancy is 8.43769498715119e-15; maximum spin-flip discrepancy is 9.769962616701378e-15.

[Numerical] The moderate run adds 64 cases at L=16,17,32,33,64,65,128,129. Its largest doubled ring is 258. Across both runs the largest logarithmic-equation residual is 2.776e-14 and the largest scaled last-Newton energy-correction indicator is 1.072e-11. All 152 solves of paired energies passed. The sampled maximum Q is 4 at the two-site Heisenberg point. The table in RESULT.md gives the largest completed even and odd lengths.

[Numerical] The zero-anisotropy doubling formula agrees to 9.24e-13 in Q. In particular the odd half-filled defect vanishes. Positive Delta has negative odd defects on the sampled grid, whereas the even defects are positive. Comparing these odd results to C129's matched-sector defects would be a convention error.

[Numerical] At the 76 nonpositive-anisotropy points in the two scans, cumulative root and density errors obey the bounds from PROOF.md. Eighteen of these lengths reach the analytic cutoff regime; their sampled derivative variation obeys its ceiling. The exact bound itself is proved analytically, since a sampled variation alone would be a lower estimate of total variation.

[Numerical] check_analytic.py additionally uses Delta=-0.9,-0.5,-0.1,0 and L=4,5,16,17,64,65, giving 24 dressed-energy identity cases and 14 cutoff cases. It finds

    maximum dressed-energy identity error       1.066e-14,
    maximum dressing-equation error             1.961e-16,
    maximum scattering-phase identity error     8.456e-16.

The checks integrate the full-line convolution independently of the energy sum. For the supersymmetric point, direct ED through L=12 gives maximum special-twist energy error 2.221e-14 and odd periodic energy error 5.330e-15. The C=2pi^2/9 defect ceiling passes. At L=12, all eight requested anisotropies give identical energies between dense Newton and a forced GMRES solve at the reported precision. The output is in analytic_checks.json and analytic_detail.json.

## Solver provenance and precision

[Numerical] bethe_scan.py and check_c129.py were copied from C129 into this folder. The only solver change is an optional return_roots argument and corresponding return path, used for diagnostics. The logarithmic equations, Newton method, analytic Jacobian, line search, and energy summation are unchanged. half_scan.py fixes the sectors for this task and adds the requested Delta=0.99 grid point. Python imports and bytecode stay in this folder.

[Numerical] Roots and transcendental functions use double precision; energy sums use NumPy long double as in C129. A small residual or last-Newton correction is not an error enclosure. At large L the extensive energies cancel, amplifying their errors in Q. The large output therefore retains residuals and the scaled correction indicator. The free-grid comparison is performed using a cancellation-free exact formula.

## Prepared, not executed

[Assumed, pending] run_typhon.slurm requests an array of eight CPU jobs, four CPUs and 12 GB per job, 24 hours. Each job handles one requested anisotropy at L=256,257,512,513,1024,1025,2048,2049,4096,4097, with paired rings through 8194. It requires Python with NumPy and SciPy. A site-required account or CPU partition can be supplied at submission. No GPU is required.

```sh
sbatch run_typhon.slurm
```

This command was not executed. No cluster connection was attempted. The script writes a failure record on an exception and exits nonzero if any case failed. A future run must inspect those records and cancellation diagnostics before quoting its maximum. No several-thousand-site result is claimed here.
