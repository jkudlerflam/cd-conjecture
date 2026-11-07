# C129 numerical checks

[Numerical] These are floating-point checks, not interval-certified inequalities. They test the formulas and numerical implementations. They do not establish an all-size bound.

## Completed runs

The commands below were run in work_C129 with NumPy and SciPy. The executable is python3. An initial invocation using python failed because that executable is not installed, after which all listed commands completed.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 check_c129.py --mode free --output free_check.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 check_c129.py --mode twist --maxL 10 --points 17 --output twist_check.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 check_c129.py --mode sanity --maxL 8 --output sanity_check.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 bethe_scan.py --small --output bethe_small.jsonl
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 check_c129.py --mode twist --maxL 12 --points 9 --output twist_check_12.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 check_variational.py
```

[Numerical] The free-fermion formula agrees with sector ED through L=10 to 3.56e-15, including odd lengths.

[Numerical] The first twist run uses every L from 2 to 10, all nontrivial n<=L/2, 17 flux displacements in [0,pi], and Delta in {-1.001,-1,-0.999999,-0.999,-0.99,-0.9,0,0.5,1,2,3}. There are 4,675 samples. Reflection in the minimizing flux covers the other half of the flux circle. Each sample is compared to the exact trigonometric upper bound, not only its quadratic relaxation. The largest ratio to the quadratic ceiling is 0.8105694691387021 at L=2,n=1,Delta=-1.001 and displacement pi. The zero- and fully-filled sectors have zero twist spread analytically and are omitted from the loop.

[Numerical] A second run extends to L=12 with nine flux displacements and the same anisotropies, for 3,564 samples. It has the same largest ratio and no violation. The two sample counts overlap and should not be interpreted as distinct cases. The sampling is coarser in flux than the original survey's 121-point grid. Since the reference minimum is known analytically, no sampled minimum can create a false pass by sitting above the true minimum.

[Numerical] The field test uses the 36 combinations J_perp=1, J_z in {-0.95,-0.5,0,0.5,1,1.5}, DM in {0,0.3}, and B in {0,0.4,0.8}, each at L=2,4,6,8. Physical field-free sector energies use uniform hopping phase theta=atan2(DM,J_perp). The XXZ anisotropy is J_z/sqrt(1+DM^2), not J_z. Sectors are minimized after adding the field. Both the unoptimized field-sector inequality and the twist-optimized inequality are checked. The smallest field slack is -4.44e-16, and the smallest slack of the quadratic twist-optimized bound is 0.15421256876701234. The largest physical Q is 4. The 144 values reproduce the supplied survey's rounded Q and M values within 4.751e-5.

[Numerical] Bethe checks use all L=2,...,10, every n=1,...,floor(L/2), and Delta in {-0.999,-0.99,-0.9,-0.5,0,0.5,0.9,1,1.5,3}. This gives 250 cases. Every length-L Bethe energy is compared to independent sector ED. Bethe energies for (2L,2n) are also calculated, so the paired Bethe length reaches 20. The maximum ED discrepancy is 7.238654120556021e-14 and the maximum sampled L(F_(2L)-2F_L)/M_0 is 4.0 at (L,n,Delta)=(2,1,1). The test does not claim ED verification of every doubled ring above L=10. Each record includes the Bethe residual, iteration count, and final Newton energy-correction indicator. No failures occurred.

[Numerical] The variational check evaluates both explicit trial-state bounds and their dilute relaxation for 72 cases with L=2,...,6, all nonzero n<=L/2, and Delta in {-1,-0.999999,-0.999,-0.9,0,1,3,10}. The largest signed violation is 3.56e-15, at floating-point scale. The ferromagnetic multiplet energy check has the same error. Six representative bonds, with two optimization starts each, verify the analytic M formula to 4.45e-16. Forcing the iterative Newton linear solver on L=8,n=4 at all ten Bethe-grid anisotropies gives the same energies as the dense Newton solve.

## Bethe equations used

[Proved, algebraic convention] With J=1 and centered integers or half-integers I_j=j-(n+1)/2, the minimizing-flux equations are

    p(lambda_j)-(1/L)sum_k theta(lambda_j-lambda_k)=2pi I_j/L.

For -1<Delta<1, let zeta=arccos Delta. The code uses

    p(lambda)=2atan[tanh(lambda)/tan(zeta/2)],
    theta(lambda)=2atan[tanh(lambda)/tan(zeta)],
    E=Delta L/4-sum_j sin(zeta)^2/[2(sinh(lambda_j)^2+sin(zeta/2)^2)].

For Delta=1 it uses p=2atan(2lambda), theta=2atan(lambda), and bare energy -1/[2(lambda^2+1/4)]. For Delta>1 it uses the trigonometric continuation with zeta=arccosh Delta and continuously unwrapped phases. The bare energy is -sinh(zeta)^2/[2(sinh(zeta/2)^2+sin(lambda)^2)]. At the analytically minimizing flux, the centered equations are the same for odd and even L. The odd-length ground-state identification is supported by the ED check, and is not claimed here as an extension of the cited literature theorem.

[Numerical] Newton steps use the analytic Jacobian. A line search decreases the logarithmic-equation residual, and roots are checked for strict order. Dense linear solves are used through 256 roots. Larger systems use GMRES with a diagonal preconditioner. The GMRES branch was exercised on small systems, not on the large requested grid. Summation uses NumPy long double, but transcendental functions and roots use double precision. The final Newton energy correction is a diagnostic, not a rigorous error bound. Cancellation in L(E_(2L)-2E_L) requires special attention at large L, particularly near Delta=-1.

## Prepared but not run

[Assumed, computational work pending] run_typhon.slurm is an array of 11 CPU jobs. Ten jobs use the requested fixed anisotropies. The eleventh uses Delta=-1+c/L^2 for c in {0.1,1,10,100}. Lengths run from 16 to 4096 by doubling, with paired rings through 8192. Densities include fixed n=1,2,3,4, n approximately sqrt(L), the half-filled sector and its neighbor, and a fixed-density subsample. This is not the entire set of sectors at those lengths.

The script requests four CPUs and 12 GB per job, and assumes python3 with NumPy and SciPy is available in the submission environment. A site-required account or CPU partition can be supplied to sbatch. There is no GPU dependency. No network command is included and no submission was made from this sandbox.

Failures are retained as records with status failed. A completed scheduler job therefore does not certify that all solves succeeded. Inspect the records and the scaled correction indicators before interpreting a growth pattern. The local checks found no counterexample, but the requested large-L falsification tests remain unperformed.
