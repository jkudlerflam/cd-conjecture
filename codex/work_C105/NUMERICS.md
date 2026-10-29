# C105 numerical record

[Numerical] Only lightweight checks ran locally. All numerical code is Python with NumPy and SciPy, stored in this folder. No stochastic sampling was used.

## Exact half-line XX

[Numerical] Run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 xx_halfline.py

Output is xx_halfline_results.json. The method is the infinite Jacobi-operator spectral transform, with a three-term polynomial recurrence and Gauss-Legendre quadrature on [0,12]. It evaluates all bonds through x = 8193 at orders 400 and 800. Runtime was 0.12 seconds. The maximum difference between the two quadratures is 2.075e-15. The analytically bounded omitted spectral tail is at most 1.799e-33. The agreement of quadratures is an empirical integration check, not a rigorous error enclosure.

[Numerical] An independent first-bond identity is checked in certificate_checks.py using the Catalan constant expressed by trigamma functions. Its value, 0.05291698652698151, differs from the recurrence quadrature by 6.1e-15.

## Small interacting and interval-decomposition checks

[Numerical] Run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 checks.py

Output is small_checks_results.json. The models are XX, Heisenberg, XXZ with anisotropy 1/2, and SU(3) swap. The largest sector dimensions are 924 for the spin-half models at N = 12 and 1680 for SU(3) at N = 9. Ground-state vectors are obtained by deterministic Lanczos initialization. The calculation uses balanced charge sectors. It does not establish the absence of a lower state in a different sector from first principles.

[Assumed, uniform reference values] The bulk energies are -1/pi for XX, 1/4-log(2) for Heisenberg, -3/8 for XXZ(1/2), and 1-log(3)-pi/(3sqrt(3)) for SU(3) swap. Their A values are 1, pi/2, 3sqrt(3)/4, and 4pi/3. The last SU(3) energy is the Bethe-ansatz value, equivalently 1-(2/3)[psi(1)-psi(1/3)]. No bulk energy was fitted from the small weighted chains.

[Numerical] Uniform block energies are computed for all lengths needed by the interval linear programs. For each quadratic or ramp weight, scipy.optimize.linprog optimizes the allowance over every interval. Assertions check

    actual defect <= optimized allowance <= nested layer allowance,
    sum of interval coefficients >= maximum weight.

All passed. The largest eigenvector residual was 6.667e-12. The largest error in sum_x w_x<h_x> = E was below 3.6e-14. The suite took 0.18 seconds. All local profiles, including the raw oscillations, are in the output file.

[Numerical] The smoothing kernels are (1,2,1)/4 for period-two oscillations and (1,2,3,2,1)/9 for SU(3). These finite filters do not remove finite-size corrections. In particular, their small-system values should not be reported as converged CFT coefficients.

## Certificate checks

[Numerical] Run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 certificate_checks.py

Output is certificate_results.json. It verifies the exact XX weight factorization and the lower-energy inequality for N = 8 through 256, computes the explicit full-rank local-marginal perturbations for r = 2 through 8, checks the first-bond formula, and tests the naive centered-energy extension of the XX certificate on all four models. That last candidate fails for Heisenberg. The failure is retained in the output and is not treated as a program error.

## Unrun large finite-chain job

[Proved, infrastructure observation] Both SSH attempts to Typhon failed before authentication. cluster_job.sbatch has not been submitted, and no remote directory was created. The prepared command is

    python3 checks.py --cluster --output large_gaussian_results.json

It evaluates finite Gaussian quadratic and ramp chains at N = 512,1024,2048,4096,8192. The code checks for SLURM_JOB_ID before allowing this path. The submission script leaves the partition unspecified because cluster discovery failed. Its memory request is 8 GiB and its time request is one hour. Those resources have not been benchmarked on Typhon.

[Assumed, limitation] Further large interacting calculations would require a separately prepared and submitted SLURM job. They are not represented by the unrun Gaussian job. The completed interacting evidence is the small suite above.
