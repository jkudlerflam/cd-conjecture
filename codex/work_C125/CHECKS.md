# C125 numerical checks

2026-11-06. [Numerical] All numerical claims here refer to the scripts and outputs in this folder. Floating-point sign tests are not certified spectral enclosures.

## Reproduction and scope

[Numerical] Run from work_C125 with one BLAS thread:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 audit_checks.py
    OPENBLAS_NUM_THREADS=1 python3 spin_certificates.py

The first two scripts write checks.json and audit.json. The direct spin certificates were also saved to spin_certificates.txt. The main script takes about four seconds on this machine; the additional audit takes about one second. The scripts import only Python standard libraries, NumPy and SciPy. No network access is used.

[Numerical] Full-spectrum ED covers spin-1/2 XXZ at Delta=-0.9,-0.5,0,0.5,1,2 and lengths 2 through 8, with ground energies extended to even L=12. SU(3) full spectra cover lengths 2 through 8, with its ground energy also computed at L=9. The attractive-XXZ monotonicity audit computes full spectra at L=4,6,8,12. Magnetization or color-occupation blocks reduce the largest local full diagonalizations to dimension 924. These are small ED calculations.

[Numerical] Small-Trotter transfer spectra use beta=0.2,2,8 and m=1,2,3,4 for spin 1/2, and m=1,2,3 for SU(3). Column dimension is d^(2m), so the largest transfer matrix has dimension 729. An independent contraction of the original length-four Trotter product verifies the spatial-transfer trace identity for XX, Delta=-0.9, Heisenberg and SU(3).

[Assumed scope limit] These small m values are not a controlled approximation to a large-beta quantum transfer spectrum. Their purpose is to test positivity and the contraction convention. Large-beta XX is instead tested against its exact infinite-column free-fermion product.

## Passes and failures

| Test | Outcome |
| :--- | :--- |
| [Numerical] XX Jordan-Wigner trace versus spin ED | Maximum absolute log-trace error 2.28e-13. |
| [Numerical] XX crossed infinite product versus Jordan-Wigner trace | Maximum absolute log-trace error 5.59e-9, at L=1024, beta=32768, where log Z is of order 10^7. Doubling the mode cutoff gives no resolved change. |
| [Numerical] XX C118 ratio bound | Passes all 48 samples, L=4 through 8192 and doubled lengths through 16384, beta/L=0.01 through 128. |
| [Numerical] Trivial and gluing ratio bounds | Pass every sampled model, length and temperature. |
| [Numerical] Positive small-Trotter spectra, Delta>=0 | Pass up to floating-point roundoff. |
| [Numerical] Positive small-Trotter spectra, Delta<0 and SU(3) | Fail. Negative and complex eigenvalues occur. |
| [Numerical] R monotone in even L at Delta=-0.9 | Fails on exact spin traces, including L=4 versus 6. |
| [Numerical] D nonnegative on all even SU(3) rings | Fails at L=4. |
| [Numerical] Theta convex and c_A nonincreasing at Delta=-0.9 | Both fail on L=4,6,8,10,12. |
| [Numerical] Theta convex at the other tested anisotropies | Passes on the sampled even lengths. |
| [Numerical] SU(3) on vacuum lengths | Both signs pass on 3,6,9. |
| [Numerical] Supplied Sym^2 energies, plus exact first point | Both signs pass for all six N families. |
| [Numerical] Global linear Matsubara spectrum condition | Fails for the exact XX locally detected branch. The corrected nonlinear LR rate passes every sample. |
| [Numerical] Countermodel cone, covariance Hessian, R derivative and rank bound | Pass in all 81 momentum-resolved audit samples. |
| [Proved finite arithmetic check] SU(2) Racah recurrence and harmonic sum | Exact rational identities pass for k=1,...,12. This checks the finite formulas, while the general proof is in PROOF section 2.3. |

[Numerical] Representative small-Trotter errors and spectral signs at beta=2 are shown below. Eigenvalues are divided by the largest absolute eigenvalue. The last column is log Z_m(4)-log Z(4), so it also records the size of the Trotter error rather than concealing it.

| Model | m | Smallest real part | Largest imaginary magnitude | Log-trace error |
| :--- | ---: | ---: | ---: | ---: |
| XX | 4 | 6.98e-12 | 1.70e-16 | 0.00587127 |
| XXZ Delta=-0.9 | 4 | -0.0105640 | 0.00644676 | 0.00134409 |
| XXZ Delta=0.5 | 4 | 2.61e-7 | 1.05e-17 | 0.0178289 |
| Heisenberg | 4 | 5.17e-6 | 1.54e-17 | 0.0395091 |
| SU(3) | 3 | -0.0354088 | 0.150133 | 0.450366 |

[Numerical] Exact-ring witnesses are

    Delta=-0.9, beta=21.9287213768:
       R_4=1.07347733534, R_6=1.13747364150,
       R_6-R_4=0.06399630616.

    SU(3): E_4=-2, E_8=-5.40634543903,
       D(4)=-1.40634543903.

    XX, beta=1, Matsubara k=100:
       omega=628.318530718,
       exact spatial rate=8.99347827408,
       proposed omega/(4eM)=115.572734979.

[Numerical] For the supplied Sym^2 data, the first c_A/M values after adding L=N are 5.092958179, 5.798189641, 6.251154306, 6.567765808, 6.982149743, 7.125295634 at N=2,3,4,5,7,8. The corresponding infrared A/M values are 3.141592654, 4.468042885, 5.235987756, 5.744626567, 6.382918407, 6.597344573. The code uses secant slopes in L^2, not ordinary equal-L second differences.

## Numerical limitations and corrections made during the run

[Numerical] Direct adaptive integration in k missed the narrow thermal region near the XX Fermi point at beta=4096 and 32768. It initially produced a log-trace discrepancy about 0.033 and an erroneous zero value for delta_beta. Rewriting that integral in y=beta cos k resolves the region and gives the final errors recorded above. The product cutoff was not the cause. The corrected implementation remains in checks.py.

[Numerical] The first compact-boson enumeration attempted beta=128 with a conformal-level cap 500 and exhausted that cap. Those cases were removed, the bound parameter was aligned with b=1/2 in the proof, and the final finite enumeration uses beta=2,8,16,32. No final case exhausts the enumeration cap. The large-n and large-beta conclusions use the analytic counting calculation, not an extrapolation from an exhausted table.

[Assumed scope limit] Bethe input energies are used as supplied and not recomputed. Their numerical error bars are not certified here. The SU(3) vacuum sequence has only three locally computed points. The all-spin and all-n statements rely on their analytic derivations, not these finite samples.

[Work record] Larger ED is available through run_typhon.slurm and checks.py --cluster. The script refuses the cluster mode without a SLURM job identifier. It has not been submitted or run. No network shell command or Typhon connection was attempted.
