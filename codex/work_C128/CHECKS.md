# C128 numerical checks

2026-11-06.

[Numerical reproduction] From work_C128 run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 checks.py

The script uses Python, NumPy and SciPy, writes checks.json in this folder, and took about half a second for the final local run. It imports no neighboring run's code and reads no external numerical data. Formulas from C125 are implemented independently. Floating-point checks are not interval-certified proofs.

## Finite-Trotter histories

[Numerical] For Delta=0,0.5,1,2, beta=0.2,2,6,8, m=1,2,3,4 and L=2,4,6, the code constructs the explicit gate, both history matrices K_0,K_1, the positive matrix U, and the original spin-space checkerboard product. There are 192 trace comparisons. The largest history matrix has dimension 256. Maximum relative error in the nonnegative contraction is 5.5722e-12; the positive-spectrum trace agrees with the nonnegative contraction to 9.1872e-15. Gate formulas and realignment signs are checked separately at four step sizes per model. The positive-spectrum sandwich log Ztilde_L<=R_L<=2 log Ztilde_L passes.

[Numerical] Representative values at beta=2,m=4,L=4 follow. The last column is the actual Trotter error, not a contraction error.

| Delta | log Ztilde_L | R_L | log Z_m minus log Z |
| :--- | ---: | ---: | ---: |
| 0 | 0.07006794 | 0.13755813 | 0.00587127 |
| 0.5 | 0.09344679 | 0.18323667 | 0.01782895 |
| 1 | 0.14619512 | 0.28607740 | 0.03950913 |
| 2 | 0.35760917 | 0.63942162 | 0.09952505 |

[Numerical scope] These small Trotter numbers do not approximate the low-temperature continuum uniformly. At beta=8 some log-trace errors exceed ten. They test the tensor contraction, positivity and candidate inequalities on the finite-m model. No inference of an interacting continuum theorem is drawn from them.

[Numerical] The generous scalar Fock candidate uses S=2, V'=4pi M and K=log 3-2 log(1-exp(-1)). It passes every sampled matrix, with minimum interacting log-margin 1.62348. This candidate is proved here only for continuum XX. The stronger XX-specific bound log Ztilde<=log 3+beta[csc(pi/L)-L/pi] fails at 1, 9 and 14 samples for Delta=0.5,1,2 respectively. For Heisenberg, beta=8,m=4,L=4 gives log Ztilde=2.74279041 versus bound 2.22640443. This finite-m failure is distinct from the exact continuum failure established by the ground-energy slope in PROOF section 3.

## Exact XX and the new estimates

[Numerical] Thirty-five samples use L=2,4,6,16,64,256,1024 and beta/L=0.01,0.1,1,4,32. The full parity-resolved Jordan-Wigner spin trace is compared with the positive infinite-column product, truncated at 20000 modes and repeated at 40000. No finite Trotter approximation enters this test. The script tests the product domination, convex midpoint bound, the new trace and ratio bounds, C119's bound, and the linear Fock majorant.

[Proved numerical error control] For the omitted odd-frequency logarithmic factors, log(1+x)<=x and exp[-asinh x]<=1/x give the conservative bound

    tail <= 4 (beta/(2pi))^L N^(1-L)/(L-1)

after N modes. The even-frequency tail is smaller. Even-occupation projection preserves the same upper bound on the logarithmic trace change, because the odd-occupation sum of the retained modes is at most the even-occupation sum. The displayed trace bound and linear Fock margins remain positive after subtracting this tail allowance. All exact-versus-crossed errors are below the allowance plus 1e-7 for floating-point cancellation at large log Z.

[Numerical] The largest raw cutoff error is 0.005188, at the slowly convergent L=2 end. It is covered by the tail estimate. After subtracting that allowance, the maximum residual is 5.59e-9, at a large log trace. Thus these data are not reported as uniformly nine-digit evaluations of the infinite product. The integral identity I_L=pi/[2 sin(pi/L)]-L/2 is independently integrated at six real lengths between 2 and 64, with maximum error 1.55e-14. Exact XX ground defects satisfy the new constants through L=1024. The Heisenberg four-site ground energy is independently diagonalized and its comparison uses the inherited exact bulk energy 1/4-log 2.

## Event count and failed independent clocks

[Numerical] Thirty-six exact four-site spectral checks evaluate the local-event fugacity generating function for all four anisotropies, beta=0.2,2,8, and fugacity z=1.25,2,4. All satisfy log F_i(z)<=beta(z-1), with minimum margin 0.04693.

[Numerical] The independent-column comparison is tested directly on the explicit gate at steps 0.1,0.01,0.001,0.0001 for all four models. At step 0.0001, c/a=4.999999996e-5 and the independent event weight (t/2)^2=2.5e-9, giving ratio 19999.99998. This ratio is the same for every Delta in the chosen normalization. For the two-state conditioned bridge at remaining time 0.001, the rate is 1000.00008333, against the unconditioned activity 0.5. These numerical witnesses illustrate the analytic divergences, not extrapolated finite-size trends.

## Frequency tests and conditional closure

[Numerical] Fifteen exact XX branch values cover beta=1,10,100 and Matsubara index k=1,2,10,100,1000. At beta=1,k=100, the exact gap is 8.993478274, the false linear claim omega/(4eM) is 115.572734979, and the corrected logarithmic rate is 5.749900072. The nonlinear rate passes all fifteen samples.

[Numerical] The nonlinear bosonic sum appearing in the conditional theorem is tested at L=2,4,8,16 and beta=0.1,1,10,100. Every truncated sum is below the proved integral bound. These are checks of the scalar summation step, not evidence that an interacting lattice has the assumed many-state majorization.

## Corrections and unperformed work

[Numerical correction] The first run stopped with an OverflowError while evaluating the explicit tail bound as a product of large powers at large L. The bound is now evaluated in logarithms. The corrected full suite passes. No failed assertion or failed candidate is suppressed from the final account.

[Work record] run_typhon.slurm runs the same checks with m=5 added, using the --cluster guard. It was not submitted. The sandbox cannot access Typhon, and no SSH or network shell command was attempted. No large local calculation or independent referee review was performed.
