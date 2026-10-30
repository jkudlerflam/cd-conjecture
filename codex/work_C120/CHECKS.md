# C120 numerical checks

2026-10-29.

[Numerical] Run in this folder with

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 checks.py

The last complete run took about ten seconds. Results are in checks.json, with a short machine-readable summary in summary.json. The largest dense physical matrix has dimension 729. All spectra and Gibbs traces in the small-ring calculations are complete. No charge-sector projection or low-energy truncation is used.

[Assumed] Normalizations are inherited from the background and C118. XX uses h = Sx Sx + Sy Sy and M = 1/2. Heisenberg uses h = S dot S and M = 1/2. SU(3) uses h = P and M = 1.

## 1. Candidate outcomes

[Numerical unless explicitly marked proved] The dimension ceilings in this table are tested as proposed replacements, not silently assumed.

| Candidate | XX | Heisenberg | SU(3) | Interpretation |
| :--- | :--- | :--- | :--- | :--- |
| Dmax(rho||sigma) <= 4 log d | Fails | Fails | Fails | Four physical seam sites do not bound this divergence |
| D_alpha(rho||sigma) <= 4 log d | Fails | Fails | Fails | Already fails at alpha = 0.9 in the sample |
| R <= beta<V>_rho + D_(1/2)(rho||sigma) | Fails | Fails | Fails | The proved inequality has the opposite direction |
| R <= beta<V>_sigma - D_alpha(sigma||rho) | Pass | Pass | Pass | Valid upper bound, with the wrong large-L scale at fixed order on XX |
| R <= 4 log d | Fails | Fails | Fails | Rejects a fixed four-site capacity without a normalization term |
| R <= trace distance | Fails | Fails | Fails | Distance is at most one |
| R <= -2 log root fidelity | Fails | Fails | Fails | Omits the unnormalized overlap |
| Support D0 = 0 at finite beta | Pass | Pass | Pass | Full rank is known from the exact Gibbs spectra |
| Conditional-min-entropy dimension witness | Pass | Pass | Pass | Feasible fixed marginal checked, not full finite-temperature optimization |
| Petz mixed-trace identity | Pass | Pass | Pass | Proved identity |
| Fuchs-van de Graaf inequalities | Pass | Pass | Pass | Root-fidelity convention |
| Signed-history inequality, PROOF (19) | Pass | Pass | Pass | Uses the actual history dimension |
| Replace history dimension by d^2 | Fails | Fails | Fails | Fails already in the finite-depth networks |

[Proved and numerical] An additional faithful qubit family tests the dimension premise for every listed Renyi order. Both Petz and sandwiched values grow with the qubit parameter. A separate completely positive qubit-map family violates the trace-ratio ceiling log 4. These are algebraic counterexamples to the auxiliary premises, not local-chain counterexamples.

## 2. Small-ring state calculations

[Numerical] Compared lengths are L = 4 and 2L = 8 for XX and Heisenberg, and L = 3 and 2L = 6 for SU(3). Temperatures are beta = 0, 0.1, 0.5, 1, 2, 4, 8, 16, 32. The six Petz orders are 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, giving 162 forward and 162 reverse evaluations. The maximum residual in the exact Petz trace identity is 1.92e-13. The original Umegaki seam identities are also checked.

[Numerical] The following beta = 16 values show where the slope resides.

| Model | R | D_(1/2), Petz | -2 log F | Trace distance | 2 log[G(X,Y)/Z1] |
| :--- | ---: | ---: | ---: | ---: | ---: |
| XX | 3.370382 | 0.986428 | 0.522305 | 0.627595 | 2.848077 |
| Heisenberg | 5.581806 | 0.896006 | 0.448427 | 0.600924 | 5.133379 |
| SU(3) | 22.311180 | 1.233411 | 0.616706 | 0.678440 | 21.694474 |

[Numerical] The false forward Petz replacement beta<V>_rho+D_(1/2) at beta = 16 gives approximately -3.133366, -7.662625 and -20.642334. Each is below R, as required by the proved reverse direction. At beta = 32, D_(0.9) equals 4.432828, 4.471519 and 6.167056, all above the respective four-site capacities.

[Numerical] The conditioned direct Dmax witnesses are

| Model | beta | Dmax | 4 log d |
| :--- | ---: | ---: | ---: |
| XX | 4 | 10.926100 | 2.772589 |
| Heisenberg | 4 | 15.654739 | 2.772589 |
| SU(3) | 2 | 15.017146 | 4.394449 |

[Numerical method] Powers of Gibbs states are constructed from the Hamiltonian eigensystems. Petz traces are nonnegative spectral-overlap sums evaluated with log-sum-exp. Root fidelity is a singular-value sum of weighted overlaps, avoiding numerical square roots of nearly singular density matrices. Dmax uses the largest singular value of the relative square-root matrix. Direct Dmax is recorded only when beta times the spectral width is at most 25. Sandwiched calculations also require that width multiplied by max(1,(1-alpha)/alpha) be at most 25. Low-order singular-value powers can amplify roundoff, so values outside this range are omitted, not regularized by an arbitrary eigenvalue floor. No omitted value is counted as a numerical success. There are 74 conditioned sandwiched evaluations, 17 direct Dmax pairs, and 15 finite-temperature conditional-min witnesses.

[Numerical scope] Finite-temperature conditional min-entropy tests compute Dmax(rho||I_C tensor rho_D) for one site and verify the analytic interval [-log d,log d]. This supplies a feasible dimension witness. It is not generally the optimized Hmin. At zero temperature the pure-state formula is evaluated without an optimization. The one-site ground-state Hmin is -log 2 for both XX and Heisenberg and -log 3 for SU(3), in both the split and joined systems. The positive doubling defects are respectively 0.215301195, 0.348906591 and 1.394448725.

## 3. Actual spatial transfer networks

[Numerical] The network implementation uses the explicit histories and idle-site constraints of PROOF Section 8. Two layers and depths n = 1, 2, 3 are used for XX and Heisenberg. Three layers and depths n = 1, 2 are used for SU(3). Each is checked at beta = 0.5, 2, 8. The largest history matrix has dimension 729. Independent temporal contractions validate both the numerator and denominator traces, with maximum absolute residuals 6.22e-15 and 6.11e-14 respectively. These tests compare exact contractions of a finite product formula. They do not assume that a small Trotter depth approximates a cold quantum trace accurately.

[Numerical] Representative beta = 8 results are

| Model | Depth n | History dimension | R^(n) | Exact ED R | Negative-pair penalty |
| :--- | ---: | ---: | ---: | ---: | ---: |
| XX | 2 | 16 | 2.225141 | 1.507576 | 0 |
| Heisenberg | 2 | 16 | 2.755986 | 2.747529 | 0 |
| SU(3) | 2 | 729 | 3.171792 | 11.155371 | 0.508770 |

[Numerical] All sampled blocked XX and Heisenberg transfer matrices are entrywise nonnegative. The SU(3) sample at depth two has negative entries and negative pair products. The signed-history bound holds in all 24 cases. The large finite-depth SU(3) error at beta = 8 is reported explicitly. No continuum-Trotter extrapolation or uniform sign-penalty estimate is inferred from these depths.

## 4. XX free fermions at larger lengths

[Numerical] Thirty-five full thermal ratios cover L = 4, 8, 32, 128, 512, 2048, 8192 and beta/L = 0.1, 1, 4, 16, 64. Both Jordan-Wigner parity sectors are retained. The code separates the known ground energy and recursively sums all even/odd excitations, then uses the cancellation-free doubling defect in PROOF (20). Agreement with the earlier direct full-trace recursion on L <= 32 is 1.31e-11. The new expression avoids subtracting two bulk quantities of order beta L at large L.

[Numerical] At L = 8192, R is 12.566363 at beta/L = 16 and 50.265484 at beta/L = 64. Both exceed the tested fixed capacities and obey the established C118 XX upper bound. These large-ring tests check scalar consequences of the proposed capacity substitutions. They do not compute finite-temperature many-body trace distance or fidelity on length 16384.

[Numerical] Separate ground-state determinant checks reach L = 512, paired length 1024. They evaluate root fidelity, trace distance, support D0, both Renyi definitions at orders 0.05, 0.25, 0.5, 0.75 and 0.9, and the pure-state conditional min-entropy of four sites. The determinant fidelity at L = 4 agrees with full ED below 1e-12. At L = 512,

    root fidelity = 0.432261760,
    D0 = -2 log F = 1.677447895,
    Hmin(four sites | rest), split = -1.596227105,
    Hmin(four sites | rest), joined = -1.596254976.

[Proved interpretation] The pure states have different supports, so Dmax is infinite and D0 need not equal the finite-temperature value zero. Their Renyi formulas are evaluated analytically from the overlap. No finite-temperature limit is inferred by numerically deleting small eigenvalues. No fitted scaling exponent is claimed.

[Unperformed] No large interacting ED, finite-temperature large-ring overlap calculation, or full finite-temperature conditional-min-entropy optimization was run. No cluster access or external research-file modification was needed. All numerical conclusions above state the actual tested sizes and temperatures.
