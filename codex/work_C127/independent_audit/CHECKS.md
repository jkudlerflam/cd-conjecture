# C127 independent numerical checks

2026-11-06.

[Numerical] Run from this subfolder with

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 audit.py

The scripts use NumPy and SciPy and write only beside themselves. Both completed successfully after isolation from concurrent parent-folder edits. No other run's output is imported. Each local script takes about one second on this machine.

## Scope and normalization

[Numerical] The transfer calculation covers XX, XXZ Delta = 0.5, Heisenberg, and XXZ Delta = 2, at beta = 0.2, 2, 8 and m = 1, 2, 3, 4. There are 48 full transfer calculations, maximum column dimension 256. Each has scalar trace comparisons at L = 4, 6, 16, 64. The two-site transfer U is positive; crossed gaps are minus one half the logarithm of normalized eigenvalues. Rotation is translation by two time-history entries, giving momentum modulo m. audit.py folds the bosonic reference momentum modulo the same m before comparison.

[Numerical] The half-spreads M are 1/2 for Delta = 0, 0.5, 1, and 3/4 for Delta = 2. The tested reference has S = 6, V' = 4eM. Three constants are compared, K = 0, K = log 3, and Kstar = log 3+8/(3e) = 2.07962413179. The last is proved for exact XX in PROOF.md. Its use for the interacting samples is a numerical test, not a theorem.

## Transfer verification and errors

| Test | Numerical outcome |
| :--- | :--- |
| Original-space versus spatial L = 4 Trotter contraction | Maximum relative error 4.03e-12. |
| Rotation commutator relative norm | Maximum 1.57e-10. |
| Minimum eigenvalue divided by top eigenvalue | -8.42e-17, consistent with roundoff. |
| log Ztilde <= R <= 2 log Ztilde | All 192 samples pass. |
| Scalar boson bound with Kstar | All 192 samples pass; minimum margin 0.7280193312. |
| Scalar boson bound with K = log 3 | Two coarse m = 1 samples fail. |
| Momentum-sector boson bound with Kstar | 56 of 480 folded-sector comparisons fail. |

[Numerical] At beta = 2, m = 4, the errors log Z_m(4)-log Z(4) and the number of distinct energies detected by one nonzero-frequency Sz source are

| Model | Trotter log-trace error | Distinct detected energies |
| :--- | ---: | ---: |
| XX | 0.00587127 | 3 |
| XXZ Delta = 0.5 | 0.01782895 | 8 |
| Heisenberg | 0.03950913 | 8 |
| XXZ Delta = 2 | 0.09952505 | 10 |

[Numerical] At beta = 8, m = 4 the corresponding Trotter log-trace errors are 0.75537, 1.75153, 2.85461, and 4.77674. These are substantial. The low-temperature finite-step spectra are not presented as converged quantum-transfer data. The exact XX checks below avoid that issue.

[Numerical] The source calculation uses the similarity-transformed diagonal insertion invsqrt(Ke) D_k sqrt(Ke), not an assumed Hermitian insertion. It counts distinct resolved energies having relative squared overlap above 1e-8 and eigenvalue above 1e-13 of the top eigenvalue. Only beta >= 2, m >= 3 is used for this diagnostic. Full source weights and gaps are in checks.json.

## Failure witnesses

[Numerical] The K = log 3 scalar candidate fails at L = 64, beta = 8, m = 1. At Delta = 2, log Ztilde = 1.36426917124 and the margin is -0.25299251189. At Heisenberg, log Ztilde = 1.32260790778 and the margin is -0.22341644914. Raising K to the explicit XX value Kstar makes these samples pass. These finite-step failures do not refute a quantum-chain bound with adjustable K.

[Numerical] The largest folded-momentum failure occurs at Heisenberg, beta = 0.2, m = 4, L = 64, k = 2. The actual sector log trace is -350.4816546; the reference sector log trace is -736.7751321. Even after adding Kstar the margin is -384.2138534. A scalar comparison can pass while a momentum-sector comparison fails. The exact XX argument in PROOF.md supplies an analytic lattice failure of the global linear momentum-preserving comparison.

[Numerical] Testing the low-frequency cone on every finite-step eigenstate gives a negative minimum margin -0.0724667252 at Heisenberg, beta = 8, m = 2, k = 1. This is not a counterexample to the conditional continuum detected-sector statement. That statement has neither a finite-Trotter guarantee nor a guarantee for all undetected sectors.

## Exact XX and auxiliary bounds

[Numerical] The exact XX spin trace is evaluated with both Jordan-Wigner parity sectors. The independent crossed product uses 20000 modes. Thirty samples cover L = 4, 6, 16, 64, 256, 1024 and beta/L = 0.01, 0.1, 1, 4, 32. The largest absolute full log-trace discrepancy is 5.59e-9, at a large full log trace. The analytic integral bound, its simpler 8 beta/(3 pi L) consequence, and the Kstar bosonic majorant pass all samples.

[Numerical] Four additional samples at L = 2 use the physical trace to evaluate log Ztilde without truncation. The product convergence is slower there. At beta = 32 the 20000-mode product is short by 0.00129689492, while its analytic tail allowance is 0.00129694357. This discrepancy is not reported as roundoff or used as the trace in the inequality check. Every L = 2 trace discrepancy lies within its explicit tail allowance and numerical tolerance.

[Proved error bound used in the computation] The omitted logarithm of each XX sector product is bounded by

    4/(L-1) [beta/(4 pi)]^L (N-1/2)^(1-L).

This follows from asinh z >= log(2z), log(1+y) <= y, and an integral bound on the decreasing tail. The normalized total trace obeys the same logarithmic allowance because its sector contributions are positive. For the parity-projected odd sector, adding omitted fermion modes multiplies each old parity sum by at most the unrestricted product tail, since the old even sum is at least the old odd sum. This justifies the bound despite the difference product in its formula.

[Numerical] The fixed k = 1 pair family in XX has 4, 20, and 80 distinct energies below 0.5 at beta = 100, 500, and 2000. Fifteen exact local branch samples check the nonlinear rate (all pass) and the global linear rate (fails at high frequency). Twenty-eight exact-density bubble samples check the full C125 low-frequency correlator estimate with ||Sz|| = 1/2; all pass. The smallest absolute bound margin in that sample is 3.64e-5. The bubble is summed over 60001 fermionic frequencies and is a numerical check rather than an interval enclosure.

[Numerical] Fifty values of t between 0.001 and 1000 test the boson-product upper bound pi^2 t/6 and the lower bound b t/2 for t >= 1. All pass. Twenty samples test the softened high-frequency reference bound with c_alpha = 3.38828937933; the omitted power-law tail is explicitly included, and all pass. Twenty-four Gaussian-charge samples cover 1, 2, 10, and 100 charge directions at fixed total stiffness 4, and beta/L from 0.001 to 100. All obey the 16 beta/L estimate.

[Numerical] The diagonal ferromagnetic witness is checked at d = 2, 3, 5 and beta J = 4, 8, 16, 24, with even L near exp(beta J/2). For d = 5 at beta J = 24, L = 162754, log Ztilde = 1.60941333577 compared with log 5 = 1.60943791243; the oscillator contribution is unresolved at this aspect ratio. This agrees with the proved necessity K(d) >= log d.

## Reproducibility and unfinished computation

[Assumed work record] A rerun against the parent-folder checks.json failed with KeyError: transfer after checks.py and its output were replaced concurrently. The parent files were preserved, and this subfolder was created to give the present run a stable implementation and schema. Both isolated scripts subsequently passed. The two runs use different Trotter grids and constants; their numerical summaries must not be mixed.

[Assumed scope limit] No large interacting computation was run locally. run_typhon.slurm prepares m up to 6, guarded by a SLURM job identifier. It has not been submitted. No Typhon connection, network command, or separate-agent review was performed. Floating-point tests are not certified spectral enclosures or proofs of continuum inequalities.
