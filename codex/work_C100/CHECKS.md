# C100 checks

[Numerical checks, not proofs of CFT matching.] Four seeded random Hermitian
two-qubit bonds, each centered and normalized to operator norm one, on six-site
rings. Temperatures sampled via beta = 0, 0.3, 2, 8. Three nonzero momenta.
The centered norm is used, not a numerical gauge optimization for M.

- Maximum continuity-equation matrix residual: 1.691e-15.
- Maximum first-moment discrepancy, double commutator versus spectral sum: 1.776e-15.
- Smallest first moment: -3.123e-17 (roundoff allowed at beta = 0).
- Largest ratio m1/(8 m^3 q^2): 0.050875, with m = 1.
- Largest local midpoint correlator divided by 4 m^4: 0.181767.
- CFT midpoint spatial integral error: 0.000e+00.
- Largest CFT midpoint Fourier-transform error: 5.551e-17.

All assertions passed. No large computation, infrared fit, or Drude-limit
extrapolation was performed. Analytic checks of scaling, chiral conventions,
and the absence of a temperature-independent conclusion are in PROOF.md.
