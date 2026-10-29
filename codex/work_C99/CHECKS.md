# C99 checks

[Numerical] All assertions passed in python3 checks.py. These are small matrix and quadrature checks of the analytic proofs, not evidence for or against the general Casimir bound.

Four-spin active-chain entropy derivative: 0.951426150714. Analytic value: 0.951426150896. Absolute error: 1.83e-10.
Periodic four-spin XX ground energy agrees with -2J/sin(pi/L) to 1.33e-15. The asymptotic Casimir expansion is analytic, not fitted.
Coherent Ising spectator: maximum correlation error 1.22e-15; maximum entropy change under onsite real-time evolution 2.29e-16.
Spectator correlation lengths at k=1,2,4,6: 3.67186093, 27.2960221, 1490.47894, 81377.3957.
Euclidean one-qubit filtering: maximum derivative error 1.66e-10.
Seven-spin uniform chain: maximum distant-MI identity error 1.3e-15; rate at t=pi/8 is 1.24645048007, predicted 1.24645048028.
Multisector local-join positive integral: recovered A=5.815, expected 5.815, error 8.88e-16. Finite-window arctangent identity also passed.

[Proved analytically, not established by finite sampling] Divergence of xi_Z r_0 and of the Euclidean rate M log N follows from their displayed formulas. Finite calculations check signs, factors, tensor assignments, and the cancellation of within-half terms.

[Not tested] No uniform continuum matching, general capacity lower bound, or microscopic violation of A<=C M log d has been established.
