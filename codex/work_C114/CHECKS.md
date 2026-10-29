# C114 ring checks

[Numerical] Full spin-ring diagonalization with sparse matrices. Continuum A values are inherited calibration inputs, not fitted results.

| model | L | E(L) | E(2L) | D_L | L D_L/M | L coherent/M | bare cost |
|:--|--:|--:|--:|--:|--:|--:|--:|
| XX | 4 | -1.414213562373 | -2.613125929753 | 0.215301195 | 1.722409560 | 2.207552878 | 0.707106781 |
| XXZ_minus | 4 | -1.186140661635 | -2.239078212171 | 0.133203111 | 1.065624889 | 1.871579174 | 0.593070331 |
| XXZ_plus | 4 | -1.686140661635 | -3.086994355137 | 0.285286968 | 2.282295745 | 2.615613310 | 0.843070331 |
| Heisenberg | 4 | -2.000000000000 | -3.651093408937 | 0.348906591 | 2.791252729 | 3.096774194 | 1.000000000 |
| SU3 | 3 | -3.000000000000 | -4.605551275464 | 1.394448725 | 4.183346174 | 4.800000000 | 2.666666667 |

[Numerical] Every positive D_L refutes the zero-cost lift and the zero-error comparison E(L)/L >= E(2L)/(2L). The inequalities D_L <= optimized coherent cost <= bare cost <= 4M pass. Using the bulk energies inherited from C98, every tested E(L) is also strictly below L e0, excluding an exact finite-range positive certificate that wraps on that size. No approximate SOS certificate has been constructed or numerically tested.

[Numerical] Maximum eigenpair residual 1.33e-13. Largest Hilbert-space dimension 729.

[Numerical] All assertions passed. The Gram matrix agrees with single-ring interval purities. No asymptotic bound follows from these finite sizes.
