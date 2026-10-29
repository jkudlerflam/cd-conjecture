# C94 numerical checks

[Numerical] These are finite logical vacuum-clock calculations, with quantum contact sources and distributed marker energies. Full C88 auxiliary traces are treated analytically in PROOF.md, not diagonalized here.

Canonical comparisons: 50. Thermal comparisons: 10. Maximum joint matrix dimensions: 13365 canonical and 880 thermal.

Largest canonical eigenpair residual: 1.2836e-11. Largest bath linear-solve residual: 4.2799e-13. All canonical defect costs are positive.

The table reports max(0,-Delta R)/(J W). The canonical case is (30) against (15,15), N=30, v=1. The thermal case is (2,2) against (4), N=4, s/t=0.65, T=J.

| J/t | Canonical observed epsilon | Thermal observed epsilon |
| :-- | --: | --: |
| 0.04 | 3.15548577e-05 | 0.000468636879 |
| 0.02 | 1.06565203e-05 | 0.00011827144 |
| 0.01 | 4.13336142e-06 | 3.83056026e-05 |
| 0.005 | 1.78199823e-06 | 1.89019317e-05 |
| 0.0025 | 8.21654969e-07 | 9.47357197e-06 |

Largest thermal residual-difference cutoff error bound: 9.1005e-10.

[Numerical] At the separate free vacancy crossing, the logical response grows as 1/T.

| T/t | D_wrong-D_good, t=1 | T(D_wrong-D_good) |
| :-- | --: | --: |
| 0.01 | 0.5114280838 | 0.005114280838 |
| 0.003 | 1.551516888 | 0.004654550663 |
| 0.001 | 4.523203005 | 0.004523203005 |
| 0.0003 | 14.92410442 | 0.004477231325 |
| 0.0001 | 44.64096559 | 0.004464096559 |

Analytic limit for the last column: 0.00445752917643.

[Numerical] Half-density coefficient checks use the exact bath resolvent. They do not form the larger joint interacting matrices at N=9.

| N=v | Defect | W | (B_wrong-B_good)/W |
| :-- | :-- | --: | --: |
| 3 | one_pin | 1 | -0.0932690230279 |
| 6 | one_pin | 1 | -0.0701028652189 |
| 6 | length_pair | 6 | -0.00307336964928 |
| 6 | long_cell | 6 | -0.00595479366352 |
| 9 | one_pin | 1 | -0.0616632510684 |
| 9 | length_pair | 6 | -0.00304738869758 |
| 9 | long_cell | 6 | -0.00572515213572 |

[Scope] The data support decreasing normalized errors on these finite systems. They do not bound a thermodynamic supremum or an all-order remainder uniform in defect length. The b=15 joint calculations have high occupancy. No half-density b=15 interacting trace is claimed.

[Reproduction] Run python3 -B checks.py, python3 -B coefficient_scaling.py, python3 -B crossing_certificate.py, then python3 -B summarize.py from this folder. The inherited modules are local copies of the read C89 implementations, with the local import adjusted.
