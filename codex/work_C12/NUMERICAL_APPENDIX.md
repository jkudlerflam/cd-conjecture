## Final numerical appendix

At L=24, all momentum sectors k=0,...,12 were checked. Reality of h gives the remaining sectors by complex conjugation. The lowest state is in k=0, with E=-14.18107312747. The largest reported eigenvector residual in this check is 6.2e-09.

The completed periodic DMRG checkpoints are listed below. These are raw finite-bond-dimension fits, not claims for the thermodynamic central charge.

| L | Largest completed chi | Energy | Half-chain S | Half-chain S2 | Shape fit, vN | Shape fit, Renyi2 | Density range |
|---:|---:|---:|---:|---:|---:|---:|---|
| 24 | 256 | -14.18095565151 | 2.454555 | 1.864485 | 2.9110 | 2.7302 | 0.3295 to 0.3359 |
| 36 | 192 | -21.27165110127 | 1.680674 | 1.014356 | 1.2272 | 0.6464 | 0.2710 to 0.4172 |
| 48 | 192 | -28.35843670102 | 1.827161 | 1.335844 | 0.8680 | 0.4691 | 0.2387 to 0.4617 |

At L=24, the last checkpoint has energy error 0.000117 and half-chain entropy deficit 0.0101 relative to ED. Density modulation and continuing bond-dimension dependence prevent a reliable large-L extrapolation. The interrupted runs and wall-clock stopping condition are recorded in PROGRESS.md.

The nearby detuning delta=0.75 is also affected by momentum and size effects.
At L=21, the lowest tested momentum is k=0, the shape estimates are 3.6497 and 3.7907, and half-chain S=2.660515.
At L=24, the lowest tested momentum is k=8, the shape estimates are 2.4480 and 1.9032, and half-chain S=2.214717.
All sectors were checked at L=21; only k=0,1,4,8,12 were checked at L=24 in this auxiliary detuning scan. Its molecular minima have already disappeared, so this does not rescue the original three-fluid argument.

Numerical environment: Python 3.12.3, NumPy 2.4.6, SciPy 1.17.1. BLAS was limited to one thread per Python process. All random scans use fixed seeds in their scripts.
