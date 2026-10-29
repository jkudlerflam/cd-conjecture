# C114 XX coherent reconnection

[Numerical] Exact finite-ring correlation matrices, of dimension at most 32. The energy q_L / sum_r tr(rho_r^2) is the optimum in the coherent prefix-swap subspace, as proved in PROOF.md. No many-body extrapolation.

| L | purity sum | optimum cost | L times cost | L^(3/4) times cost |
|--:|--:|--:|--:|--:|
| 4 | 2.562500000 | 0.275944110 | 1.103776439 | 0.780487805 |
| 8 | 4.338423436 | 0.150580388 | 1.204643101 | 0.716285073 |
| 16 | 7.316960282 | 0.087567629 | 1.401082061 | 0.700541030 |
| 32 | 12.321928365 | 0.051748684 | 1.655957888 | 0.696244526 |
| 64 | 20.737246578 | 0.030711672 | 1.965546985 | 0.694925801 |

[Numerical] Three arbitrary complex Hermitian bonds at d=2, L=4 also pass the Gram identity, Re K=q I, and the uniform coherent upper bound. No reality or on-site symmetry was imposed for these tests.

[Proved] For XX, the number-variance inequality gives tr(rho_r^2) <= min(r,L-r)^(-2/pi^2) for 1 <= r < L. It proves that L times the optimized cost diverges. The observed L^(-3/4) scaling is a sharper CFT calibration, not required for that proof.

[Numerical] All assertions passed.
