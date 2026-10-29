# C103 checks

[Numerical] Formula and normalization checks, not a general proof.

## Blocked free chain

A = 0.466666666667; predicted L^2 correction = -0.144997891818.

| L | A_ring(L) | L^2(A_ring-A) |
| :-- | --: | --: |
| 8 | 0.464344037743 | -0.148648251100 |
| 12 | 0.465648573607 | -0.146605400651 |
| 16 | 0.466096747830 | -0.145899222085 |
| 24 | 0.466414240331 | -0.145397569522 |
| 32 | 0.466524847789 | -0.145222530797 |
| 48 | 0.466603690245 | -0.145097675182 |
| 64 | 0.466631253088 | -0.145054016243 |
| 96 | 0.466650930676 | -0.145022890408 |

Ring energy residual <= 5.68e-14; blocked band residual <= 1.39e-15.
Gauge-optimized M lies in [1.119803902719, 2.108271469549].
Therefore Q <= 0.044217448309; the example obeys the free theorem.

## XX thermal current

| T/t | R_J(T)/t | (R_J-2t)t/T^2 |
| :-- | --: | --: |
| 0.1 | 1.964725324301 | -3.527467569871 |
| 0.05 | 1.991321195735 | -3.471521706007 |
| 0.025 | 1.997838379335 | -3.458593063204 |
| 0.0125 | 1.999460091267 | -3.455415888797 |

Predicted coefficient = -3.454361540381; independent k integral residual = 2.22e-16.
R_J(10t)/t = 0.000475087413; high-temperature leading value = 0.000477464829.

## Evaluations of inherited C101 formulas

[Numerical] These evaluate established formulas and do not repeat Bethe ansatz or many-body simulations.

| Family | Q |
| :-- | --: |
| XX | 0.636619772368 |
| Heisenberg | 1.000000000000 |
| SU(3), symmetric k=1 | 0.841239671429 |
| SU(4), symmetric k=1 | 0.750000000000 |
| SU(2), symmetric k=2 | 0.630929753571 |
| SU(3), symmetric k=2 | 0.550190659178 |
| SU(8), symmetric k=2 | 0.406195447596 |
| Sp(4), assumed sectors | 0.666666666667 |
| Sp(6), assumed sectors | 0.612516944788 |
| Sp(8), assumed sectors | 0.566666666667 |
| Sp(16), assumed sectors | 0.465277777778 |
| Hubbard repulsive, upper bound | 0.477464829276 |
| Hubbard attractive, upper bound | 0.159154943092 |
| t-J q=2, upper bound | 0.946394630357 |
| t-J q=3, upper bound | 1.000000000000 |

## Optical-moment identity on small XXZ rings

[Numerical] Derivatives are taken on the finite-ring ground branch before any size limit.

| L | Delta | K | D | R | 2s0^2/s1 |
| :-- | :-- | --: | --: | --: | --: |
| 6 | 0.3 | 0.3312673390 | 0.3309958267 | 0.0002715123 | 0.0002715123 |
| 6 | 1.0 | 0.3114195153 | 0.3066247547 | 0.0047947606 | 0.0047947606 |
| 8 | 1.0 | 0.3042577841 | 0.2961924395 | 0.0080653446 | 0.0080200399 |

Double-commutator moment residual <= 4.86e-17; finite-difference stiffness residual <= 5.23e-08.

[Numerical] All assertions passed. No lattice counterexample to the proposed A bound was found.
