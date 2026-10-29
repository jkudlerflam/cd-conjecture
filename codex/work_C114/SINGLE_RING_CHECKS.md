# C114 additional single-ring checks

[Numerical] The coherent cost is evaluated using the proved formula q_L/S_L. This does not diagonalize the doubled ring. All single-ring Hilbert spaces have dimension at most 729.

| model | L | q_L | S_L | q_L/S_L | L q_L/(M S_L) |
|:--|--:|--:|--:|--:|--:|
| XX | 4 | 0.707106781 | 2.562500000 | 0.275944110 | 2.207552878 |
| XX | 6 | 0.666666667 | 3.489197531 | 0.191065900 | 2.292790801 |
| XX | 8 | 0.653281482 | 4.338423436 | 0.150580388 | 2.409286201 |
| XXZ_minus | 4 | 0.593070331 | 2.535058475 | 0.233947397 | 1.871579174 |
| XXZ_minus | 6 | 0.568252661 | 3.423359504 | 0.165992692 | 1.991912308 |
| XXZ_minus | 8 | 0.559769553 | 4.235635635 | 0.132157155 | 2.114514472 |
| XXZ_plus | 4 | 0.843070331 | 2.578577889 | 0.326951664 | 2.615613310 |
| XXZ_plus | 6 | 0.789344663 | 3.518685922 | 0.224329389 | 2.691952668 |
| XXZ_plus | 8 | 0.771748589 | 4.380109580 | 0.176193900 | 2.819102398 |
| Heisenberg | 4 | 1.000000000 | 2.583333333 | 0.387096774 | 3.096774194 |
| Heisenberg | 6 | 0.934258546 | 3.526685184 | 0.264911240 | 3.178934882 |
| Heisenberg | 8 | 0.912773352 | 4.390519431 | 0.207896438 | 3.326343014 |
| SU3 | 3 | 2.666666667 | 1.666666667 | 1.600000000 | 4.800000000 |
| SU3 | 6 | 2.201850425 | 2.509417731 | 0.877434792 | 5.264608752 |

[Numerical] Maximum eigenpair residual 1.33e-13. All dimension and norm inequalities passed. These sizes do not establish asymptotic powers for the interacting models.
