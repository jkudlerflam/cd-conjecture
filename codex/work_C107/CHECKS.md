# C107 checks

[Numerical] Small spectral and quadrature checks of the analytic derivations. Run with `python3 checks.py`. These are not extrapolation evidence for the general conjecture.

[Numerical] Nine full spin matrices, L = 4, 6, 8 and gamma = 0.25, 0.5, 1, agree with the antiperiodic fermion sum. Largest dimension 256, largest energy discrepancy 1.95e-14.

| gamma | L | L^2 (Q_L - gamma) | predicted coefficient |
|:--|--:|--:|--:|
| 0.1 | 64 | -8.846746963 | -8.520758466 |
| 0.1 | 128 | -8.591680669 | -8.520758466 |
| 0.1 | 256 | -8.538040401 | -8.520758466 |
| 0.25 | 64 | -3.184327373 | -3.166498079 |
| 0.25 | 128 | -3.170889109 | -3.166498079 |
| 0.25 | 256 | -3.167591627 | -3.166498079 |
| 0.5 | 64 | -1.153509695 | -1.151453847 |
| 0.5 | 128 | -1.151966438 | -1.151453847 |
| 0.5 | 256 | -1.151583304 | -1.151453847 |
| 1 | 64 | 0.287881747 | 0.287863462 |
| 1 | 128 | 0.287868124 | 0.287863462 |
| 1 | 256 | 0.287866564 | 0.287863462 |

[Numerical] The sign and large-L coefficient agree with the Euler-Maclaurin derivation for all four anisotropies. The smallest anisotropy has the longest crossover.

| L | half-filled ferromagnetic XX Q_L, t = 1 |
|--:|--:|
| 63 | 0.500036266 |
| 64 | 2.000562376 |
| 127 | 0.500008924 |
| 128 | 2.000140567 |

[Numerical] Even and odd XX lengths approach 2 and 1/2 respectively, in agreement with the exact formulas.

[Numerical] The abstract countermodel derivative and additivity checks pass on the stated grids. Their proofs do not rely on sampling.

[Proved scope] No large diagonalization, finite-size fit, or cluster job was used. No physical bond family with unbounded A was constructed.
