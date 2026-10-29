# Two-seam response check

[Numerical] Two three-site XX rings, six qubits total. These odd rings test the finite-temperature matrix identity, not a vacuum Casimir coefficient. A 24 by 24 Gauss-Legendre quadrature integrates the Kubo susceptibility.

| beta | F(2L)-2F(L) | half integrated chi | error | min sampled chi | max sampled chi |
|--:|--:|--:|--:|--:|--:|
| 2 | -0.182863659037 | -0.182863659037 | 1.7e-16 | -0.398973119 | -0.263645645 |
| 8 | -0.677778013882 | -0.677778013882 | 7.8e-16 | -2.165430123 | -0.089115365 |

[Numerical outcome] All assertions passed. Susceptibility need not be positive on non-vacuum-compatible finite rings.
