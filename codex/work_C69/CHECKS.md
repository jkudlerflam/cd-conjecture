# C69 elementary checks

[Numerical] All assertions in checks.py passed. These are formula checks, not a microscopic CFT simulation.

Maximum centered-difference error for the strip entropy derivative: 4.30618e-10.
Multisector late-time entropy rate: 0.51611879309.
Finite-time lower fraction at q=4: 0.889682199924; the integrated strip rate satisfies the bound.
Weighted harmonic speed: 1.1724137931, between 0.2 and 4.
Maximum sampled exponential-tail / LR-majorant ratio: 0.31711.

| Constant | Value |
|:--|--:|
| SIE conditional coefficient 216/pi | 68.7549354157 |
| Thermal coefficient 3/pi | 0.954929658551 |
| Thermal coefficient at alpha=0.9 | 1.06103295395 |
| Literal pi/6 coefficient 6/pi | 1.9098593171 |
| Spatial doubling coefficient 3/log 2 | 4.32808512267 |
| Heisenberg-required coefficient pi/log 2 | 4.53236014183 |
| Conservative LR coefficient 6e | 16.3096909708 |

[Proved in RESULT.md, tested here] Equal spatial scaling leaves entropy increments invariant and divides rates by the scale. The pi/6 convention with u=v/2 agrees with the ordinary nonchiral pi/3 formula. The one-site spatial bound is at most 3(n+1) log d/alpha.

[Assumed/open] These checks do not verify uniform lattice preparation, dynamical CFT convergence, or the conjectured universal bound on A.
