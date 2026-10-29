# C110 checks

[Numerical] All matrices have dimension at most 256. These checks test the explicit derivations, not the general bound.

[Numerical] The interacting gapped bond has trace zero and equal partial traces, so it lies in the Hilbert-Schmidt canonical gauge.

| L | uniform energy | gap | rigorous gap lower bound |
| :-- | --: | --: | --: |
| 4 | -4.0000000000 | 1.7021749414 | 1.2000000000 |
| 6 | -6.0000000000 | 1.7017318080 | 1.2000000000 |
| 8 | -8.0000000000 | 1.7017317653 | 1.2000000000 |

| L | L times exact weighted deficit | L times trial deficit | rigorous upper bound times L |
| :-- | --: | --: | --: |
| 4 | 0.0556592398 | 0.0535568808 | 0.0888888889 |
| 6 | 0.0607816130 | 0.0595840690 | 0.0991648352 |
| 8 | 0.0629728622 | 0.0622321813 | 0.1036475762 |

[Proved formula, numerical evaluation] For f(x)=1+0.4 cos(2 pi x), the integral of f prime squared divided by f is 3.295850206548. The rigorous asymptotic lower and upper coefficients are 0.065917004131 and 0.109861673552.

[Numerical] Scalar sums at L=64,128,256,512 verify the product trial coefficient and the exactly solvable transverse-coboundary coefficient. No large spin matrices were formed.

| L | flux | XX exact energy | finite-sum energy |
| :-- | --: | --: | --: |
| 4 | 0.00000000 | -2.8284271247 | -2.8284271247 |
| 4 | 1.04719755 | -2.7320508076 | -2.7320508076 |
| 4 | 3.14159265 | -2.0000000000 | -2.0000000000 |
| 6 | 0.00000000 | -4.0000000000 | -4.0000000000 |
| 6 | 1.04719755 | -3.9392310120 | -3.9392310120 |
| 6 | 3.14159265 | -3.4641016151 | -3.4641016151 |
| 8 | 0.00000000 | -5.2262518595 | -5.2262518595 |
| 8 | 1.04719755 | -5.1815405504 | -5.1815405504 |
| 8 | 3.14159265 | -4.8284271247 | -4.8284271247 |

[Numerical] Nine XX flux matrices agree with the even-length formula. At flux pi, Q_L approaches -4 for hopping t=1.

[Numerical] The one-site twist optimization agrees with max(1,2|r|), including the gauge-dependent linear divergence at r=1.

[Work record] All checks passed. No heavy numerics, network shell commands, or cluster jobs were used.
