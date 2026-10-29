# C106 checks

[Numerical] Small matrix identities and one-dimensional quadratures only. These check the analytic proof; no infrared extrapolation of an interacting chain is used.

[Numerical] delta=0.1: bond half-spread error 0, extremal-state equal-marginal residual 1.42e-16, v/M=0.25.
[Numerical] delta=0.01: bond half-spread error 4.44e-16, extremal-state equal-marginal residual 2.48e-16, v/M=0.25.
[Numerical] delta=0.001: bond half-spread error 0, extremal-state equal-marginal residual 2.56e-16, v/M=0.25.

| delta | T/delta | s(T) / [2 pi T/(3v)] |
| :-- | --: | --: |
| 0.01 | 0.0001 | 1.0000155457 |
| 0.01 | 0.001 | 1.0015748618 |
| 0.01 | 0.01 | 1.0971395071 |
| 0.01 | 1 | 0.2269454517 |
| 0.01 | 10 | 0.2511222839 |
| 0.001 | 0.0001 | 1.0000155466 |
| 0.001 | 0.001 | 1.0015749534 |
| 0.001 | 0.01 | 1.0971378707 |
| 0.001 | 1 | 0.2269265456 |
| 0.001 | 10 | 0.2496366164 |
| 0.0001 | 0.0001 | 1.0000155466 |
| 0.0001 | 0.001 | 1.0015749543 |
| 0.0001 | 0.01 | 1.0971378543 |
| 0.0001 | 1 | 0.2269263566 |
| 0.0001 | 10 | 0.2496223136 |

[Proved check] A/M = 2r = 1/2 in this family, so the free bound A <= 2M is respected.

[Numerical] All assertions passed.
