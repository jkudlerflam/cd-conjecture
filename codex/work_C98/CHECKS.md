# C98 checks

[Numerical] Small-ring ED, exact-formula evaluations and finite-matrix identities. No large-size interacting extrapolation or new Bethe-ansatz computation.

## XXZ ring calibration

h=SxSx+SySy+Delta SzSz, M=1/2, d=2. Q_L=6L(e0 L-E_L)/(pi^2 M).

| Delta | L | E_L | Q_L | residual |
|:--|--:|--:|--:|--:|
| -0.5 | 4 | -1.186140661635 | 0.428294128 | 1.2e-15 |
| -0.5 | 6 | -1.704757983078 | 0.420517762 | 3.9e-15 |
| -0.5 | 8 | -2.239078212171 | 0.417532013 | 1.6e-14 |
| -0.5 | 10 | -2.779413246380 | 0.416098355 | 5.5e-13 |
| -0.5 | 12 | -3.322693494020 | 0.415309436 | 4.6e-14 |
| 0 | 4 | -1.414213562373 | 0.685615408 | 1.2e-15 |
| 0 | 6 | -2.000000000000 | 0.657587569 | 3.3e-15 |
| 0 | 8 | -2.613125929753 | 0.648262727 | 3.6e-14 |
| 0 | 10 | -3.236067977500 | 0.644027219 | 1.8e-12 |
| 0 | 12 | -3.863703305156 | 0.641747365 | 5.4e-12 |
| 0.5 | 4 | -1.686140661635 | 0.905279624 | 8.5e-16 |
| 0.5 | 6 | -2.368033988750 | 0.861072728 | 3.2e-15 |
| 0.5 | 8 | -3.086994355137 | 0.846179619 | 4.3e-14 |
| 0.5 | 10 | -3.819032780516 | 0.839337964 | 2.4e-14 |
| 0.5 | 12 | -4.557272440830 | 0.835619255 | 7.3e-13 |
| 1 | 4 | -2.000000000000 | 1.105995832 | 9.6e-16 |
| 1 | 6 | -2.802775637732 | 1.049714203 | 1.6e-14 |
| 1 | 8 | -3.651093408937 | 1.030226965 | 2.9e-13 |
| 1 | 10 | -4.515446354492 | 1.021008083 | 4e-14 |
| 1 | 12 | -5.387390917445 | 1.015842550 | 4e-13 |

## Finite sharp-bound failures

Heisenberg L=4 deficit 0.227411277760, proposed no-remainder ceiling 0.205616758356.
Heisenberg E(8)-2E(4)=0.348906591063, proposed no-remainder ceiling 0.308425137534.

## XX open-chain and rank checks

| L | E_open exact | ED error | half-chain Schmidt rank |
|--:|--:|--:|--:|
| 4 | -1.118033988750 | 6.7e-16 | 4 |
| 6 | -1.746979603717 | 2.2e-15 | 8 |
| 8 | -2.379385241572 | 4.9e-15 | 16 |
| 10 | -3.013337091666 | 6.2e-15 | 32 |
| 12 | -3.648114905279 | 2.7e-15 | 64 |

Open-minus-ring limiting cost = 1/2-1/pi = 0.181690113816.
Product-of-two-ring seam cost tends to 2/pi = 0.636619772368.

## Fundamental SU(N) ring calibration

h=P, M=1. Vacuum-compatible sizes only. Q_L uses the known Bethe bulk energy.

| N | L | E_L | Q_L | Q_infinity |
|--:|--:|--:|--:|--:|
| 2 | 4 | -2.000000000000 | 1.105995832 | 1.000000000 |
| 2 | 6 | -2.605551275464 | 1.049714203 | 1.000000000 |
| 2 | 8 | -3.302186817874 | 1.030226965 | 1.000000000 |
| 3 | 3 | -3.000000000000 | 1.024521904 | 0.841239671 |
| 3 | 6 | -4.605551275464 | 0.888964983 | 0.841239671 |
| 3 | 9 | -6.579435759086 | 0.864827875 | 0.841239671 |
| 4 | 4 | -4.000000000000 | 0.850519516 | 0.750000000 |
| 4 | 8 | -6.921009740038 | 0.778288377 | 0.750000000 |

## Symmetric SU(N) inherited continuum formulas

Formula checks only. The all-N,k proof and physical inputs are inherited from C14/C35.

| N | k | log_2 d | A/(2pi M) | A/(pi M log_2 d) |
|--:|--:|--:|--:|--:|
| 2 | 1 | 1.000000000 | 0.500000000 | 1.000000000 |
| 2 | 2 | 1.584962501 | 0.500000000 | 0.630929754 |
| 3 | 1 | 1.584962501 | 0.666666667 | 0.841239671 |
| 3 | 2 | 2.584962501 | 0.711111111 | 0.550190659 |
| 7 | 2 | 4.807354922 | 1.015873016 | 0.422632833 |
| 8 | 2 | 5.169925001 | 1.050000000 | 0.406195448 |
| 10 | 10 | 16.495261691 | 1.690014903 | 0.204909135 |
| 100 | 100 | 194.850520479 | 9.637817921 | 0.098925247 |

## Uniform partial-copy swap

Four-site density-matrix checks use the thermodynamic nearest-neighbor marginal.

| model | theta | measured cost per doubled bond | formula | error |
|:--|--:|--:|--:|--:|
| Heisenberg | 0.000000000 | 0.000000000000 | 0.000000000000 | 0 |
| Heisenberg | 0.392699082 | 0.352493206039 | 0.352493206039 | 5.6e-17 |
| Heisenberg | 0.785398163 | 0.704986412078 | 0.704986412078 | 0 |
| Heisenberg | 1.570796327 | 0.000000000000 | 0.000000000000 | 1.1e-32 |
| XX | 0.000000000 | 0.000000000000 | 0.000000000000 | 0 |
| XX | 0.392699082 | 0.223658011958 | 0.223658011958 | 1.1e-16 |
| XX | 0.785398163 | 0.447316023917 | 0.447316023917 | 0 |
| XX | 1.570796327 | 0.000000000000 | 0.000000000000 | 6.7e-33 |

## Dyadic identity and concavity obstruction

xi=1: dyadic identity error 0; deficit <=1, Casimir-like coefficient K=1.
xi=10: dyadic identity error 0; deficit <=1, Casimir-like coefficient K=10.
xi=1000: dyadic identity error 0; deficit <=1, Casimir-like coefficient K=1000.

Maximum reported ring eigenpair residual 5.38e-12.

[Numerical outcome] All assertions passed. The two no-remainder Heisenberg inequalities are violated as expected. The asymptotic logarithmic conjecture was not falsified.
