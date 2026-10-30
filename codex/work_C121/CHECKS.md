# C121 numerical checks

2026-10-29.

[Numerical] Reproduce in this folder with

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 checks.py

The final run takes about two seconds. The complete output is checks.json and the compact output is summary.json. All dense Hamiltonians have dimension at most 729. Every small-ring trace uses the full spectrum. No charge-sector restriction, spectral truncation or eigenvalue floor is used.

## 1. Definitions and scope

[Assumed conventions] XX has h = Sx Sx + Sy Sy, d = 2, M = 1/2. Heisenberg has h = S dot S, d = 2, M = 1/2. SU(3) has h = P, d = 3, M = 1. The compared lengths are 4 and 8 for the qubit chains and 3 and 6 for SU(3). The measured seam basis is the standard Z or color basis.

[Numerical] There are 48 model, spacing and insertion-count cases. Each uses tau = a L/M and beta = N tau, with a = 0.25, 1, 2, 4 and N = 1, 2, 4, 8. Q_L pinches one site, while Q_D pinches the two sites 0 and L on the doubled ring at the same times. Both are computed by positive density-operator propagation. At N <= 4, every individual one-ring history is also enumerated independently. The maximum difference between the channel trace and summed positive histories is 1.12e-16. The independent coherent-history sum agrees with the unmeasured trace to 8.89e-16.

[Numerical] All 48 cases obey the proved bounds

    0 <= log Z_L-log Q_L <= N log d,
    0 <= log Z_(2L)-log Q_D <= 2N log d.

All enumerated cases also obey Q_L^2 <= d^N S_pin and R <= 3N log d+G. Here G = log(S_pin/Q_D). At N = 8 the code records only the lower witness G >= log(Q_L^2/Q_D)-N log d, rather than enumerating d^8 histories. A positive lower witness certifies a failure; a negative one does not establish G <= 0.

## 2. The proposed positivity step fails

[Numerical] The following one-insertion values use tau = L/M. The two G columns would have to be nonpositive for the respective no-cost gluing claims.

| Model | L | beta | R | log Z_L-log Q_L | G = log(S_pin/Q_D) | log(S_coh/Z_(2L)) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| XX | 4 | 8 | 1.507576 | 0.668660 | 0.762551 | 0.814429 |
| Heisenberg | 4 | 8 | 2.747529 | 0.668650 | 1.888216 | 2.054382 |
| SU(3) | 3 | 3 | 3.923285 | 1.052837 | 2.292663 | 2.824673 |

[Numerical] There are 37 certified violations of S_pin <= Q_D among the 48 cases, including positive lower witnesses where the history sum was not enumerated. These are failures of an auxiliary inequality, not failures of the proposed universal bound on R or A. Several tested smaller spacings do not violate it, and no conclusion about all sufficiently small universal spacings follows.

[Proved relation to stoquasticity] The alternating Z rotation makes even XX and Heisenberg rings stoquastic and commutes with every measured projector. All the listed quantities are unchanged by that rotation. The failed second step therefore persists with nonnegative Euclidean matrix elements. SU(3) is tested without a stoquastic assumption.

## 3. Large XX rings and a direct measurement calculation

[Proved formula, numerical implementation] PROOF Section 7 derives the one-pinch echo formula retaining both Jordan-Wigner parity sectors. Fifteen independent comparisons with full ED at L = 4, 6, 8 and beta = 0.01, 0.2, 1, 4, 16, plus the colder one-pinch XX cases above, give a maximum discrepancy 3.22e-15. This checks both the full partition trace and the one-pinch loss.

[Numerical] Thirty-five free-fermion cases use L = 4, 8, 32, 128, 512, 2048, 8192 and a = beta M/L = 0.25, 1, 2, 4, 16. R uses the full trace at both L and 2L, reaching length 16384. A ground-energy subtraction and the cancellation-free doubling defect from C120 avoid subtracting bulk free energies. All cases obey C118's proved XX bound. Representative largest-ring values are

| L | beta/L | R | One-pinch loss | log(S_pin/Z_(2L)) |
| ---: | ---: | ---: | ---: | ---: |
| 8192 | 0.5 | 0.162103 | 0.693146961 | -1.917338 |
| 8192 | 2 | 1.378852 | 0.693147175 | -0.700589 |
| 8192 | 4 | 3.066133 | 0.693147180 | 0.986692 |
| 8192 | 8 | 6.279471 | 0.693147181 | 4.200030 |
| 8192 | 32 | 25.132742 | 0.693147181 | 23.053300 |

[Proved interpretation] Q_D <= Z_(2L), so a positive last column disproves S_pin <= Q_D without computing the large doubled-ring channel. Since N = 1 and spin-flip symmetry gives q_0 = q_1, the last column is computed from the exact identity S_pin = Q_L^2/2. The measurements at large L are actual one-pinch free-fermion evaluations, not inferred from the dimension bound alone. No large-L multi-pinch or two-seam channel determinant was computed.

## 4. Continuous pinching

[Numerical] At fixed beta = 4, N = 1, 4, 16, 64, 256 tests approach the predicted limit Tr exp[-beta D(H)] for each small model. At N = 256 the observed loss and the exact limiting loss are

| Model | Loss at N = 256 | Limiting loss |
| :--- | ---: | ---: |
| XX | 1.665344 | 1.675508 |
| Heisenberg | 1.803147 | 1.816281 |
| SU(3) | 6.803988 | 6.864527 |

[Numerical scope] These differences are finite-N errors, not a failure of the limit. No convergence rate is inferred from this short sequence. The XX open-chain formula for the limiting loss is checked independently against ED.

[Proved and numerical] The XX seam energy cost is E_Zeno-E_L = 1/2+(1/2)tan(pi/(2L)). It tends to M = 1/2. At L = 8192 and beta = 2L, continuous pinching has log-trace loss 8192.263679, while one coarse pinching has loss 0.693147175. This compares two different operations and displays the cost of making the seam classical throughout time.

## 5. Algebraic controls and limitations

[Numerical] Eighteen further cases test the channel bounds on a general complex Hermitian matrix and a commuting diagonal matrix, each of dimension nine, at three temperatures and N = 1, 2, 4. Explicit positive-history sums agree with the direct channel trace. The commuting cases satisfy Q = Z. These checks use no locality or stoquasticity assumption and calibrate the general finite-dimensional theorem.

[Work record] The first numerical run passed the mathematical assertions and wrote checks.json, but summary serialization failed because a count was a numpy integer. The count was explicitly converted to a Python integer. The final run completes successfully and writes both files. No failed mathematical check was discarded.

[Unperformed] No heavy local diagonalization, Typhon job, large interacting extrapolation, finite-temperature large-ring multi-pinch determinant or independent referee review was performed. The general reconnection lemma remains open.
