# C121 result

2026-10-29.

[Proved] Coarse seam pinching has the required small error at arbitrary spacing. For the positive imaginary-time measurement functional Q defined in PROOF Section 3, N pinchings on a d-dimensional seam obey

    d^(-N) Z_L <= Q_L <= Z_L.

On the doubled ring, pinching both seam sites gives d^(-2N) Z_(2L) <= Q_D <= Z_(2L). The proof uses the operator pinching inequality for the lower bound and Schatten Holder after a unitary-average expansion for the upper bound. It needs neither a finite-size gap nor an approximate Markov property. Taking N = max(1, ceil(beta M/(a L))) gives a log-trace error at most (1+beta M/(a L)) log d per seam, for any fixed a > 0.

[Proved] Both formulations of the history Cauchy-Schwarz step work. Coherent Euclidean histories a_s sum to Z_L and satisfy Z_L^2 <= d^N sum |a_s|^2. Positive recorded measurement histories q_s sum to Q_L and satisfy Q_L^2 <= d^N sum q_s^2. With

    G = log[(sum_s q_s^2)/Q_D],

one obtains the exact reduction

    R_L <= 3N log d + G.

[Disproved] The proposed no-cost reconnection step G <= 0 does not follow from positivity and fails on XX, Heisenberg and SU(3), already for one insertion at tau = L/M. The respective measured values are 0.762551, 1.888216 and 2.292663 on L = 4, 4 and 3. The coherent-history gluing step fails there too. An analytic finite XX witness is given in PROOF Section 6. At L = 8192 and beta = 4L, even log[(sum q_s^2)/Z_(2L)] = 0.986692 > 0. Thus this failure persists on large rings and does not depend on measuring the doubled-ring denominator.

[Proved mechanism] A coarse record does not specify the intervening seam evolution. Conditioning the seam at two times leaves virtual couplings between its two sides. In a positive Trotter network, grouping microscopic histories introduces cross terms B_xx B_yy that are not covered by the classical comparison with B_xy B_yx. Stoquasticity removes signs but does not control these terms. The even XX and Heisenberg tests are unitarily equivalent to stoquastic tests with the same seam projectors.

[Proved] Continuously pinching the seam instead converges to Tr exp[-beta D(H)]. On XX this removes the two incident hopping bonds, giving

    E_Zeno-E_L = 1/2+(1/2)tan(pi/(2L)) -> M.

Its low-temperature log-trace cost is order beta M. Sparse measurements are cheap because they leave between-slice propagation, while continuous measurements remove that propagation and change the ground energy by order M.

[Open precise lemma] For some explicit a > 0, prove

    sum_s q_(L,s)^2 <= exp[b_0(d)+b_1(d)N] Q_D,
    N = max(1, ceil(beta M/(a L))),

with finite constants independent of the Hamiltonian, ring size and crossover scales. If this holds, the target follows with

    K(d) = b_0(d)+3 log d+b_1(d),
    k f(d) = [3 log d+b_1(d)]/a,
    A <= [4/(pi a)] [3 log d+b_1(d)] M.

In fact R_L-3N log d <= G <= R_L+2N log d, so this lemma is comparable in strength to the target, not a separate solution of its infrared content. The zero-cost version is false at a = 1 in the tested examples. These results do not exclude a sufficiently small universal a or a controlled finite cost per slice. A coherent-history alternative with a smaller formal coefficient is stated in PROOF Section 8.

[Numerical scope] Full ED covers XX and Heisenberg through length eight and SU(3) through length six. Forty-eight coarse-channel cases, direct history enumeration, eighteen general algebraic controls and continuous-pinching checks all satisfy the proved inequalities. An exact one-pinch XX formula is calibrated against ED to 3.22e-15 and evaluated through L = 8192, with full doubled-ring traces through length 16384. No large interacting or multi-pinch free-fermion extrapolation is claimed.

[Open] No general bound on R_L or A beyond the inherited results is proved, and no counterexample to the target bound is found. The task is reduced to the spatial reconnection lemma rather than the cost of temporal pinching. Full derivations are in [PROOF.md](PROOF.md), tests and failed candidates in [CHECKS.md](CHECKS.md), and reproduction in [checks.py](checks.py).
