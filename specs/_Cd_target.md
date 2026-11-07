## Shared target for the C(d) round (C126 onward, 2026-11-06)

Conventions are those of `_bound_background.md` (read it with all updates) and `_reconnection_lemma.md`. Read REPORT.md sections 9.76 to 9.88 before starting.

Target. Define

    C(d) = sup over Hermitian bonds h on C^d (x) C^d with M(h) > 0, and over vacuum-compatible ring lengths L, of L (E_{2L} - 2 E_L) / M.

Prove that C(d) is finite for every d, with any explicit function of d, however large (towers of exponentials are acceptable). On the CFT subsequence L (E_{2L} - 2 E_L) -> pi A/4, so this gives A <= (4/pi) C(d) M. Two variants also count. (a) A Casimir lower bound E_L >= e_0 L - K(d) M/L on vacuum-compatible L, which gives A <= (6/pi) K(d) M directly. (b) The finite-temperature ratio bound R_L(beta) <= K(d) + kappa(d) M beta/L of C116. The zero-temperature target is variational, it asks for a state on the 2L ring, built from the L ring, with energy at most 2 E_L + C(d) M/L. The proof must give a precise meaning to "vacuum-compatible" (for example lengths on which the ground state carries no forced defect) and say which one it uses.

Known. Product gluing gives E_{2L} - 2 E_L <= 2M, so any finite range L <= L_0 is trivial with C = 2 L_0. C = 0 for classical and for commuting chains (C118, C124). C <= 10 pi/3 for XX, independent XX stacks and the C106 two-band family (C119). Heisenberg gives L (E_{2L} - 2 E_L)/M = 4, 2.79, ..., pi^2/4 at L = 2, 4, ..., infinity, so C(2) >= 4, and some growth in d is necessary (SU(N) symmetric-representation chains). For chains with a positive spatial transfer operator, D(L) = E_{2L} - 2 E_L is nonnegative and nonincreasing in L (C125).

Lessons that any proof must respect (do not repeat these routes).
1. Abstract crossed-channel data do not suffice (C125). Positivity of the transfer operator, unitary Trotter-circle rotations, a crossed cone condition, an entire tilt, the trivial bounds and the infrared CFT form are all satisfied by a positive spectral family with A unbounded at fixed d and M. Its multiplicity of low crossed states per crossed momentum is unbounded. The proof has to use the fixed-d locality of the microscopic chain beyond these properties.
2. Local energy functionals (modulations, gradients, Rindler weights, fixed-support bootstraps) see near-critical massive modes at the order of the anomaly (sections 9.78, 9.79).
3. Unfiltered counting of low-lying ring states fails (C115), and the number of infrared sectors is not bounded in terms of d. Fast sectors can have arbitrarily long crossover lengths (C106, C125), so no proof may assume bounded crossovers.
4. The coarse-seam reconnection lemma is equivalent to the target (section 9.86). Operator-order comparison of Gibbs operators loses a factor L (C122). Araki-type constants are doubly exponential in beta M. Extensive comparisons (Berezin–Lieb type) lose order L beta.
5. Convexity of Theta(q) = -L E_L in q = L^2 beyond a fixed L_0(d) would suffice but is false in general (XXZ near Delta = -1, C125).

Heuristic orientation (not proved). If every sector entered its CFT regime by temperature of order v_a (crossover length of order one), the entropy bound s(T) <= log d applied at T = V 2^{-k} would allow at most (6/pi) log d modes per dyadic velocity shell [V 2^{-k-1}, V 2^{-k}], with V = 4eM, and A <= (12/pi) V log d would follow. Every obstruction above is a long-crossover effect. A successful proof therefore has to control the crossover regime, for instance by bounding the doubling defect at scale L using only data at scale L.

Route considered and set aside. The identity sum_x sin^2(pi (j - x)/L) = L/2 gives H_L = (2/L) sum_x H_SSD^(x), a sum of translated sine-square-deformed rings, so E_L >= 2 lambda_min(H_SSD) for every chain, with equality within CFT and, by Katsura's theorem, exactly for free fermions. But lambda_min(H_SSD) is a local energy functional of the envelope, so lesson 2 is expected to spoil the needed bound lambda_min(H_SSD) >= e_0 L/2 - K M/L (heuristic, not checked numerically).

Routes in this round, one per agent, are listed in the spec of each run.
