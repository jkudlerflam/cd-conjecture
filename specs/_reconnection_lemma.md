## The reconnection lemma for interacting chains

Setting: one-site translation-invariant nearest-neighbour chain H = sum_x h_{x,x+1} on (C^d)^L, M the gauge-optimized half-spread of h. Ring partition functions Z(L,beta) = Tr exp(-beta H_L). R_L(beta) = log[Z(L,beta)^2 / Z(2L,beta)]; a bound R_L(beta) <= K(d) + k M f(d) beta/L for all vacuum-compatible L and beta implies A = sum_a c_a v_a <= (4k/pi) M f(d) (C116).

C121 (`../work_C121/`) proved: inserting N pinchings P(.) = sum_s |s><s| . |s><s| of a seam site, at imaginary times spaced tau = beta/N, gives positive recorded histories q_{L,s}, s in [d]^N, with d^(-N) Z_L <= Q_L = sum_s q_{L,s} <= Z_L, and similarly on the doubled ring with two seam sites (pinched at both reconnection points, Q_D); the history Cauchy–Schwarz gives Q_L^2 <= d^N sum_s q_{L,s}^2; hence the exact reduction

    R_L <= 3 N log d + G,      G = log[ sum_s q_{L,s}^2 / Q_D ].

Classically (diagonal bond) G <= 0 because fixing the seam state cuts the ring and sum_s q_{L,s}^2 is a sum of products of the two halves; quantum mechanically a coarse record does not cut the ring between measurements, and G <= 0 is false (C121). The lemma to prove:

    G <= b_0(d) + b_1(d) N,    N = max(1, ceil(beta M/(a L))) for some fixed a > 0,

which gives A <= [4/(pi a)][3 log d + b_1(d)] M. Known: it holds for XX and independent XX stacks (C123, via C118); exact XX computations to L = 8192 give G/N bounded and saturating in N (0.69 to 0.88 at a = 1); interacting exact diagonalization (Heisenberg, XXZ, SU(3) on rings to L = 6) gives G/N of order one with no growth with L. Known obstructions (do not repeat): operator-order comparison of the two Gibbs operators loses a factor L already on XX (C122); Araki-type imaginary-time locality constants grow like exp(exp(2 beta M)); local energy functionals, gradient norms, fixed-support bootstraps and dimension-bounded divergences fail (REPORT §9.68 to §9.85).

Suggested first case: stoquastic chains (nonpositive off-diagonal matrix elements of h in a product basis, possibly after a local unitary), where the Trotterized path integral has nonnegative weights and the classical argument can be attempted on the seam column with the history structure intact; identify exactly which cross terms (C121 calls them B_xx B_yy versus B_xy B_yx) survive and bound them per interval using that the interval length tau is of order L/M and the ring's low-energy level spacing is of order v/L. Then general chains.
