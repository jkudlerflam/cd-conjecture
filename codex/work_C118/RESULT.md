# C118 result

2026-10-29.

[Open] The requested bound for arbitrary quantum chains remains unproved. I obtained explicit all-temperature bounds for XX and classical chains, and identified failures in the direct seam, convexity and information-capacity arguments. I did not find a counterexample to the general ratio estimate or to A <= C M f(d).

[Proved] For the periodic XX chain h = Sx Sx + Sy Sy, with d = 2 and M = 1/2, every even L >= 4 and every beta >= 0 satisfy

    R_L(beta) <= 2 log 2 + (5 pi/3) beta/L.

This has the requested form with K(2) = 2 log 2, f(2) = 1, and k = 10 pi/3. The proof includes both Jordan-Wigner parity sectors. It bounds the shifted-grid quadrature error of log[1+exp(-beta cos k)] using integral |f''| <= 4 beta, then cancels the bulk pressure. It extends to independent XX tensor products with K(d) = 2 log d and the same k and f = 1. It does not improve C101's settled bound on A for free chains.

[Proved] For every diagonal classical nearest-neighbor bond, including oriented nonsymmetric bonds,

    R_L(beta) <= log d.

The spatial transfer matrix is entrywise nonnegative. Cauchy-Schwarz on the diagonal of T^L gives (Tr T^L)^2 <= d Tr T^(2L), without assuming that T is Hermitian or positive semidefinite.

[Proved] The general target at all beta would force L(E_(2L)-2E_L) <= k M f(d) at each compatible L. Thus the finite-size ground-energy problem remains necessary even after keeping the full partition function.

[Proved] Averaging a reconnection over L seam positions gives a valid Jensen bound, but does not make its norm O(M/L). Explicit product states give an averaged-seam norm at least 2M for XX, Heisenberg and SU(3), independent of L. Constant partial swaps also fail to define a continuous symmetry, with nonzero commutator matrix elements exhibited in PROOF.md. Local rank alone bounds neither the norm nor the inverse overlap of a Euclidean insertion.

[Proved] R_L is neither globally convex nor globally concave in beta already for XX L = 4. The first excitation gaps establish negative curvature at low temperature, while R'_L grows from zero to a positive limit. A continuum family of n slow c = 1 sectors also shows that the CFT scaling form and a bound at a single tau do not, by themselves, bound the large-tau slope. This is a limitation of that inference, not a microscopic counterexample at fixed d.

[Proved] The exact information identity is

    R_L = beta <V>_rho + I(A:B)_rho
          + D(rho_A || gamma_L) + D(rho_B || gamma_L),

where rho is the doubled-ring Gibbs state, gamma_L is the L-ring Gibbs state, and V reconnects the two rings. The marginal relative entropies have no bound from Hilbert-space dimension alone. The alternative identity R_L = beta<V>_sigma-D(sigma||rho) gives an upper bound only after retaining the seam energy. Neither identity turns R into a capacity-limited mutual information.

[Numerical] All completed checks pass where the associated theorem applies. The proposed convexity, concavity, small averaged-seam norm, and R <= mutual information inequalities fail on the requested models. Full-spectrum ED covers XX/Heisenberg through length 8 and SU(3) through length 6. Linear-cost XX full traces reach paired length 8192. The XX trace agrees with small ED to 2.14e-14, and the information identities agree to 2.60e-12. The trial interacting extension K(d) = 2 log d, k = 10 pi/3, f(d) = log_2 d passes the finite sample, but no general theorem follows from this check.

[Unperformed] Larger interacting ED is prepared in run_typhon.slurm, including SU(3) lengths 6 and 12. It was not submitted, and no Typhon connection was attempted. No large computation ran locally. No file outside work_C118 was modified.

[Work record] Full derivations are in [PROOF.md](PROOF.md), numerical scope and failed candidates in [CHECKS.md](CHECKS.md), and the reproducible implementation in [checks.py](checks.py). The remaining task is a quantum estimate controlling the finite-size free-energy correction or the cancellation between seam energy and relative entropy at order beta/L. No independent referee review was performed.
