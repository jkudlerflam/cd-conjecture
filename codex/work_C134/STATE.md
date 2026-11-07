# C134 handoff

2026-11-07. Read RESULT.md, then PROOF.md. All files are confined to work_C134, as required. No project-level memory edit, cluster submission, or independent referee pass was made.

[Proved] New compact theorem for every integer L>=2 and Delta in [-1+delta,0], 0<delta<=1:

    F_L >= L e_0 - 1200 J/(delta^3 L),
    |F_L-L e_0| <= 1200 J/(delta^3 L) on even L,
    L(F_(2L)-2F_L) <= 3000 J/delta^3.

PROOF Section 1 gives smaller explicit pointwise and compact constants. Proof uses the negative scattering kernel, a global increasing finite counting function, Neumann inversion of cumulative root discrepancy, exact full-line dressed quadrature, and clipping where rho=2D/L. The clipped derivative has bounded variation uniformly in L. Odd roots are justified by convex Yang-Yang functional plus local coordinate Bethe verification and analytic continuation of the characteristic polynomial, then Perron-Frobenius ground simplicity. This is new work in the run, not independently reviewed.

[Proved using literature input] Delta=+1/2 in the user's AF convention is the supersymmetric Delta=-1/2 of the negative-transverse convention. Odd minimizing-flux energy is -3JL/8; even energy at a flux 2pi/3 away from the minimum is -3JL/8. C129 twist cost then gives C=2pi^2/9 in units of J for all L. References and normalization are in PROOF Section 11.

[Proved, inherited] Delta=0 has C=4-2sqrt(2) in units of J. For odd half-filled doubling its defect is zero, because F_(2L) uses n=L, not 2floor(L/2). Delta=-1 has zero defect. C129's (Delta+1)L^2 band is size-dependent and does not cover a fixed neighborhood.

[Assumed, unresolved] No uniform constant over (-1,1], no compact interval crossing zero, and no separate Delta=1 proof. Uniform estimates as Delta approaches 1 would extend to the endpoint by finite-dimensional continuity. C125 monotonicity of D alone is insufficient; its positive transfer construction is only for Delta>=0 on even rings. The new compact constant diverges at -1 and cannot yet be patched to the crossover band.

[Proved, literature audit] Kozlowski half-filled result is rate-free convergence, whereas Corollary 4.3 has the energy term below half filling. DKKMT Theorems 2 and 22 exclude exact half filling. DDV Section 7.1 performs the shifted-tail scaling limit without the explicit norm estimates needed for this target. Later rigorous QTM result 1811.12020 is high-temperature only. Several older full texts could not be accessed; LITERATURE.md distinguishes full-text checks from abstracts and indexed passages. No exhaustive no-theorem claim is made.

[Numerical] 152 paired Bethe cases on all eight requested anisotropies through L=129, paired length 258. ED for every L=2,...,12, both odd neighboring sectors, max energy discrepancy 8.44e-15. Max Q=4 at L=2, Delta=1. Root, quadrature, and special-twist checks pass; see CHECKS.md. All checks are floating point.

[Assumed, pending computation] run_typhon.slurm prepares eight jobs through L=4097, paired 8194, not submitted. Large solver convergence and cancellation still need inspection. No large local computation was run.
