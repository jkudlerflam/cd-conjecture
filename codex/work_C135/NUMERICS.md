# C135 numerical methods and coverage

[Numerical] The scripts use installed NumPy 2.4.6, SciPy 1.17.1 and CVXPY 1.9.2 with Clarabel and SCS. Nothing was installed. All computations and generated files are inside work_C135. Large continuations are prepared but were not submitted, since this sandbox has no cluster access. These are floating-point bounds with explicit algebraic repairs, not rigorous interval enclosures.

## Models and finite coverage

[Numerical] models.json contains 62 bonds, with full real and imaginary matrix entries. The list contains Heisenberg, XXZ at 0.5, 2, -0.9, -0.99 and -0.9999, XX, XX with field b=1/2, XXZ at Delta=0.5 with DM_z=0.7, the chiral bond SzSz+DM_z-I/4, critical TFIM, and the 21 TFIM bonds g=1+c/s for c=-2,0.5,5 and s=2,...,8. The latter are separate fixed Hamiltonians. Reading their entries at hierarchy level n=s gives the requested near-critical diagonal sequence without silently changing h in a two-level LTI difference.

[Numerical] The seed for the 30 random Hermitian bonds is 1352026. Bonds 00 through 14 are real, of which 00 through 07 are also swap symmetric. Bonds 15 through 29 are complex. A scalar is removed, then each bond is radially divided by a computed feasible gauge half-spread. Thus M=1 to numerical accuracy. This is not a uniform distribution on the quotient sphere. The field and DM components are unrestricted for the nonsymmetric sample.

[Numerical] All 62 bonds have ring energies through L=10. Heisenberg, XX, XX with field and critical TFIM were diagonalized through L=12. The 24 free models, XX, XX with field and 22 TFIM bonds, have energy spectra extended through L=20 using their finite momentum formulas, checked against every available diagonalization. This adds energies, not new Petz dictionary points at large L. Interacting and generic rings beyond the stated sizes have not been computed.

[Numerical] Symmetry-reduced interacting bonds have LTI solves through n=8. Real random bonds have solves through n=6. The complex random set has solves through n=4, with one extra n=5 solve for random_15. The first generic complex n=6 run was stopped after it exceeded the useful local benchmark. The real formulation of the exact quadratic covariance SDP supplies free-model LTI energies through n=10. Full spin SDP records are preserved in full_spin_lti for those models, and their dual matrices remain in matrices/.

[Numerical] All 62 bonds have hull upper distances at L=3,...,8. Lower distances were computed with n=L+1 through L=5 for real target marginals and through L=3 for complex targets. The earlier standard-model run also reached L=7 for XX, Heisenberg, XX with field and critical TFIM. The exact per-model coverage and every missing entry are in TABLES.md and summary.json. No missing distance is assigned zero.

## Bounds, residuals and repairs

[Proved] For any Hermitian numerical X, the number

    lambda_min(h tensor I+X tensor I-I tensor X)

is an admissible LTI dual objective in exact arithmetic, regardless of solver status. core.py computes it afresh and subtracts a small floating-point safety allowance. That allowance is not an interval-arithmetic error proof. The files retain the primal objective, corrected lower value, solver status, LTI equality residual, minimum primal eigenvalue, trace residual, complementary-slackness residual and the SCS residuals where available.

[Proved] The orthogonal projection of a density matrix onto the LTI affine space can be performed in its Pauli expansion. The endpoint-trace equation equates coefficients of strings that differ by translating the same nonidentity core through identity padding. Average each such orbit and set the identity coefficient to 2^(-n). This gives a Hermitian trace-one LTI matrix q. If its minimum eigenvalue is at least -delta, then

    q_feasible=(q+delta I)/(1+2^n delta)

is positive, trace one and LTI. The script uses an additional small positive allowance in delta. Its energy supplies the upper bound on E_LTI(n) used when subtracting levels. Specifically the reported lower bound on n²(e0-E_LTI(n))/M uses

    n²[dual_lower(n+2)-primal_upper(n)]/M.

Subtracting two raw solver objectives would not have this guarantee.

[Proved] The gauge-scale primal gives a feasible half-spread upper bound on M. For a lower bound, project the difference of its two cone dual matrices onto the annihilator of the scalar and three divergence directions. Divide the result by the larger of one and its nuclear norm. The resulting D satisfies Tr D=Tr(DG_a)=0 and ||D||_1<=1, so Tr(hD)<=M. Records include both endpoints. report.py uses the appropriate endpoint for the sign and direction of each normalized bound. These remain floating-point checks.

[Proved] For quadratic bonds the full spin LTI optimization equals optimization over physical block-Toeplitz generalized covariances. One direction follows by taking the covariance of any parity-averaged LTI state. Conversely, every physical covariance defines a Gaussian state. Identical overlapping covariance blocks give identical Gaussian reduced states; parity removes odd spin moments, so the corresponding spin state is LTI. The energy is quadratic and therefore unchanged. This proves that quadratic_lti.py is an exact reformulation, not an additional relaxation.

[Proved] In the covariance SDP, 0<=G<=I, overlapping blocks agree, and G+JG^T J=I. The objective is Tr(AG), where A is obtained by resolving the actual two-spin bond into quadratic CAR monomials. The coefficient reconstruction is checked in matrix norm. Let X,Y be arbitrary equality multipliers, and let P,Q>=0 be repaired cone multipliers. Define

    R=A+D(X)+Y+JY^T J-P+Q.

For every feasible G,

    Tr(AG)>=-Tr Y-Tr Q+Tr(RG)
            >=-Tr Y-Tr Q+sum_j min(lambda_j(R),0).

This is the corrected dual lower bound used in quadratic_lti.py. The covariance primal is repaired by averaging block diagonals, imposing the particle-hole affine relation, and mixing with I/2. Its energy supplies an upper bound. The early complex formulation had poor numerical dual conditioning; its outputs were replaced with the real formulation, which is valid for all free bonds in this sample. No failed complex-formulation values enter the final tables.

[Numerical] A first full-spin SCS solve at TFIM g=3.5, n=8 returned optimal_inaccurate and a repaired primal-dual gap about 0.0019. It is retained in full_spin_lti but superseded for the energy tables by the exact covariance reformulation. A real random_13 n=6 solve was refined with Clarabel. Some generic solves still return optimal_inaccurate. Their corrected bounds are retained as exploratory. No convergence claim is inferred from a small equality residual alone.

[Numerical, audit trail] quadratic_lti.txt contains the early, poorly conditioned complex-covariance formulation and is not used for final claims. quadratic_lti_real.txt is an intermediate real-formulation log. quadratic_lti_final.txt and the version-2 covariance records contain the final energies. optimize_lti.txt records a gauge-solver failure; optimize_lti_finish.txt contains the resumed search after adding the SCS fallback. gauge_refinement.txt records independent direct-dual improvements of loose gauge bounds. The failures are retained rather than silently deleted.

[Proved] For the trace-distance lower bound, the primal minimizes Tr S with S>=(rho_target-tau) and S>=-(rho_target-tau), tau the first two-site marginal of an LTI_n state. From its dual take a Hermitian Y and an equality multiplier X. Normalize both by max(1,||Y||). Then

    Tr(Y rho_target)
      +lambda_min(-Y tensor I+X tensor I-I tensor X)

is a lower bound on the distance to LTI_n. The script also computes the distance to a repaired primal marginal as an upper bound on that finite relaxation, for diagnostics only. It is not an upper bound on distance to Omega_inf.

[Proved] The hull program minimizes the trace norm over convex weights on the dictionary. Clipping negative numerical weights and renormalizing their sum gives a convex combination and a direct feasible upper objective. The dictionary consists of Petz-periodized ring marginals for every computed ring and central window, shared across all test bonds, plus the six pure product states polarized along the coordinate axes. The latter are uniform MPS of bond dimension one. Each exact Petz point belongs to Omega_inf by PROOF.md section 2. Numerical support cutoffs and trace residuals are recorded; points with trace residual above 1e-8 are excluded.

[Numerical limitation] Sparse diagonalization supplies a Rayleigh value and an eigenvector residual. A small residual does not independently certify that no lower eigenvalue was missed. These entries are numerical ground-energy estimates, not rigorous spectral enclosures. Translation averaging of any selected normalized ground vector gives a valid TI state, even if a degenerate eigensolver misses part of the ground space. The metadata distinguish sampled multiplicity and whether all returned eigenvalues were degenerate. No unique exposed point is asserted in a degenerate case.

## Free ring formulas and checks

[Proved] For XX with field b, the spin ring energy is the minimum over occupations m of

    sum_(m lowest k) cos(k)+b(m-L/2),
    k=2pi(j+phi)/L, phi=0 for m odd and phi=1/2 for m even.

Both Jordan-Wigner parity sectors are included. Its infinite density for |b|<=1 is -[sqrt(1-b²)+b arcsin b]/pi. For the TFIM convention h=-XX-g(ZI+IZ)/2 and g>=0, the lowest even sector gives

    E_L=-sum_(j=0)^(L-1) sqrt(1+g²-2g cos(2pi(j+1/2)/L)),
    e0=-(1/pi) integral_0^pi sqrt(1+g²-2g cos k) dk.

At g=0 the degenerate sectors have the same ground energy. The sector selection follows by diagonalizing the quadratic periodic and antiperiodic Majorana matrices with the occupation-parity constraint. The antiperiodic vacuum is admissible and minimizes for g>=0. free_rings.py checks the finite formula against all full spin diagonalizations available in this run, including odd L. The integral is evaluated with the complete elliptic integral, without fitting a ring extrapolation.

[Proved, sector check] For this TFIM finite-spectrum formula choose down-spin fermion occupation (1-Z)/2, whose parity is product Z. This differs by particle-hole convention from the occupation covariance used in the counterexample and code's coefficient reconstruction. For g>0 rotate the spin axes to write the TFIM as -sum ZZ-g sum X. Its matrix in the Z basis is irreducible with nonpositive off-diagonal entries, so its ground vector is unique and strictly positive. Global spin flip commutes with the Hamiltonian and preserves positivity, forcing that vector into the even spin-flip sector. Rotating back identifies the even Jordan-Wigner sector in this convention. At g=0 the result follows by continuity. This supplies the ground-sector selection used in the formula, independently of the numerical comparisons.

[Numerical] math_checks.json records independent full-spin checks of the exact three-site dual, the seam-window counterexample and the rational four-site Gaussian state. It also checks covariance-symbol positivity, its particle-hole identity and the first two Fourier coefficients for XX, paired TFIM and an interacting Heisenberg input covariance. The latter checks positivity of the covariance construction, not preservation of the interacting energy. analytic_checks.json verifies the explicit XX certificate through n=32 in one-particle form and through n=8 in full spin form.

## Searches and unrun continuation

[Numerical, exploratory] optimize.py uses the eight-angle C131 coordinates (Jx,Jy,Jz,Dx,Dy,Dz,bx,by,bz), followed by radial gauge normalization. Four local starts were run for the n=2 to n=4 corrected LTI difference and for the Petz upper endpoint of the directional bracket at L=6,8,10. Budgets were 100 evaluations per LTI start and 300 per ring-upper start. An upper-endpoint maximization locates difficult periodization trials, not a lower witness against (K). These limited searches are not global maxima. Their convergence flags and full maximizing coordinates are in search/.

[Numerical limitation] The requested multistart optimization at the largest symmetry-reduced n cannot be substituted for optimization over all eight angles, since a generic perturbation destroys that symmetry. The local all-angle search therefore used the smaller n=2 to 4 pair. The n=6 to 8 all-angle search is implemented in search.slurm and was not run. No conclusion about drift with n follows from the local n=2 search alone.

[Assumed continuation, not executed] energy_scan.slurm runs each of the 62 bonds separately through ring L=20 and LTI n=8, resuming existing files. deeper_lti.slurm requests n=9,10 for the two-level differences at n=7,8, but generic complex n=10 has not been benchmarked and may exceed that allocation. distance_scan.slurm uses the completed ring dictionary and requests n=9 lower distances for L=8. search.slurm requests eight starts and 300 evaluations for the n=6 to n=8 all-angle LTI difference, plus the three ring-upper searches. These jobs use no package installation and must be run from a copy of this folder in the cluster environment. No cluster job was submitted.

[Proved limitation] With n<=8, positive directional lower witnesses at L=8 or 10 are unavailable. Increasing ring sizes alone cannot fix this, since F1 forces E_LTI(n)<=E_L/L whenever n<=L. The upper endpoint can still be optimized there. A future lower-endpoint search at L=10 needs n>=11 or a different certified infinite-density lower bound.
