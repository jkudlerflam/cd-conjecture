# Independent audit of the branch-source construction

Status: **[checked algebra and conditional source matching]**. This is an
independent AI review by the exact-alternative subagent, 2026-10-28. No human
review, infrared theorem, or central-charge determination is claimed.

I read [BRANCH_SOURCES.md](../sources/BRANCH_SOURCES.md) and the full
[finite-cell matching theorem](../../direct_attack_20261028/infrared/FINITE_CELL_MATCHING.md).
The source result is correct conditional on that theorem and its strengthened
coefficient norms. I found no issue affecting its stated scope.

## Algebraic checks

1. Rewriting (1) using sin^2((k-k_b)/2) cancels all factors of two. Each
   numerator factor has Fourier coefficients 1/2 and
   -exp(plus-or-minus i k_b)/4, so its coefficient l1 norm is exactly one.
   Convolution proves (2), including the separation bound. Degree is N-1.
2. Hermiticity of (3) follows after taking the adjoint and replacing r by -r,
   using p_(-r)=conjugate(p_r). With c_x=sum_k exp(ikx)c_k/sqrt(L) and source
   Fourier factor exp(-iqx), its row k, column k+q vertex is exactly
   [p(k)+p(k+q)]/2. The common scalar orbital matrix removes dependence on
   the Fermi eigenvector in the diagonal zero-transfer vertex.
3. Every individual fermion bilinear is even and has an ordinary spin
   representation supported on its endpoint interval, including the needed
   parity factors. A parity-preserving code isometry preserves this property.
   The operator G_I O G_I^* is reducing for the all-code projection, and has
   norm ||O||. A single bilinear of range R uses R+1 cells; the complete
   symmetrized density uses at most 2R+1 cells. The 46-cell and 91-cell counts
   for N=46 are correct.
4. The derivative estimate |p'_a| <= R/D_a follows directly by multiplying
   each Fourier coefficient by |r|<=R. The stated matrix row-sum estimate
   and Neumann condition therefore hold. They require the polynomials to be
   held at their reference crossing momenta while evaluating the shifted
   crossings, as the note does.

A presentation detail worth retaining explicitly: to preserve the shorter
46-cell support of *each* bilinear, lift each bilinear on its own endpoint
interval. Their sum has the norm bound obtained by the triangle inequality.
Alternatively lift the entire 91-cell density at once, which preserves its
exact norm but introduces code projectors over the entire 91-cell interval.
Both are valid local reducing probes; they should not be conflated if exact
termwise support or exact norm equality of the sum is used later.

## Analytic source check

The source remainder D_j in the auxiliary embedding vanishes on the
auxiliary vacuum on both sides. This remains true for complex linear
combinations D(z). After normal ordering, every nonzero monomial has an
auxiliary creator and an auxiliary annihilator. The zero-temperature
rightmost-vertex argument thus applies even with external logical fields
and complex source parameters. Self-adjointness at complex z is not needed
for the convergent analytic expansion.

On the specified polydisk, k+sum rho_j s_j+r <= 1/(4 omega). The finite-cell
theorem's Lipschitz argument can be applied with K+D(z) as the reference:
its further disk in the leakage direction stays inside 1/(2 omega), giving
||F(z)||<=4r. The Banach-valued multivariable Cauchy formula then gives

    ||partial^alpha F(0)|| <= 4r alpha! product_j rho_j^(-alpha_j).

This is a bound on the *matching correction* to source coefficients. The
pure logical source factors out and is not claimed to be O(r). Bounded
time-dependent test functions are admissible as instantaneous vertices
when their supremum norms are included in the stated input coefficient
norms. Scalar derivatives need the separate marker norm and radius reserve,
which the note explicitly records.

The source statement does not repair an incorrect Hamiltonian matching
hypothesis, establish pointwise correlators, or supply emergent current
Ward identities. Its explicit warning about additional leakage for an
untransformed source is necessary and correct.

## Finite independent check

Run `python3 contrib/decisive_attack_20261028/exact/source_checks.py` with
numpy available. The script constructs three filters on an 11-site ring,
checks their interpolation and l1 bounds, builds the one-particle matrices
of the local densities, and verifies the complete Fourier vertex formula
for every allowed transferred momentum. Output is saved in
`source_results.json`. This supplements the algebra; it is not an RG test.
