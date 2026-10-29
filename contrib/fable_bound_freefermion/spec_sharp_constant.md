# Problem: sharp bound on the sum of Fermi velocities of a nearest-neighbour Bloch Hamiltonian

Self-contained. No access to any other material is assumed.

## Setting

Let n >= 1. Let A_0 be an n x n complex Hermitian matrix and B an arbitrary n x n complex matrix.
For k in [-pi, pi) put z = e^{ik} and define the Hermitian matrix

    h(k) = A_0 + z B + conj(z) B^dagger .

Assume det h(k) is not identically zero in k (so the polynomial z^n det h(k) = det(B^dagger + z A_0 + z^2 B),
of degree at most 2n in z, has finitely many roots on the unit circle).

Fermi points: the pairs (k_F, u) with h(k_F) u = 0, u a unit vector in ker h(k_F). At each such k_F pick an
orthonormal basis of ker h(k_F) that diagonalizes the compressed matrix P_F h'(k_F) P_F, where
h'(k) = i (z B - conj(z) B^dagger) and P_F is the orthogonal projector on ker h(k_F). The velocities at k_F
are the eigenvalues v of P_F h'(k_F) P_F (equivalently v = <u| h'(k_F) |u> for the chosen basis vectors u;
these are the slopes at k_F of the analytic eigenvalue branches of h that vanish at k_F).

Define

    S(A_0, B) = sum over all Fermi points k_F and all kernel basis vectors u of |v| .

Known facts (you may use them, or reprove them):
(F1) The number of terms in S is at most 2 rank(B) <= 2n (multiplicity of a root of det counts the kernel dimension).
(F2) For every term, |v| <= 2 |<u|B|u>| <= 2 ||B||_op.  Indeed f(k) = <u|h(k)|u> = alpha + 2 Re(z beta), alpha = <u|A_0|u>,
     beta = <u|B|u>, vanishes at k_F and v = f'(k_F), so v^2 = 4|beta|^2 - alpha^2.
(F3) Drude/f-sum bound: S <= 2 pi ||B||_1, where ||B||_1 is the trace norm (sum of singular values).
     Proof sketch: E(k) = sum of negative eigenvalues of h(k) = tr(P(k) h(k)) with P the negative-energy projector;
     E' = tr(P h') jumps by -|v| at each Fermi term; E'' = tr(P h'') + tr(P' h') with tr(P' h') <= 0 and
     h'' = -(zB + conj(z) B^dagger); integrating over the circle gives S <= - int dk tr(P (zB + conj(z)B^dagger)) <= 2 pi ||B||_1.
(F4) Equality cases known: A_0 = 0 and B normal give S = 4 ||B||_1 (each eigenvalue b_s of B gives two Fermi points
     with |v| = 2|b_s|). A_0 = 0 and B = sigma |u><v| with u orthogonal to v gives flat bands, S = 0.
(F5) Numerical optimization over (A_0, B) for n <= 3 never produced S > 4 ||B||_1.

## Question

Prove or disprove:  S(A_0, B) <= 4 ||B||_1  for all n, all Hermitian A_0, all B.

If true, give a complete proof with all steps. If false, give an explicit counterexample (A_0, B) with n as small as
possible, listing the Fermi points, kernel vectors and velocities, and the value of S / ||B||_1.

Possibly useful reformulations (not required):
(R1) The relation for distinct Fermi points F, G: <u_G| B^dagger |u_F> = z_F z_G <u_G| B |u_F>, obtained by
     subtracting <u_G|h(k_F)|u_F> = 0 and <u_G|h(k_G)|u_F> = 0. Hence the matrix Gamma = Z U^dagger B U
     (U = [u_F], Z = diag(z_F)) is Hermitian off the diagonal, and its diagonal is -alpha_F/2 - i v_F/2.
(R2) When A_0 = 0 the Fermi vectors come in pairs (z, -z) with the same u, distinct u's are B-orthogonal
     (<u_s|B|u_t> = 0 for s != t) and satisfy B^dagger u_s = lambda_s B u_s with |lambda_s| = 1; then S = 4 sum_s |<u_s|B|u_s>|,
     and the claim becomes sum_s |<u_s|B|u_s>| <= ||B||_1 for such a B-orthogonal family of unit vectors.
(R3) Convex geometry (A_0 = 0): let W = { tr(P B) : 0 <= P <= 1 } in the complex plane (a convex set). Then
     E(k) = min_P tr(P h(k)) = 2 min_{w in W} Re(z w) is minus twice a planar support function, Fermi points are
     directions in which W has a flat edge, and S = 2 x (total length of all edges of W). The claim is then:
     the total edge length of W is at most 2 ||B||_1. For normal B, W is the zonotope sum_s [0, b_s] and this is equality.
(R4) In general, with A_0, the same picture holds in R x C for K = { (tr(P A_0), tr(P B)) } with directions (1, 2z).

Also of interest: if the inequality is true, characterize all equality cases.
