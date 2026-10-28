# A local Clifford cursor with exact quadratic exchange

## Scope

[proved] This is a fixed-alphabet, nearest-neighbor isolated-cell
Hamiltonian with an exact product code, exact parity, and an unbounded-rank
Gaussian second-order exchange. It supplies a nonscalar source interface.
It does not produce C19's independently displaced flavour bands. Section 7
proves a restriction on the resulting hopping matrix. Cell cuts are
prescribed, not selected by a translation-invariant ring Hamiltonian.

The construction and the proofs below are new in C21. C19 supplies the
Majorana convention and the entrywise second-order contraction. C18 already
showed that a syndrome-conditioned work energy can evade the restricted
bare-Bell-check obstruction. The present cursor makes that mechanism
quadratic and of growing rank, without a history-state computation.

## 1. Complete isolated-cell rule table

[proved, definition] A site has a cursor register with eight states

    0, 0X, 0Y, 0Z, 1X, 1Y, 1Z, f

and a data qubit. The first 0 is the vacuum; the prefixes 0 and 1 in
0a and 1a label lanes. The filler f completes the port Pauli algebras.
The site dimension is 8*2=16. This is a module count, not d0 for the
requested completed chain.

Let Q project onto all seven nonvacuum cursor states. In the three-state
syndrome space X,Y,Z, define

    A = |X><Z| + |Z><X|,
    D = diag(-1,-1,+1).

Every site has the identical on-site term o. Every internal bond has the
identical transport term b. The following list specifies all matrix
elements; unlisted matrix elements vanish.

| Support | Term | Coefficient |
| :-- | :-- | :-- |
| One site | Q tensor I_data | Delta |
| One site, lane 0 | A tensor X_data | g |
| One site, lane 1 | A tensor Y_data | g |
| One site | sum_a (|1a><0a|+|0a><1a|) tensor I_data | -t |
| Two sites | T+T^dagger, defined below | -t |

The transport partial isometry is

    T = sum_(a=X,Y,Z) |0,0a><1a,0| tensor Z_left^epsilon(a),
    epsilon(X)=epsilon(Y)=1, epsilon(Z)=0.                 (1)

The identity on the right data qubit is understood. The tensor notation
in (1) groups the two cursor registers before the two data registers;
the physical operator is still on two sites. It moves a lane-1 cursor
to lane 0 of the next site, only if that site is vacant. Two occupied
sites are not moved, and neither is f.

For n sites,

    H_cell = sum_(i=1)^n o_i + sum_(i=1)^(n-1) b_(i,i+1).  (2)

Equivalently, the repeated internal bond is
h_int=(o tensor I+I tensor o)/2+b, with endpoint counterterms o/2.
For n=1, use o alone. These endpoint terms are part of this isolated-cell
specification. No site-dependent coefficient or instruction occurs.

The tested constants are Delta=3, t=1, g=3/5. The all-n proofs require
Delta>|g|+2|t|, with g and t nonzero for full hopping rank.

## 2. Initialization, all sectors, and parity

[proved] The complete code is

    G psi = |0>_cursor^tensor n tensor psi_data.             (3)

No data bit is initialized or copied. Every nonvacuum sector is excluded
by energy, including multiple cursors and any number of fillers.

To prove this, let N_c=sum_i Q_i. The syndrome potential has norm at most
|g| on each occupied site. Internal lane hopping is bounded below by
-|t| Q_i. For an intersite transport T, its initial and final spaces are
orthogonal, so

    -(T+T^dagger) >= -(T^dagger T+TT^dagger).

The projectors on the right are bounded by the lane-1 occupation on the
left plus lane-0 occupation on the right. Summing over bonds counts each
lane at most once. Therefore

    H_cell >= (Delta-|g|-2|t|) N_c.                        (4)

Every term annihilates (3). Equation (4) proves that its dimension is
exactly 2^n, its energy is exactly zero, and its complete gap is at least
Delta-|g|-2|t|. The unperturbed degeneracy is not a numerical tuning.

There are no hidden clock, program, or particle-number sectors in (2).
The ground code is the entire vacuum sector. This exclusion theorem
does not address unwanted arrangements of cell markers on a ring,
because no marker selector has been attached.

[proved] Give cursor states the parity signs

    q_c = diag(+1,-1,-1,+1,-1,-1,+1,+1)

in the alphabet order above. Physical site parity is q=q_c tensor Z.
Both syndrome mixing and data X or Y are odd, so their product is even.
Lane hopping preserves syndrome parity. The Z in (1) is even. Thus

    [H_cell, product_i q_i]=0,
    (product_i q_i)G = G P,  P=product_i Z_i.               (5)

## 3. The two port algebras, with every completion specified

[proved, definition] Start with the Bell basis on two auxiliary qubits,

    |Omega>, X_p|Omega>, Y_p|Omega>, Z_p|Omega>,
    |Omega>=(|00>+|11>)/sqrt(2).

Let B be the matrix of these four columns, and define the numerical
4-by-4 matrices F_a=B^dagger(sigma_a tensor I)B. This formula fixes all
entries and phases, including the imaginary Y-syndrome convention.

The left port acts by the same F_a on each of the ordered quadruples

    (0,0X,0Y,0Z),       (f,1X,1Y,1Z),

and as identity on data. For the unrotated right port, use instead

    (0,1X,1Y,1Z),       (f,0X,0Y,0Z).

These are direct sums of Pauli representations. Thus each port is a
genuine qubit algebra on the 16-dimensional physical site, not merely a
specified action on code vectors.

Define the on-site diagonal unitary V_R to act as i Z_data on cursor
states 1X and 1Y, and as identity on every other cursor state. Conjugate
the entire unrotated right Pauli algebra by V_R. The resulting operators
are F_a^R. The left operators are F_a^L. V_R commutes with q.

Use F^L on the first physical site and F^R on the last. At n=1 the two
algebras live on the same site and need not commute with one another.
The directed bond between cells is still the Hermitian operator

    v = sum_(a=X,Y,Z) F_a^R tensor F_a^L.                   (6)

Each summand has the correct product parity. Each port is a Pauli
representation, so ||v||=3. Directly from (3),

    G^dagger F_a^s G=0,
    G^dagger F_a^s F_b^s G=delta_ab I.                     (7)

There are no unspecified states in this definition. The role of f is
only to make the port action complete outside the source subspace.

## 4. Exact source resolvent

[proved] A port applied to (3) creates exactly one cursor and no filler.
The one-cursor sector is invariant under H_cell. Its ordered positions
are

    (site 1,lane 0), (site 1,lane 1), ...,
    (site n,lane 0), (site n,lane 1),

a path of length m=2n. Its dimension, including data and syndromes, is
3m 2^n.

Gauge away the intersite controlled Z's, with

    W_i = diag(P_(<i), P_(<i), I),  P_(<i)=Z_1...Z_(i-1),

on either lane of site i. Internal lane hopping already has identity
decoration. In this gauge all hopping is -t I and the diagonal blocks are

    Delta I + g A tensor gamma_j,
    gamma_(2i-1)=P_(<i) X_i, gamma_(2i)=P_(<i) Y_i.        (8)

The Y-syndrome channel is a scalar path. On the X,Z block, A^2=I, so

    Gamma_j=A tensor gamma_j,
    {Gamma_i,Gamma_j}=2 delta_ij I.                        (9)

Here and below the identities involving Gamma refer to that X,Z block.

[proved] The root resolvent of a scalar path with the potentials (9)
is scalar plus a Clifford vector. This follows without an expansion in
g. Let R_i(z) be the endpoint resolvent of the tail beginning at i.
For z>0, Schur complementation gives

    R_i(z)=[(Delta+z)I+g Gamma_i-t^2 R_(i+1)(z)]^(-1).

The final tail has this form directly. If R_(i+1) is a scalar plus a
linear combination of Gamma_j for j>i, the expression to invert is
a I+sum_j b_j Gamma_j. The Clifford relations give

    (aI+b.Gamma)^(-1)=(aI-b.Gamma)/(a^2-sum_j b_j^2).       (10)

The denominator is nonzero because the tail Hamiltonian is positive.
Induction proves the claim for every tail. This also shows that no
Clifford bilinear or higher product is generated in the return resolvent.
Since the space is finite dimensional, the same span contains every
spectral weight and heat kernel, by taking residues or using their
spectral decompositions.

Consequently the left source heat kernel has the exact form

    M_L(u)=a(u) diag(1,0,1) tensor I
           +y(u)|Y><Y| tensor I
           +A tensor sum_(j=1)^m f_j(u) gamma_j.           (11)

All scalar functions in this equation are real. They decay
exponentially by (4).

[proved] The right-port rotation in Section 3 gives the required
fermionic sign, not an ungraded spin exchange. At the right source its
one-cursor action is

    K=diag(iZ_n,iZ_n,I).

On the X,Z block, K^dagger W_n=diag(-iP,I). Conjugating
A tensor gamma_j by this matrix gives

    [[0,-iP gamma_j],[i gamma_j P,0]] = A tensor eta_j,
    eta_j=i gamma_j P.                                   (12)

Reflection of the scalar path reverses the coefficient list. Therefore

    M_R(u)=a(u) diag(1,0,1) tensor I
           +y(u)|Y><Y| tensor I
           +A tensor sum_j f_(m+1-j)(u) eta_j.            (13)

Equations (11)-(13) prove nonscalar, parity-compatible source resolvents
using only the local terms in the table. They are not C19's Jacobi weights.

## 5. Growing poles and full hopping rank

[proved] The Y channel alone has m distinct accessible energies

    E_l=Delta-2t cos(pi l/(m+1)), l=1,...,m.

Its endpoint weights are 2 sin^2(pi l/(m+1))/(m+1), all positive.
Thus geometry carries a number of poles growing with n on a fixed
alphabet. The nonscalar X,Z channel adds further poles.

An upper bound on the total number is 3m=6n. To see this, each tail
Schur complement on X,Z has determinant
(a^2-sum b_j^2)^N, where N=2^n. Taking the product over tails expresses
the full X,Z characteristic polynomial as r(z)^N for a rational r.
Because the result is a polynomial, a reduced denominator of r cannot
remain. Its degree is 2m. Hence this sector has at most 2m distinct
energies, and Y has m. At a root, an eigenvector vanishing in its entire
syndrome/data block would vanish at the next vertex by invertible
hopping, then at every vertex. Every eigenvalue is therefore accessible
as a matrix-valued source weight. Degeneracies and coincidences can lower
the count; no all-n claim that it is exactly 6n is required.

[proved] Define the Gram matrix

    Gf_(jk) = integral_0^infinity f_j(u) f_k(u) du.         (14)

It is strictly positive definite, not merely positive semidefinite.
The first appearance of Gamma_j in a root moment requires walking to
vertex j, using its potential once, and returning. Its first order is
2j-1. That unique shortest walk has coefficient g t^(2j-2). Thus

    f_j(u)=-g t^(2j-2) u^(2j-1)/(2j-1)! + higher powers.  (15)

A linear relation among the f_j is eliminated successively by its
u, u^3, ..., u^(2m-1) coefficients. They are linearly independent.
Integrability and analyticity then imply strict positivity in (14).
This argument proves full rank even when its smallest singular value
is too small to resolve reliably with double precision.

There are at least m+1 total accessible energies. Indeed, the m
independent coefficient functions have zero value at u=0, so they lie
in a space of dimension at most k-1 when there are k energies. This
bound is compatible with m+1<=k<=3m and with the numerical k=3m examples.

## 6. Second-order exchange

[proved] Write R for the m-by-m reversal matrix, R_(j,k)=delta_(j,m+1-k).
The entrywise contraction of (11) and (13), with the overall minus sign
from the positive resolvent, gives

    K2 = a0 I + sum_(j,k) J_(jk) eta_j tensor gamma_k,
    J = -2 R Gf,
    a0 = -integral_0^infinity [2a(u)^2+y(u)^2] du.          (16)

The factor two is sum_ab A_ab A_ab. Linear terms vanish because A has
zero diagonal, and the Y channel does not meet A. Hence (16) is exactly
Gaussian, not just its leading term in a second small coupling g.
Its odd operator Schmidt rank is m=2n by (14).

For a chain of prescribed cells, first order vanishes by (7). At second
order, different bonds have zero cross term: an outer cell is acted on
once, remains orthogonal to its ground code during the unperturbed
resolvent, and is killed by the final projection. Thus

    H_eff=lambda^2 sum_bonds K2+O(lambda^3).               (17)

This is a finite-cluster expansion. It does not state a uniform
thermodynamic error bound. Higher orders are not asserted to be Gaussian.

## 7. Why this still fails the requested target

[proved restriction] Full hopping rank does not establish n gapless
flavours. For (16),

    J^(-1) J^T = Gf^(-1) R Gf R.

This matrix is similar to
Gf^(-1/2)(R Gf R)Gf^(-1/2), which is positive definite. The zero-field
Majorana Bloch symbol, up to an overall convention-dependent nonzero
factor, is

    B(k)=i[J exp(ik)-J^T exp(-ik)].                        (18)

If B(k) has a zero, exp(2ik) is an eigenvalue of J^(-1)J^T. It must be
positive real. On the unit circle it is therefore 1, so k is 0 or pi.
There are no displaced Fermi points at zero field in this family.
Absence of eigenvalue 1 gives a gap.

For n=1 write Gf=[[a,b],[b,c]]. Strict positivity implies
(a+c)^2-4b^2>0. The determinant test for a zero of (18) reduces to

    (a-c)^2 cos^2(k)+[(a+c)^2-4b^2] sin^2(k)=0.            (19)

It is gapped whenever a differs from c. Adding the one-flavour chemical
potential changes the first squared term into the square of a shifted
multiple of cos(k); the strictly positive sin^2 coefficient remains.
It still cannot give a pair of displaced simple Fermi points.

[numerical] At the stated parameters,
a=0.00450765772747812, b=0.000301369017429837,
c=0.0000411764055123314. The difference a-c is about 0.00446648,
well above numerical errors. This example is gapped, not a c=1 chain.

[open] Separate source paths, different boundary mixing, or further
syndrome channels might change the form R Gf while preserving the
Clifford-vector closure. No such repair with independent hopping phases
and chemical potentials is supplied here. Adding arbitrary nonlocal
logical operators by hand would abandon the local rule table.

## 8. Numerical interpretation

[numerical] `clifford_cursor_checks.py` independently builds the physical
one-cursor Hamiltonian, with its controlled-Z bonds, then diagonalizes it
and contracts the two-cell resolvent. It does not construct H2 by starting
with (16). At n=1,2,3 the sector dimensions are 12,48,144, and there are
6,12,18 distinct accessible energies. The error after projection onto
identity and all eta_j tensor gamma_k is at most 5.0e-15. The source
resolvent span is separately checked at z=0.73 to below 5.6e-17.

The smallest hopping singular values are 4.19e-5, 1.55e-9, and 5.75e-14.
The last is poorly conditioned for a numerical rank claim. Equation (15)
is the all-n rank proof. Exact code degeneracy and parity follow from
(4)-(5), and are also checked directly for the elementary matrices.

Full finite-coupling ED includes fillers, multiple cursors, and all data
states. It compares the canonical effective matrix obtained from the
polar factor of the overlap with G to lambda^2 sum K2. Small eigenvalue
errors alone would not be that test. See RESULT.md and the JSON record
for the completed sizes and residuals. This ED tests the bond (16), which
fails the target, not an implementation of C19's intended bond.
