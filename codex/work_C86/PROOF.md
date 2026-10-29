# C86. Obstructions to replacing the fold

[Result] No design with at most ten letters and a certified member was
obtained. The new results are analytic obstructions to several proposed
mechanisms, a larger class with a restricted eleven-letter lower bound,
and exact rejection certificates for single-rail alternatives. They do
not establish a universal eleven-letter lower bound for route A.

## 1. The retained reduction and its scope

[Definition] This note uses the two-species, common-scalar-carrier
reduction of C75, C80, and C84. The species charges differ by one payload
charge. After removing the Jordan-Wigner strings, both species have the
same finite Hermitian positive matrix H. Conversion at vertex i changes
the species and annihilates or creates logical fermion i. A port creates
the opposite-charge pair in the same antisymmetric combination used in
C84. Its scalar right and left source vectors are u and v.

Let K(t)=exp(-tH), let P_i=|i><i| be the conversion projectors, and set

    R=(H tensor I+I tensor H)^(-1),
    x=R(u tensor v),
    U_i=(P_i tensor I)x,  V_i=(I tensor P_i)x,
    C_ij=U_i^dagger R V_j,
    A_ij=U_i^dagger R U_j+V_i^dagger R V_j.              (1)

[Inherited proved reduction] Up to a common positive scale, C84 gives

    T=-2 Re C,                 M=-2i Im A.              (2)

The omitted scale is alpha^4 g^2, with lambda^2 restored in physical
units. Re and Im in (2) are entrywise. Consequently T is real and M is
Hermitian. A common phase of the intercell hopping can also be allowed.
The Bloch pencil is

    h(k)=M+epsilon I+T exp(ik)+T^T exp(-ik).             (3)

[Scope] Equation (2) is a restriction on the mechanism, not a theorem
for arbitrary ports. Unequal scalar carriers, payload-dependent source
operators, extra conversion channels, and higher-order leading processes
require their own reduction. A two-site source here means an arbitrary
scalar vector on two carrier vertices. This is a relaxation of physical
port locality, not a claim to have implemented every such vector with a
nearest-neighbor microscopic bond.

## 2. A filter identity valid for complex carriers

[Proved] For a vector w define the real function

    f_(w,i)(t)=integral_0^t
       <w|K(s)|i><i|K(t-s)|w> ds.                       (4)

It is real because conjugation of the integrand is the substitution
s -> t-s. It need not be pointwise positive for a general complex w or H.
The identity needed below is

    Re C_ij = (1/2) integral_0^infinity
                    f_(u,i)(t) f_(v,j)(t) dt.          (5)

Here is the derivation. Expand the three resolvents in (1) as Laplace
integrals, with variables a,b,c>=0. The sums over the two unobserved
carrier vertices compose adjacent semigroups. This gives

    C_ij = integral da db dc
       <u|K(a)|i><i|K(b+c)|u>
       <v|K(a+b)|j><j|K(c)|v>.                         (6)

Put t=a+b+c. The integration domain in a,c is the triangle
a,c>=0, a+c<=t. The integrand becomes

    q_i(a,t) q_j(c,t)^*,
    q_i(a,t)=<u|K(a)|i><i|K(t-a)|u>,
    q_j(c,t)=<v|K(c)|j><j|K(t-c)|v>.

The substitution (a,c)->(t-a,t-c) maps this triangle to its complement
in the square [0,t]^2 and conjugates the integrand. Thus twice the real
part of the triangle integral equals the square integral. The latter
factorizes into the two integrals (4), proving (5). Positivity of H
ensures convergence of every integral used here.

[Consequence] Complex hopping phases do not invalidate this identity.
A proof restricted to real heat kernels would leave out the phase
suggestion in the task. Equation (5) includes it.

## 3. Reflection gives a signed Gram matrix

[Definition] Suppose there is a unitary involution P on the carrier
space satisfying

    P H P^dagger=H,  P u=v,  P v=u,
    P P_i P^dagger=P_(r(i)).                            (7)

The active conversion set is closed under the involution r. Let J be
its permutation matrix on logical flavors. P may include vertex gauge
phases. It is unitary, not antiunitary. No reality assumption on H is
made. Sources can have support on more than one vertex.

[Proved] Equation (7) implies
f_(v,j)=f_(u,r(j)). Define

    G_ij=(1/2) integral_0^infinity
                      f_(u,i)(t) f_(u,j)(t) dt.

Then G is real symmetric positive semidefinite and

    F:=Re C=GJ,
    Herm(F)=(GJ+JG)/2.                                 (8)

Use an orthonormal real basis with J=diag(I_p,-I_q). In this basis

    G = [[G_++, G_+-], [G_-+, G_--]],
    Herm(F) = [[G_++,0],[0,-G_--]].                     (9)

If p,q>0, this matrix cannot be positive definite or negative definite.
If the filters in (4) are linearly independent, G is positive definite
and the inertia in (9) is exactly (p,q). This conclusion does not depend
on diagonal masses, the shape of the reflected graph, or its complex
loop amplitudes.

[Proved] A common port phase cannot repair the obstruction. For
F_theta=exp(i theta)F,

    Herm(F_theta)=
      [[cos(theta)G_++, -i sin(theta)G_+-],
       [i sin(theta)G_-+, -cos(theta)G_--]].             (10)

When cos(theta) is nonzero, the two diagonal blocks have opposite signs
on the two parity subspaces, strictly so when G is positive definite.
When cos(theta)=0, each parity subspace has zero quadratic form.
Definiteness is impossible in both cases. For G>0 and cos(theta)>0,
the inertia is (p,q), by a Schur complement. For cos(theta)<0 it is
(q,p). At cos(theta)=0 there can also be zero eigenvalues.

[Proved special case] For a connected real Jacobi path with all vertices
active, the filters are independent. If vertex i is distance d_i from
endpoint u, the first nonzero term of K_(u,i)(t) is

    (product of the positive edge magnitudes) t^(d_i)/d_i!.

Substitution in (4) gives leading term

    (product of edge magnitudes)^2
                    t^(2d_i+1)/(2d_i+1)!.

The distances are distinct. Successively taking the smallest Taylor
power in any putative linear dependence proves that all coefficients
vanish. Thus for a reflected m-vertex path the inertia of Herm(F) is
(ceil(m/2),floor(m/2)), for every m. This extends the earlier two-flavor
obstruction to an unbounded family.

[Proved robustness] If G>0, put
gamma=min(lambda_min(G_++),lambda_min(G_--)). An added Hermitian
correction of norm less than gamma cannot remove either sign in (9).
More generally a Hermitian correction of rank less than q cannot make
(9) positive definite. Its kernel intersects the q-dimensional negative
subspace nontrivially. A correction making it negative definite needs
rank at least p. These are bounds on logical matrix corrections, not
on the number of physical port terms. One physical pair port can
generate a full-rank matrix, as C84 already demonstrates.

## 4. Reflection also forces reciprocal spectra

[Proved] Reflection obstructs more than C84's positivity entry test.
It also makes the spectrum of (3) reciprocal, even when H has flux.
Define S=(P tensor P) Swap on the product-carrier space. It commutes
with R and fixes u tensor v, hence fixes x. Moreover

    S U_i=V_(r(i)),  S V_i=U_(r(i)).

Taking matrix elements of R yields

    J C J=C^dagger,       J A J=A.

Since J is real, equations (2) imply

    J T J=T^T,           J M J=M,
    h(-k)=J h(k) J.                                    (11)

With a common phase theta in T, the reflection center is shifted to
k=-theta. Internal time-reversal breaking does not break (11) while
the unitary reflection (7) remains true. Rephasing the logical fermions
only conjugates the pencil and its reflection operator. It cannot
alter this spectral conclusion.

[Proved] If this reciprocal pencil has at least four distinct Fermi
points, its unordered pair-sum margin with replacement is zero. Every
nonfixed root k has partner -k and their sum is zero modulo 2pi. A
fixed root is 0 or pi and its repeated pair also sums to zero. Four
distinct roots always supply at least two distinct such unordered
pairs. A common phase changes the common sum to -2theta and does not
remove the collision. Thus a 2n-root certificate with n>=2 cannot pass
C84's pair test under (7), regardless of the onsite term M or a common
chemical potential.

[Proved] Restricting the logical space to one reflection parity does
not fix this test. The compressed F is +/-G_++ or +/-G_-- and is
symmetric. Its Bloch pencil depends on k only through cos(k), after a
possible common shift. Its pair-sum margin again vanishes whenever it
has at least four distinct roots. An arbitrary subspace mixing the
two parities is outside this statement.

[Proved] Co-locating the two sources, u=v, gives a related obstruction
without assuming reflection of H. Swap fixes x and interchanges U_i
and V_i, so C=C^dagger. Consequently Re C is symmetric and the spectrum
of (3) is even in k. Making both ends access the same scalar source
would supply a Gram matrix but would not by itself supply the required
nonresonance, even with a complex H.

[Proved application of the inherited tree observation] If H is supported
on a tree and each source occupies one vertex, gauge every edge real
by choosing vertex phases successively from a root. Gauge both species
in the same way. Their on-vertex conversion is unchanged. The two point
sources acquire only an overall phase in x, which cancels in C and A.
Both matrices can then be taken real, so M=0 and T is real. Hence
h(-k)=h(k)^* and the same pair-sum obstruction follows, without spatial
reflection or equal masses. C84 already noted that trees have no flux;
the locality argument in Section 7 explains which redesigned bonds
are forced to have this tree and point-source structure.

[Scope] These are statements about the retained Gaussian coefficient.
They are not claims that arbitrary microscopic corrections respect
(11). A scheme using a higher-order symmetry-breaking coefficient would
need new band and interaction bounds at that order. A phase that breaks
(7) is also outside the reflection theorem. For independent filters,
the indefinite hopping is additionally stable under the norm bound
following (10).

## 5. Arbitrary masses in the two-vertex single rail

[Proved] The old two-vertex no-fold obstruction does not need equal
masses. Gauge the edge real and scale its magnitude to one. Let

    H=[[a,-1],[-1,b]],  a,b>0,  z=ab>1,
    m=a+b,  s=a^2+b^2,  q=1/[2m(z-1)].

In product order 00,01,10,11 the inverse in (1) is

    R=q [[mb-1, b, b, 1],
         [b, 2z-1, 1, a],
         [b, 1, 2z-1, a],
         [1, a, a, ma-1]].                             (12)

Multiplication by H tensor I+I tensor H verifies (12). With u=e_1,
v=e_0, x=q(b,1,2z-1,a)^T. Restricting rows and columns as in (1) gives

    C_00/q^3=b^2(b^2+3z-1)+2z-1,
    C_11/q^3=a^2(a^2+3z-1)+2z-1,
    C_01/q^3=s+3z-1,
    C_10/q^3=(2z-1)[s+(2z-1)^2]+z.                    (13)

Writing S=Herm(C), direct expansion gives

    det S = -[s^2+(8z^2+z-1)s
                 +16z^4-16z^3+18z^2-6z]
                 /[64m^6(z-1)^4] < 0.                (14)

Every bracketed term is positive when z>1. For the constant term,
use 16z^3(z-1)+6z(3z-1)>0. Also

    C_10-C_01=[s+2z(2z-1)]/[4m^3(z-1)^2]>0.          (15)

For any common port phase,

    det Herm(exp(i theta)C)
       = cos(theta)^2 det S
          -sin(theta)^2 (C_10-C_01)^2/4 <0.            (16)

For edge magnitude t>0, substitute a/t,b/t in these expressions.
C acquires t^-3 and its determinant t^-6. An edge phase is removed
by a vertex gauge and cannot change the conclusion. The disconnected
t=0 limit is excluded by the connected-rail hypothesis and is not a
source of propagating flavor channels.

[Exact finite checks] The direct rational inverse and equations
(14)-(15) agree at all twelve pairs a in {3/2,2,4,9},
b in {1,3,5}. The proof is the displayed identity for all allowed a,b,
not an extrapolation from those checks.

[Scope] This theorem concerns a two-vertex carrier with both vertices
active. It is not asserted for every principal two-flavor block of a
longer carrier. Longer asymmetric rails are tested separately below.

## 6. Single species and charge-changing conversion

[Proved, restricted] Suppose one excited species is used, its particle
number is preserved by the intra-cell conversion, and its charge is
independent of position and the neighboring data bit. Let the payload
charge increment be delta!=0. A conversion that changes one bit and
leaves this carrier charge unchanged changes the cell charge by
plus or minus delta. It cannot occur in a charge-conserving H_0.
Moving the head as part of the conversion does not alter this equation.

A charge-neutral exchange of two data bits is allowed, but it is not
the linear f_i or f_i^dagger conversion in (1). If all intra-cell
operations and a payload-independent port preserve each cell's payload
number, the full model has an independent U(1) for each cell. An
intercell f_left^dagger f_right term then vanishes to every order.

[Proved, restricted] If the carrier's data bit itself supplies the
charge change, the carrier must have two orthogonal charge eigenstates.
Naming both states one species does not reduce the local Hilbert-space
dimension. In a faithful one-token product representation this still
costs at least two carrier letters per transport type. If the
distinction is stored in another site, the encoding has changed and
must be audited there. C84's seven-letter charge count cannot be
applied unchanged to a port that explicitly consumes a payload bit.

[Proved alternative obstruction] Allowing conversion between a charged
carrier and the vacuum can instead lead to an enlarged quadratic
hybridization model. There is an independent bound on that proposed
escape. For any finite-dimensional quadratic Bloch Hamiltonian

    H_full(z)=A+z B+z^(-1) B^dagger,  rank B=r,          (17)

the Laurent polynomial det H_full(z) has powers only from -r to r.
To prove it, expand the determinant multilinearly in its columns. A
term using more than r columns from B vanishes. The same holds for
B^dagger. Each surviving term has exponent p-q with p,q<=r, hence
absolute exponent at most r. Multiplication by z^r therefore gives
a polynomial of degree at most 2r.

Unless the determinant is identically zero, there are at most 2r
distinct zero-energy crossings on the unit circle. The identically
zero case has a zero mode at every momentum and cannot be a certificate
of finitely many simple roots. Thus n massless Dirac pairs of this
quadratic type require r>=n. One independent quadratic transfer channel
allows at most one pair; two allow at most two.

[Proved] Integrating out gapped quadratic auxiliaries does not evade
this count. The Schur determinant identity is

    det H_full(z)=det H_aux(z) det H_eff(z).

If the auxiliary block is nonsingular on the unit circle, its
elimination creates no additional zero there. The frequency-dependent
effective action has the same zero-energy root count.

[Scope] Equation (17) does not apply to C84's interacting auxiliary
pair process. It also does not turn a two-site microscopic port into
a rank-two B by definition. It applies when a proposed single-species
replacement actually has the quadratic form (17) with r specified.
Nonquadratic data-dependent ports remain open.

## 7. A broader, still restricted alphabet lower bound

[Proved, restricted] Consider the following replacement class. It
retains the independent binary L and R backgrounds, a frozen on-site
marker with neutral and opposite-charge source states, and the
charge equations of C84. Carrier transport types are encoded by
one on-site token, with matched A and B states for every type.
The g=0 scalar motion is the same for the two species after the
payload gauge, and leaves the payload unchanged in that gauge.
There is no contextual storage of an unmatched carrier type. The
physical Hamiltonian is nearest-neighbor. Each cell is a contiguous
block ending at its marker, and a port creates one token in each of
the two incident vacuum cells without an additional payload operation.
It seeks the same leading Gaussian certificate.

There must be four background letters and at least three marker
letters. The charge equations are

    q(background)=u or u+delta,
    q(carrier A,B)=u-delta/2 or u+delta/2,
    q(marker)=v, v-delta/2, or v+delta/2,
    delta!=0.

Frozen marker roles keep the marker set disjoint from nonmarkers.
The carrier charges keep each carrier letter disjoint from the
backgrounds and from the other species. If r is the number of
on-site carrier transport types, the count is therefore

    a>=4+3+2r=7+2r.                                   (18)

For a<=10 there is room for only r=1. At fixed payload, each data
position then has one carrier vertex of a given species, together
with its one terminal charged marker vertex. A two-site microscopic
term can move this single token only to a neighboring physical
position. It cannot distinguish two histories ending at the same
position without another state in the stated product representation.
Thus its scalar graph is a path or a disconnected subgraph of a path.
Parallel terms on one edge combine into one complex matrix element
and do not make a flux loop.

The physical intercell bond touches only the last site of the left
cell and the first site of the right cell. Starting from zero tokens,
its two created tokens must occupy these two sites. Since each has
only one excited letter of each species, the scalar sources are point
vectors. The superposition over vacuum clock positions changes their
coefficients, not the excited-token positions. A source coherently
supported on two distinct carrier positions in the same cell cannot
be obtained this way. This is a locality conclusion for this product
cell class, not a definition of all possible two-site ports.

The tree and point-source result in Section 4 now forces reciprocal
spectra, so a certificate with n>=2 and positive pair-sum margin is
impossible. Since a>=9 in this class and the task needs n>=a, this
excludes every desired member at a<=10. Hence

    a>=11                                             (19)

in this class. This does not assume the express fold, reflection,
equal masses, or its mass schedule. An extra unmatched tenth letter
cannot supply a second common charged transport type under the stated
representation assumptions. Adding a fourth neutral marker state
does not add a charged carrier vertex. A new charged marker vertex
that is used by both common species would require a matched pair of
new marker states, which also costs two letters.

[Proved extension] Even if the point-source locality restriction is
relaxed and arbitrary extended scalar sources are granted, Sections
3 and 4 exclude every reflected replacement satisfying (7). This
includes complex reflection-invariant loops and reflected two-vertex
sources. If there are at least four distinct roots, it is the spectral
reciprocity itself that fails the pair test, not only the sign of the
hopping. Thus the suggested reflected extension does not defeat this
obstruction within the retained reduction.

[Scope] Compared with C84's faithful-fold count, (19) allows the
transport mechanism to change and derives the one-rail obstruction
from locality and reciprocity. It does not say that all transport
requires two rails. Unequal species propagators with controlled
onsite terms, a data-dependent port, contextual encoding, a different
cell decomposition, or an interaction-generated leading quadratic
term can evade the class. Asymmetric extended sources are outside
the reflected extension and are only tested at the finite points in
Section 8. No universal bound above C84's seven charge-role letters
is claimed after dropping those hypotheses.

## 8. Exact rejection certificates using C84's first stage

[Proved, computer-assisted] `certify_obstructions.py` adapts the
floor-rounded product-resolvent stage of C84. The single-rail matrix
is a path with n+1 vertices, diagonal 4 except the last diagonal
4+D, nearest-neighbor entries -1, endpoint sources, and active
vertices 0,...,n-1. The last vertex is inactive, as a marker would be.
All twenty cases

    n in {2,3,6,9,10},  D in {0,1,10,1000}

have a strictly negative principal determinant of Herm(Re C), using
active indices 0,n-1. These are exact finite exclusions of the proposed
scalar matrices, including members large enough to have been relevant
to a nine- or ten-letter attempt. They are not full microscopic bonds.

[Proved enclosure] Use scale S=2^256 and N=256+24 Jacobi iterations.
Every product-graph diagonal is at least 8 and the sum of off-diagonal
magnitudes is at most 4. Thus the Jacobi contraction in infinity norm
is at most 1/2, and the inverse norm is at most 1/4. Floor-rounding
contributes less than two integer units in a solve. The remaining
iteration tail is below one unit. The first source solution has
component error at most 3, and the second solution error at most 4,
including the propagated first-source error.

For nonnegative sources the integer iterates are lower bounds. The
upper product uses (x+3)(y+4). For signed real and imaginary components
use instead the absolute radius

    4(|x_real|+|x_imag|)+3(|y_real|+|y_imag|)+24.       (20)

Equation (20) follows by bounding each product error by
4|x|+3|y|+12. All sums and final comparisons use integers or Fractions.

[Proved, computer-assisted] A second exact scan uses the uniform
path with n=9 or 10, active vertices 0,...,n-1, and sources

    v=e_0+c e_1,  u=e_n+d e_(n-1),
    c,d in {0,1,-1,i,-i,(1+i)/2,(1-i)/2}.              (21)

All 98 cases have a strictly negative endpoint principal determinant.
For an interval S_ij in [s_ij-e_ij,s_ij+e_ij], the determinant upper
bound is the largest of the four endpoint products for the diagonal
entries, minus max(|s_01|-e_01,0)^2. It is strictly negative in every
case. This handles signed sources without an entrywise-positivity
assumption.

[Proved] Every real two-by-two principal block S=Herm(F) with
det S<0 also rejects every common phase of F, since

    det Herm(exp(i theta)F)_block
       =cos(theta)^2 det S
         -sin(theta)^2(F_01-F_10)^2/4 <=0.             (22)

At a possible equality the block is singular and still cannot be
definite. Thus the certificates reject a common chiral port phase,
not just the phase zero used in their computation.

[Proved quantitative stability] If a certified principal determinant
is at most -d<0 and its norm is at most B, both signs have magnitude
at least d/B. Each certificate stores this positive rational bound.
A Hermitian perturbation with norm below d/B preserves a positive
and a negative Rayleigh direction at phase zero. This specifies how
small finite-g or finite-K corrections would have to be to preserve
the rejection; no finite microscopic schedule is asserted.

[Rejected checker, corrected] The first two-site-source run stopped
at an assertion that the diagonal and off-diagonal lower endpoints
were positive. This need not hold for signed or complex sources.
The corrected acceptance is the general interval determinant bound
above. `rejected_checks.json` retains the failed assumption. No
conclusion uses the failed test.

[Scope of pipeline use] These candidates fail before C84's positive
LDL congruence and root-bracketing stage. Running its positive-metric
whitening on an indefinite matrix would be invalid. No root count,
positive pair margin, finite-K schedule, or new MEMBER file is
reported. Indefinite hopping alone does not exclude all possible
multiband spectra. The stronger pair obstruction of Section 4 applies
to the reflected class; outside it the finite certificates only reject
the requested positive-Hermitian mechanism.

## 9. Consistency checks and remaining scope

[Numerical] `consistency_checks.py` independently integrates (4)-(5)
using the spectral heat kernel for complex reflected graphs of sizes
4,5,6 and extended complex sources. Relative errors are below 1.7e-15.
The separate checks of J T J=T^T, J M J=M and (11) are below 2.5e-15.
A non-reflected five-vertex complex graph with equal sources gives
symmetric Re C to relative error below 6.2e-17. These are numerical
consistency checks, not replacements for Sections 2 through 4.

[Numerical] `explore.py` also checks 24 single-rail cases, the 98
extended-source cases, and reflected complex graphs with 4,5,6,10
vertices. Some tiny eigenvalues of the Gram matrix are below floating
precision. In particular its ten-vertex numerical minimum is about
-2.5e-21. It is not used to assert failure of Gram positivity. The
analytic integral proves semidefiniteness; the exact interval tests
prove the claimed finite rejections.

[Open] No all-parameter theorem for an asymmetric long single rail is
proved here. No completed one-species data-dependent port or ten-letter
bond was obtained. A successful next construction must supply both
the positive-Hermitian or alternative crossing argument and a mechanism
that avoids reciprocal pair sums, without merely assigning a common
phase to all hoppings. The inherited nonzero-port defect comparison
and infrared transfer would still be needed after any finite band
certificate. No central charge was measured or inferred for a new model.
