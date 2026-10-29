# C87. Extensions beyond a single product token

[Result] No bond with at most ten letters and a certified member n>=a
was obtained. This note proves obstructions for specified extensions,
including an arbitrary-number-of-tokens phase theorem and a classification
of diagonal port dressings that preserve a Gaussian coefficient. It also
gives an interacting port whose second-order coefficient is Gaussian.
That example is an abstract finite gadget, not a microscopic candidate.

[Convention] Gaussian operator means scalar plus number-conserving
fermion quadratic on the full logical Fock space. All cells below are
open intervals between markers. A theorem about an open cell need not
hold for tokens winding around a physical ring. Charge means an additive
on-site conserved charge unless another meaning is stated.

[Inherited] C84 has an eleven-letter finite certificate at n=11, with
22 simple Fermi roots and a positive margin between every two distinct
unordered root pairs, including repeated roots. Its Gaussian identity
is at order lambda^2 g^2. The full second-order coefficient in lambda
at finite g generally has interactions. Its infrared conclusion remains
conditional. Neither distinction is changed here.

## 1. Homogeneous token exchanges cannot generate a flux

[Definition] Fix a multiset of letters on an open path of b sites.
One letter may represent an empty site. The off-diagonal moves are
nearest-neighbor exchanges of unequal letters. For s<t, the directed
move st -> ts has matrix element

    -t_e exp(i chi_st),  t_e>0.                         (1)

The phase chi_st depends only on the two exchanged species, not on the
bond position or any surrounding letters. The magnitude t_e may vary.
The reverse move is the Hermitian conjugate. Arbitrary real diagonal
energies, including interactions, are allowed. Some edges or vertices
may be removed by constraints.

[Proved] Every phase in (1) is removable by a diagonal unitary on the
configuration space. Choose an ordering of the species and define

    I_st(w)=number of pairs x<y with w_x=t and w_y=s,
    Phi(w)=sum_(s<t) chi_st I_st(w),
    D|w>=exp(i Phi(w))|w>.                             (2)

An adjacent exchange st -> ts changes I_st by one. The counts involving
any third letter do not change, because that letter is either to the
left of both exchanged sites or to the right of both. Counts of other
occurrences of s or t also do not change in total. Therefore

    Phi(ts)-Phi(st)=chi_st,
    (D^dagger H D)_(ts,st)=-t_e.                       (3)

The diagonal is unchanged. This proves the assertion on the full graph
and on every constrained subgraph. Every directed cycle has zero phase
after factoring out the conventional minus sign per edge. Equivalently,
all complex phases are gauge removable. The configuration graph need
not be a tree. It can have many independent cycles.

[Consequence, restricted] Adding two or more homogeneous hopping tokens
and allowing their adjacent exchange does not by itself replace C84's
flux loop. The statement includes arbitrary token numbers and species,
with fixed multiplicities. Diagonal collision energies cannot change it.
It does not include species-changing reactions, context-dependent
exchange phases, extended coherent sources, or non-diagonal conversion
operators whose phases fail to transform compatibly with (2).

[Proved] Without exchange, the order of distinguishable hard-core tokens
on a path is conserved. For two tokens the two orders are disconnected.
If they represent the two rails, hopping cannot turn one rail into the
other. Allowing a homogeneous adjacent exchange connects the orders but
still leaves the phases removable by (2). For one token on each of two
independent open rails, the product of the two one-particle gauges does
the same job, including in the presence of real diagonal interactions.

[Proved sharper two-token test] Let a and b denote two token species.
Allow right-hop phases theta_a(j), theta_b(j) across bond j, and an
exchange phase chi(j) for ab -> ba there. Remove the single-particle
phases by their cumulative position gauges. The remaining exchange
phase is

    eta(j)=chi(j)-theta_a(j)+theta_b(j).                (4)

If eta(j) is constant, an additional phase eta times the indicator
that a is to the right of b removes every exchange phase. More generally
the scalar Hamiltonian can be made real if eta(j) is constant modulo pi.
When all the relevant hops and exchanges are present, this condition
is also necessary. The six-step cycle

    (j,j+1) -> (j+1,j) -> (j+2,j) -> (j+2,j+1)
      -> (j+1,j+2) -> (j,j+2) -> (j,j+1)               (5)

has phase eta(j)-eta(j+1). A real matrix has only real cycle products.
Thus a difference not in pi Z is a genuine obstruction to a real gauge.
Sufficiency follows from the explicit order gauge, with real signs left
on the exchange edges. A gauge making every edge negative requires the
stronger condition modulo 2pi.

[Scope] A nonuniform eta is a concrete way around the theorem. It needs
position dependence or contextual states in a one-site invariant bond.
No implementation of that dependence with at most ten letters is given.

## 2. Consequence for C84's product-resolvent coefficient

[Definition] Suppose a proposed multi-token cell still admits the common
scalar reduction, with H>0, two point sources u,v, and real diagonal
conversion projectors P_i in the gauge where H is real. Their ranks may
exceed one. Set

    R=(H tensor I+I tensor H)^(-1), x=R(u tensor v),
    U_i=(P_i tensor I)x, V_i=(I tensor P_i)x,
    C_ij=U_i^dagger R V_j,
    A_ij=U_i^dagger R U_j+V_i^dagger R V_j,
    T=-2 Re C, M=-2i Im A.                            (6)

The common positive normalization from C84 is suppressed. Equation (6)
is an explicit assumption on the conversion and port reduction. It is
not automatic for every two-token microscopic Hamiltonian.

[Proved] If H has the gauge (2), point sources acquire only individual
scalar phases. Their product is a common phase of x and cancels from
C and A. Diagonal P_i commute with D. Thus every remaining quantity
in (6) can be calculated using real matrices, and M=0, T real. The Bloch
pencil, including a common field, obeys

    h(k)=epsilon I+T exp(ik)+T^T exp(-ik),
    h(-k)=h(k)^*.                                     (7)

The spectrum is reciprocal even if T is nonsymmetric, the token
potential is asymmetric, or the configuration graph has cycles. A common
hopping phase moves the reflection center in k and does not remove it.

[Proved] A reciprocal spectrum with at least four distinct roots has a
zero unordered pair-sum margin. Each nonfixed root k has partner -k,
and each fixed root, 0 or pi, has a repeated pair with sum zero modulo
2pi. Four distinct roots supply at least two different unordered pairs
with that same sum. Consequently every desired n>=2 member of (6)-(7)
fails C84's pair-sum stage. This is stronger than a failure of positive
Hermitian hopping and is independent of its inertia.

[Proved, computer-assisted finite check] The six configurations of two
distinguishable tokens on three sites give a 36-dimensional product
resolvent. In the real gauge take hopping and exchange amplitudes -1,
diagonal 10+2a+3b+|a-b| for positions (a,b), sources (2,1) and (0,1),
and P_i the rank-two projector onto a=i. The diagonal is greater than
the off-diagonal row sum, so H>0. An exact Fraction inverse verifies
the inverse identity and computes the rational C stored in checks.json.
A is real symmetric, so M vanishes exactly. This tests a carrier with
a configuration-space cycle and rank-two conversion projectors, not a
single-particle path. It is a scalar rejection certificate, not a bond.

[Scope] In particular, a species-changing conversion can acquire a
spectator-dependent phase under (2). If it does, the real-diagonal
conversion hypothesis of (6) has failed. Such a design must be reduced
afresh. Applying (7) to it without checking conversion would be wrong.

## 3. Storing data only in a bounded number of token positions

[Proved dimension bound] With k indistinguishable hard-core tokens on b
sites, q internal states per token, and a marker space of dimension d_M,
the available dimension is at most

    d_M q^k binomial(b,k).

An injective encoding of n qubits therefore requires

    2^n <= d_M q^k binomial(b,k).                       (8)

Allowing distinguishable tokens gives the more permissive upper bound
d_M q^k b^k. Hence fixed k,q,d_M imply n=O(log b), not n=b-O(1).
This argument allows arbitrary coherent encoding into the position
space, not just assigning bit strings to individual configurations.

[Proved finite consequences] For two distinguishable tokens on the same
chain, with no extra data space, the dimension is b(b-1). For identical
tokens it is binomial(b,2). For separate rails of b sites each it is b^2.

| n | Distinguishable, same chain | Identical, same chain | Sites per separate rail |
| :-- | --: | --: | --: |
| 9 | 24 | 33 | 23 |
| 10 | 33 | 46 | 32 |
| 11 | 46 | 65 | 46 |

These are necessary sizes, not constructions. For example, C84's
geometry n=10, b=12 would have only 132 two-token configurations, below
1024. Making b much larger avoids this dimension obstruction.

[Proved charge obstruction] Suppose the data configurations have a fixed
multiset of on-site charge eigenletters, including fixed token species,
an otherwise uniform background, and a fixed marker state. Then

    Q|code>=(sum_s k_s q_s)|code>.                     (9)

The additive physical charge is scalar on the whole encoded space,
irrespective of token positions. It cannot equal q_0+delta sum_i n_i
on a full n-mode Fock code for delta nonzero and n>=1. Thus this form of
position-only encoding cannot inherit C84's nontrivial logical U(1)
from the retained additive charge. Increasing b does not fix (9).
Variable token number, varying charged internal species, charged marker
data, or a different exact conservation law can evade it.

[Proved ground-space obstruction] For a connected component of a finite
configuration graph whose off-diagonal entries can all be made strictly
negative, choose c greater than every diagonal entry. Then cI-H is
nonnegative and irreducible. Its Perron eigenvector is unique, so H has
a unique ground state in that component. With two distinguishable tokens
and all ordinary vacancy hops, the two fixed-order components are
connected. Without exchange there are at most two ground states in
total, and with exchanges connecting the orders there is one. Arbitrary
real diagonal interactions do not change this conclusion. Two separate
connected rails also give a connected product graph and one ground
state. These cannot furnish an exactly degenerate 2^n-dimensional
vacuum code for n>=2 using positions alone.

[Scope] This last argument concerns exact ground degeneracy. It does
not exclude engineering a narrow low-energy band with very small hops,
many disconnected wells, or extra conserved data. Such a band would
need its own separation and perturbation budget. The dimension bound,
the charge obstruction, and the ground-space obstruction are separate
statements with separate hypotheses.

## 4. A complete test for a second-order Gaussian operator

[Proved] Let the code energy be zero, let QH_0Q=A>0, let PVP=0, and
write S=QVP in a fixed orthonormal code basis. Eliminating Q gives

    H_eff(E)=-lambda^2 S^dagger(A+lambda QVQ-E)^(-1)S,
    K_2=-S^dagger A^(-1)S.                            (10)

For multiple orthogonal excited channels this is the sum over channels.
No claim of Gaussianity follows from the negative Gram form alone.

[Proved exact criterion] On m logical fermions, suppose K_2 conserves
number. Let |i>=f_i^dagger|0>, with a fixed mode ordering, and define

    kappa=<0|K_2|0>,
    h_ij=<i|K_2|j>-kappa delta_ij,
    K_G=kappa I+sum_ij h_ij f_i^dagger f_j.             (11)

Then K_2 is Gaussian if and only if K_2=K_G on the entire Fock space.
Necessity follows by evaluating a quadratic on vacuum and one-particle
states. These fix every coefficient uniquely. Sufficiency is its
displayed form. Computing only vacuum and one-particle matrix elements
without checking the rest cannot test Gaussianity.

[Proved robust witnesses] For a Gaussian operator, the hopping matrix
element divided by its fermionic sign is independent of all spectator
occupations. Two such values a,b give distance at least |a-b|/2 from
every Gaussian in operator norm. Indeed both matrix-element errors are
at most the operator norm, and their difference is a-b. Also every
mixed second finite difference of its diagonal energies vanishes. If
that difference is d on four basis states, the operator-norm distance
from any Gaussian is at least |d|/4. These bounds do not assume (11) is
the best approximation in operator norm. Pairing terms, if allowed,
do not alter these particular witnesses.

## 5. Classification of diagonal phase dressings

[Definition] Let K=kappa I+sum_ij h_ij f_i^dagger f_j be Hermitian and
number conserving on the full m-mode Fock space. The graph whose edges
are nonzero off-diagonal h_ij is connected. Let U be diagonal in the
occupation basis, U|X>=u_X|X>, with |u_X|=1. Consider a port change
S -> SU that leaves A unchanged. It gives K -> U^dagger K U in (10).

[Proved] U^dagger K U is Gaussian if and only if

    U=exp(i F(N)) exp(i sum_j theta_j n_j), N=sum_j n_j, (12)

where F is any function on {0,...,m}. Thus its only effect on K is an
ordinary one-particle gauge. The possibly nonlinear F(N) commutes with K.

Proof. A diagonal U preserves every diagonal entry of K, including the
vacuum energy. In the one-particle sector it changes h_ij to
u_{i}^* u_{j} h_ij. Suppose the full transformed operator is Gaussian.
For any occupied j and empty i joined by an edge, put
Y=X minus {j} plus {i}. Dividing by the nonzero h_ij and the common
fermionic sign gives

    u_Y^* u_X=u_i^* u_j.                               (13)

Define v_X=u_X/product_(j in X) u_j. Equation (13) says v_Y=v_X.
The configuration graph of r indistinguishable exclusion particles on
a connected graph is connected for every 0<=r<=m. One proof restricts
to a spanning tree, fixes occupancy at a leaf by sliding a particle or
hole along the path to it, and removes that leaf. Repeating matches any
target r-subset. Hence v_X depends only on r=|X|, which proves (12).
Conversely F(N) commutes with K and the remaining factor transforms
f_i^dagger f_j by the one-particle phases in (13). It is Gaussian.

[Proved implication for identical cell dressings] If U is a product of
identical cell-diagonal unitaries and the intercell hopping connects
their mode graphs, Gaussianity on the full two-cell Fock space permits
only linear on-cell phases, up to a scalar. To see this, use (12) and
remove the one-mode phases. Write the residual cell phases as g(X),
normalized to g(empty)=1. Two-cell configurations imply
g(X)g(Y)=v_(|X|+|Y|). Taking Y empty shows g depends only on |X|.
The one-particle normalization gives g(1)=1. Taking |Y|=1 gives
g(r+1)=g(r) for r<m, so every g(r)=1. This includes the full-cell
occupation, reached from r=m-1. The global result then fixes all totals.
Thus a cell charge-dependent nonlinear phase cannot change the bands
while preserving Gaussian intercell hopping in this dressing class.

[Proved quantitative example] Take four modes, with 0,1 in the left
cell and 2,3 in the right cell, and

    L_j=f_j+f_(j+1), j=0,1,2,
    K=-sum_j L_j^dagger L_j,
    U=(-1)^(n_0 n_1).                                 (14)

The hopping 1 -> 2 has coefficient -1 without spectator 0 and +1 with
it, after dividing by the fermionic sign. Therefore U^dagger K U is
at operator-norm distance at least one from every Gaussian. A correction
of norm less than one cannot make it Gaussian. This failure occurs in
the leading second-order coefficient, not at an unbudgeted higher order.

[Scope] The theorem classifies input phase dressings of one fixed
effective operator. It does not classify an arbitrary sum of excited
channels, an occupation-dependent resolvent, changes of code, or a
new port acting differently on the excited space.

## 6. Amplitude dressings and a cancellation that survives

[Proved restricted amplitude obstruction] Suppose the undressed
second-order operator at g=0 is -kappa I, kappa>0, and a single input
dressing factorizes as D_L tensor D_R, with D_L,D_R diagonal. Put
a(X)=|D_L(X)|^2 and b(Y)=|D_R(Y)|^2. Its diagonal energy is

    E(X,Y)=-kappa a(X)b(Y).                            (15)

For any independent changes X -> X', Y -> Y', its mixed difference is

    -kappa [a(X')-a(X)] [b(Y')-b(Y)].                  (16)

Every Gaussian diagonal energy is affine in occupations, so (16) must
vanish. Unless one factor is identically zero, at least one of a,b must
be constant; the remaining factor must itself be affine in occupation
bits. In particular equal nonconstant factors on the two cells fail.
This is an obstruction before the small g^2 hopping is even generated.

[Proved example] D=(1+n_0)(1+n_2) gives

    K_2=-D^dagger D=-1-3n_0-3n_2-9n_0 n_2.            (17)

Its mixed difference is -9, so its distance from every Gaussian is at
least 9/4. This does not rely on a numerical fit or on discarding a
scalar energy shift.

[Proved counterexample to a broader obstruction] Interacting ports can
have an exactly Gaussian second-order operator. Define

    L_1=f_0+f_1, L_2=f_2+f_3, P=n_2, Q=n_0,
    S_1=P L_1, S_2=(1-P)L_1,
    S_3=Q L_2, S_4=(1-Q)L_2.                          (18)

Give these four channels mutually orthogonal excited spaces with the
same gap one. Since P^2=P and Q^2=Q,

    K_2=-sum_a S_a^dagger S_a
       =-L_1^dagger L_1-L_2^dagger L_2.               (19)

Each individual channel in (18) has a non-Gaussian Gram contribution,
but the sum is quadratic. The cancellation is an operator identity on
all 16 code states. Each source lowers logical number by one; assigning
one additional auxiliary charge to its excited space makes the finite
gadget charge conserving. This assignment does not specify an on-site
microscopic alphabet.

[Scope] Orthogonal equal-gap channels are essential to this displayed
identity. Unequal inverse gaps leave a term proportional to L^dagger P L.
More general coordinated channels might also cancel, so the failure of
a single diagonal dressing is not a no-go for all data-dependent ports.
This example has no spatial cell grammar, no length selector, no literal
nearest-neighbor bond, and no alphabet count. It is not a <=10-letter
candidate, and its four logical modes do not meet the requested n>=a.

[Proved and numerical second-order check] For any stacked source S in
these examples, the full finite Hamiltonian is

    H(lambda)=[[0,lambda S^dagger],[lambda S,I]].       (20)

If B=S^dagger S has eigenvalue beta, its low eigenvalue is

    E_beta(lambda)=(1-sqrt(1+4lambda^2 beta))/2
       =-lambda^2 beta+lambda^4 beta^2+O(lambda^6).    (21)

This follows by diagonalizing the two-dimensional block with entries
0,lambda sqrt(beta),lambda sqrt(beta),1. Zero singular values give zero
code energies. Thus direct diagonalization tests all code levels, not
only a one-particle sector. For (19), B=2n_a+2n_b in two orthogonal
normalized modes, and

    B^2=4n_a+4n_b+8n_a n_b.                           (22)

The second-order Gaussian identity survives, but the fourth-order term
is interacting. No all-orders Gaussian claim is made.

## 7. Frozen markers that also store charged data

[Definition] Keep a frozen on-site marker subspace, disjoint from the
nonmarker subspace. Suppose it stores k binary logical bits with charge
increment delta!=0. Its neutral code has charge v+r delta and
multiplicity m_r=binomial(k,r). Assume all these marker states are
available at the port with the same other local input, and that the
port retains their information on the marker while producing either
charged marker source.
Thus each of its two marker maps is injective on the neutral code and
shifts its charge by +delta/2 or -delta/2, respectively. Allow the two
image spaces to overlap whenever their charges coincide. This is more
permissive than retaining separately labeled A and B states.

[Proved sharp charge-space bound] The neutral dimension is 2^k. At
excited charge v+(r+1/2)delta the source from neutral r and the other
source from neutral r+1 require dimensions m_r and m_(r+1). Injectivity
therefore requires at least max(m_r,m_(r+1)) states at that charge.
Neutral and excited charge lattices are disjoint. With m_r=0 outside
0<=r<=k,

    d_M >= 2^k + sum_(r=-1)^k max(m_r,m_(r+1))
        = 2^(k+1)+binomial(k,floor(k/2)).              (23)

For the equality, use max(x,y)=(x+y+|x-y|)/2. The sum of x+y is 2^(k+1).
The binomial multiplicities are unimodal, so their total variation from
zero back to zero is twice their maximum. This proves (23).
The dimension bound is attainable as a charge-space statement by
nesting the smaller source image inside the larger at each excited
charge. It does not prove compatibility with carrier propagation.

[Proved consequences] One marker bit requires at least five marker
letters, two neutral and three charged. In units delta=2, the charges
can be 0,2 and -1,1,3. The two source images share the charge-one state.
If the two images must remain orthogonal for all payloads, the stronger
bound is 3 times 2^k, hence six for one bit. Two bits already require
at least ten marker letters under the permissive bound (23).

[Proved restricted alphabet implication] Retaining four binary L/R
backgrounds and at least two charge-distinct nonmarker carrier states,
putting one charged bit in the marker gives

    a >= 4+2+5=11.                                    (24)

The marker remains a disjoint subspace even if some charges overlap
with nonmarkers. Equation (24) leaves the single-token class by allowing
marker payload and overlapping charged source images, but retains the
specified background and port assumptions. It does not exclude reducing
the background alphabet at the same time, non-injective ports that move
the marker bit elsewhere, or contextual markers lacking a frozen
on-site projector.

[Consistency with the parent] C80's four-state marker is not a
counterexample to (23). Its extra E_1 stores a displaced bit when the
token is in the interior. The terminal vacuum source uses E_0, with
the payload relocated into the L prefix. Its pair port is not injective
on an independently variable two-dimensional marker payload at fixed
other input. Contextual storage of that kind is explicitly outside the
new five-state count. C84 removes E_1 by fixing the displaced boundary
bit. No stronger count is being imposed retroactively on either parent.

## 8. Checks, rejected checker, and remaining work

[Proved, computer-assisted] checks.py verifies (2) by integer phase
arithmetic for three tokens plus holes on paths of 3 through 7 sites.
The largest graph has 210 vertices, 450 edges and 241 independent cycles.
It also verifies the two-token gauge and connected-component count for
lengths 3 through 9, with and without exchange. The deliberately
nonuniform exchange example has cycle phase pi/2, demonstrating the
boundary of the theorem. The 36-dimensional Fraction inverse described
after (7) checks the scalar reduction without floating acceptance tests.

[Proved, computer-assisted] All 256 diagonal sign unitaries on three
modes are checked against (11) for connected path hopping. Exactly 64
preserve Gaussianity, and exactly those 64 have form (12). The phase
and amplitude obstructions (14),(17), the full identity (19), and the
non-Gaussian fourth-order coefficient are checked with integer matrices.
Position counts and marker multiplicities through k=8 are exact integers.

[Numerical, sizes stated] Full diagonalization of the two-token carrier
matrices uses dimensions 6 through 72. Each predicted connected component
has a unique ground state, with all eigenpair matrix residuals below
8.6e-14. These are diagnostic potentials, not a length-selector study.
The three port gadgets (14),(17),(18) have total dimensions 64,32,80 and
16 code states each. At lambda=0.1,0.05,0.025 all low levels agree with
(21) within 1e-12; eigenpair residuals are below 3.3e-15. Maximum errors
relative to lambda^2 K_2 divided by lambda^4 are

| Port | 0.1 | 0.05 | 0.025 | Limit from (21) |
| :-- | --: | --: | --: | --: |
| Conditional phase | 32.23564 | 34.95887 | 35.73250 | 36 |
| Conditional amplitude | 196.87576 | 237.36309 | 251.00452 | 256 |
| Complementary channels | 14.83519 | 15.68778 | 15.92050 | 16 |

[Rejected checker, repaired] The first Fraction inverse failed its
inverse-identity assertion and emitted fixed-width overflow warnings.
The input conversion Fraction(numpy.int64) retained NumPy integer
arithmetic. The repaired code converts to Python int before Fraction
and uses Python object integers for the identity product. The exact
identity then passes. rejected_checks.json preserves the failure. None
of the claims use the rejected calculation.

[Pipeline scope] The ordered-token scalar extension fails C84's pair
margin identically, before any useful root bracket certificate. The
single-channel conditional phase and amplitude examples fail Gaussianity
at second order. The complementary-channel example passes Gaussianity
but lacks a microscopic encoding and spectrum with n>=a. Therefore
there is no new MEMBER, finite-K budget, all-length selector, or central
charge claim. C84's existing eleven-letter certificate was read, not
rerun as a purported C87 construction.

[Open] The surviving possibilities require a new ingredient such as
conversion-dependent configuration phases, context-dependent exchange,
multiple interacting channels with coordinated resolvents, or a new
marker and charge encoding. None is universally excluded. These results
are obstructions to explicitly stated classes, not an eleven-letter
lower bound for all the proposed extensions.
