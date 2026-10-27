# Boundary Schmidt rank and Clifford compression

Status: proved finite-dimensional statements. These exclude a class of
first-order block constructions, not the full c-d conjecture.

## 1. A single low-Schmidt-rank code vector limits all channels

Let V:C^D -> C^q tensor E be an isometry, P=VV^dagger, and let F_a be
Hermitian operators on C^q. Suppose

    V^dagger (F_a tensor I) V = alpha_a Gamma_a,
    alpha_a real and != 0,  {Gamma_a,Gamma_b}=2 delta_ab I.          (1)

If the code contains a normalized vector psi of Schmidt rank s across
this boundary cut, then

    m <= 2 q s - s^2.                                    (2)

The scales alpha_a may differ and may be arbitrarily small.

Proof. Regard the physical Hilbert space as a real inner-product space
with (u,v)_R=Re<u,v>. Write phi=V^dagger psi. The m vectors

    v_a=P(F_a tensor I)psi/alpha_a=V Gamma_a phi

are orthonormal in this real inner product, because

    Re<Gamma_a phi,Gamma_b phi>
      = <phi|{Gamma_a,Gamma_b}/2|phi>=delta_ab.

They belong to the real image of the linear map

    Herm(q) -> range(P),  F -> P(F tensor I)psi.

Projection cannot increase its real dimension. Before projection, its
kernel consists of Hermitian F annihilating the s-dimensional Schmidt
support S in C^q. In an orthogonal decomposition S direct-sum S^perp,
only the (q-s) by (q-s) bottom block of F may be nonzero. The kernel
therefore has real dimension (q-s)^2. The domain has dimension q^2,
so the image dimension is q^2-(q-s)^2=2qs-s^2. This proves (2).

In particular, a code containing even ONE vector that factors across
the active port satisfies m<=2q-1. For an even number m of Majorana
channels, m<=2q-2, and its first-order target has

    c_target=m/2 <= q-1.                                 (3)

This conclusion does not assume a local parent, frustration freeness,
a gap, or exact invariance of the code under the boundary bond. It
applies to scaled compressions, which are the setting of the proposed
loophole. It is only a bound for the specified Clifford target.

For a proposed m-channel code every nonzero code vector must have

    Schmidt rank >= ceil(q-sqrt(q^2-m))                   (4)

at EACH port that realizes the channels. The integer ceiling is taken
on the exact expression. In particular m>2q-1 requires a completely
entangled code subspace across each boundary cut.

If there is also a Hermitian physical lift F_0 of the logical parity,
anticommuting with every Gamma_a and squaring to I after normalization,
apply (2) to the m+1 generators. For s=1 this still gives m<=2q-2.

## 2. Uniform positive bonds of rank less than q

Let h>=0 act on C^q tensor C^q and have matrix rank r<q. Let G_N be the
COMPLETE kernel of the open-chain sum of h, without boundary penalties.
For every unit vector a at the first site, G_N contains a product state
with that first-site vector. The same statement holds at the last site.

Proof. Choose r vectors w_mu spanning range(h). For a fixed nonzero a,
the conditions

    <w_mu|a tensor b>=0,  mu=1,...,r,

are r homogeneous complex linear equations for q entries of b. There
is a nonzero solution b. Repeat this step to choose the next site. Every
bond of the resulting normalized product state is in ker(h). For a
specified last site, repeat the argument in reverse. QED.

Consequently (3) excludes a first-order Majorana violation from every
such complete ground code. There is a stronger statement.

Suppose a Hermitian F has compressed action A on G_N and A^2=alpha^2 I
with alpha>0. Every eigenvector of F can be extended to a ground product
state, so the maximum and minimum eigenvalues of F occur as Rayleigh
quotients of A. Hence ||F||<=alpha. Conversely ||A||<=||F||, so
||F||=alpha and F^2<=alpha^2 I. For psi in G_N,

    alpha^2 ||psi||^2 = ||P F psi||^2
       <= ||F psi||^2 <= alpha^2 ||psi||^2.

Equality forces (I-P)F psi=0. Thus every such F preserves the code.
Because every boundary vector occurs in a code product state, a boundary
operator vanishing on the whole code is the zero physical operator.
The Clifford identities of A_a therefore imply the same identities for
F_a/alpha_a themselves on C^q. For even m this implies

    2^(m/2) divides q,   m/2 <= nu_2(q) <= log_2 q.        (5)

Here the representation fact used is that the complex even Clifford
algebra is M_(2^(m/2)), whose finite-dimensional unital representations
are direct sums of its defining representation.

Furthermore both endpoint lifts preserve G_N. Since the physical
endpoint operators commute when N>=2, their logical restrictions
commute. Thus this subclass cannot realize even one pair of the desired
joint maps Gamma_a and i Gamma_a P_logical, whose commutator is nonzero.
This conclusion does not assume that the eventual inter-cell Hamiltonian
preserves the code exactly. It is forced already by the low-rank parent
and the exact involution compressions.

These claims concern the entire ground kernel. Selecting a smaller
subcode by new terms can destroy the product-extension property. Extra
boundary penalties, higher-rank bulk constraints, frustrated parents,
and approximate rather than exact Clifford maps are not excluded.

## 3. A quantitative version for approximate Clifford relations

Let A_a=V^dagger F_a V/alpha_a be Hermitian. For any unit code vector phi
of boundary Schmidt rank s, form its real Gram matrix

    M_ab=Re<A_a phi,A_b phi>.

It is positive semidefinite and rank(M)<=R=2qs-s^2. If
||{A_a,A_b}/2-delta_ab I||<=epsilon entrywise, then
|M_ab-delta_ab|<=epsilon and ||M-I||_op<=m epsilon.
For m>R, M has a zero eigenvalue, so

    epsilon >= 1/m.                                      (6)

A sharper, readily tested formulation is ||M-I||_op>=1.
Thus arbitrarily accurate normalized Clifford relations cannot bypass
the rank obstruction at fixed m. This does not give a lower bound on
unnormalized errors independent of the small alpha_a.

For the low-rank-parent setting, if ||A^2-alpha^2 I||<=eta for the
UNNORMALIZED compression A=V^dagger F V, extension of extremal physical
eigenvectors gives ||F||^2<=alpha^2+eta. Since A^2>=alpha^2-eta,

    ||(I-P) F V||^2 = ||V^dagger F^2 V - A^2|| <= 2 eta.   (7)

For disjoint ports F,B, their normalized compressions obey
||[A_F,A_B]||<=2 ||Q F V|| ||Q B V|| after using [F,B]=0.
If their square defects are eta_F, eta_B, this gives
||[A_F,A_B]||<=4 sqrt(eta_F eta_B). This is an unnormalized
inequality; the target commutator has size 2|alpha_F alpha_B|.
Therefore a parity-dressed pair requires

    sqrt(eta_F eta_B) >= |alpha_F alpha_B|/2.              (8)

The required relative error is of order one. The low-rank FF-parent
obstruction is robust against improving the normalized accuracy.

## 4. Two additional local-parent classes

A commuting-projector parent cannot supply the joint noncommuting ports.
Let P be the product of its commuting ground-space projectors. Assume
no individual projector touches both active endpoints. For endpoint
operators F and B, partition the projectors into those touching only
the first endpoint, only the last endpoint, and neither. The last
family commutes with both operators; the first family commutes with B;
the second commutes with F. These commutations and idempotence give

    P F P B P = P F B P = P B P F P.                      (9)

Hence the compressed endpoint operators commute, ruling out the
parity-dressed Clifford pair. The argument does not require F or B to
preserve P, or exact invariance under the later inter-cell perturbation.
For a nearest-neighbor commuting-projector parent, three or more sites
ensure that no parent term touches both endpoints. This statement also
covers commuting positive bonds after replacing them by kernel
projectors. It does not cover a noncommuting frustration-free parent.

A separate dimension obstruction applies to a standard injective MPS
code with a nonzero nearest-neighbor frustration-free parent. State
its assumptions explicitly. The virtual bond dimension is D_v, the
two-site MPS map from M_(D_v) to C^q tensor C^q is injective, its image
is annihilated by a nonzero positive two-site parent term, and the
complete open-chain code for the cell is the D_v^2-dimensional MPS
edge space. Then

    D_v^2 < q^2,  hence D_v<=q-1.                         (10)

The strict inequality follows because otherwise the injective image
would fill the two-site physical space and no nonzero positive term
could annihilate it. An even m=2n Clifford family on the complete edge
space requires 2^n to divide D_v^2, so

    n <= 2 nu_2(D_v) <= q-1.                             (11)

For the last inequality put k=nu_2(D_v). If k=0 it is immediate. If
k>=1, 2k<=2^k<=D_v<=q-1. Thus this MPS-parent subclass cannot provide
c_target=n>q-1, even without inspecting its endpoint maps. The result
uses the complete edge-space dimension; a differently selected subcode
or a parent with larger interaction range is not covered by (10).

These exclusions identify features needed by the unresolved
frustration-free route. A useful complete code must have no product
vector at either active boundary, use noncommuting constraints, avoid
the low-rank product-extension setting, and go beyond the specified
injective two-site MPS parent class. None of these requirements proves
that such a code cannot exist.


## 5. Consistency with the earlier joint Choi code

The earlier unrestricted dilation is not excluded by the new theorem.
For a unit logical state phi, its left-port reduced state is

    rho_L=(I+alpha sum_a <Gamma_a>_phi F_a)/q.

The Clifford expectation vector obeys sum_a <Gamma_a>^2<=1. Indeed,
for real u, (sum_a u_a Gamma_a)^2=||u||^2 I; taking u equal to the
expectation vector proves the claim. Since tr(F_a F_b)=q delta_ab,

    ||sum_a <Gamma_a> F_a|| <= ||sum_a <Gamma_a> F_a||_HS
                             <=sqrt(q).

Therefore every eigenvalue of rho_L lies between
(1-alpha sqrt(q))/q and (1+alpha sqrt(q))/q. The same proof applies
to the right-port Clifford family. At the earlier
alpha=1/(4m sqrt(q)), every nonzero code vector has full Schmidt rank
q at BOTH ports, with entropy at least log(q)-log(1+1/(4m)).
This is an analytic check on the scope of (2), not a new numerical
claim or a construction of a local parent. For the q=9,m=80 resource
example, (4) requires rank at least 8 at each port for every code
vector; the Choi code has rank 9.
