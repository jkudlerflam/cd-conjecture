# The sharp frustration-free witness length for qubits

Status: proved. This concerns a positive uniform bond on an OPEN chain
with no boundary terms. It does not classify conformal frustrated chains.

The finite-witness integer from `finite_ff_witness.md` has the sharp value

    K_FF(2)=4.                                           (1)

Thus a positive two-qubit bond that is frustration free on four open
sites is frustration free at every length. A finite longest FF segment
can have two or three sites, but cannot be longer than three.

## 1. Reduce to a two-dimensional allowed space

Put S=ker(h) subset C^2 tensor C^2. The zero space of H_N is

    G_N = intersection_j C^2 tensor ... tensor S_(j,j+1)
                              tensor ... tensor C^2.     (2)

Only S matters, not the nonzero eigenvalues of h.

If dim(S)=0, there is no zero vector at N=2. If dim(S)>=3,
rank(h)<=1<2, so the product-extension proof in
`boundary_rank_obstruction.md` supplies a zero vector at every length.
Alternatively a single homogeneous quadratic equation <w|a,a>=0 has a
nonzero solution over C, giving the repeated product a^tensor N.

If dim(S)=1, write S=span(v). An extension to three sites requires
v_12 tensor x_3=y_1 tensor v_23 for nonzero vectors x,y. Taking Schmidt
ranks across site 1 versus sites 2,3 forces v to have rank one. Thus
v=a tensor b. The same equation then forces a proportional to b. This
case has the repeated product a^tensor N at all lengths. Otherwise
G_3=0.

It remains to examine dim(S)=2.

## 2. The product vectors in a two-dimensional two-qubit space

Identify a vector with its 2 by 2 coefficient matrix. Product vectors
are the zeros of the determinant. On a projective line P(S), determinant
is a homogeneous quadratic polynomial. There are three cases over C.

(a) It vanishes identically. Then S is a line entirely of product
vectors, necessarily a tensor C^2 or C^2 tensor b. To verify the last
statement, take two independent products a tensor b and c tensor d.
If both a,c and b,d were independent, their sum would have matrix rank
two. Hence one pair is dependent. The subspace contains a tensor a or
b tensor b and is FF at all lengths.

(b) There are two distinct product directions. Write

    S=span(a_0 tensor b_0, a_1 tensor b_1).               (3)

Both a_0,a_1 and b_0,b_1 are bases, since otherwise case (a) applies.
For a three-site vector in S_12 tensor C^2, expand the third site in
the b basis:

    psi=sum_(i,j) c_ij a_i tensor b_i tensor b_j.

Membership in C^2 tensor S_23 says that the slice with last vector b_j
has middle factor a_j. Independence of the first-site a_i then gives

    c_ij !=0 only if b_i is proportional to a_j.          (4)

Define a directed graph on the two bond labels, with edge i->j if
b_i is proportional to a_j. Each vertex has at most one incoming and
one outgoing edge. Equation (4), iterated along the chain, says that
G_N is spanned by the product vectors associated to graph paths of
N-1 bond labels. The iteration is valid because distinct paths have
distinct first labels and hence independent first-site factors.

If a directed cycle exists, following it gives a product ground state
at every open length. If no cycle exists, a path visits at most the two
vertices. Therefore G_N=0 for N>=4. In particular G_4 nonzero forces
FF at every length.

(c) There is one double product direction. The line is tangent to the
rank-one quadric. It can be written

    S=span(a tensor b, a tensor d+c tensor b),            (5)

where a,c and b,d are bases. To see the normal form, transform the
unique product vector to |00>. For a second matrix M, the determinant
of x|00>+y M has a double root at y=0 only if M_11=0. Subtract its
00 entry, leaving the nonzero 01 and 10 entries in (5). Neither can
vanish, because the line is not all-product.

If b is proportional to a, the repeated product a^tensor N works at
all lengths. Otherwise any element of S_12 tensor C^2 is

    psi=a tensor b tensor x
         +(a tensor d+c tensor b) tensor y.

Its slice with first vector c is b tensor y. If psi also satisfies
the 23 constraint, this product belongs to S. The only nonzero product
direction in S is a tensor b, and b is not proportional to a, so y=0.
The first-site a slice now forces b tensor x in S, giving x=0. Thus
G_3=0.

Every case proves the upper bound K_FF(2)<=4.

## 3. Sharpness and a general lower bound

Let |+>=(|0>+|1>)/sqrt(2), and use the explicit rank-two projector

    h=I-|0,1><0,1|-|1,+><1,+|.                           (6)

The two subtracted vectors are orthonormal. Its allowed-space graph
is the single path from |0,1> to |1,+>. The three-site ground space is
span(|0,1,+>), but G_4=0. Hence K_FF(2)>3, proving (1).

The same idea works at any q>=2. Choose v=q^(-1/2) sum_j |j>, and let
S be the span of the q orthonormal products

    |0,1>, |1,2>, ..., |q-2,q-1>, |q-1,v>.

Put h=I-P_S. The left factors are a basis, and the right factors
|1>,...,|q-1>,v are another basis. The graph proof of case (b) applies
unchanged. Its only maximal path has q bond labels, so

    G_(q+1)=span(|0,1,...,q-1,v>),  G_(q+2)=0.

Consequently the general finite-witness length obeys

    K_FF(q)>=q+2, q>=2.                                  (7)

The nonorthogonal final local letter matters. A directed graph built
only from an orthonormal physical alphabet has a smaller acyclic path
bound. There is no contradiction, since S here is a quantum subspace.
No matching upper bound for q>=3 is proved.

## 4. Scope for the block construction

This makes the earlier noetherian obstruction explicit for qubits. A
longest-FF-segment selector using only a uniform positive qubit bond
can select at most three occupied sites. Boundary penalties and
frustrated kinetic selectors remain outside the theorem. Even within
this subclass, (1) is a statement about exact zero modes, not about
central charge or conformal onset lengths of arbitrary chains.
