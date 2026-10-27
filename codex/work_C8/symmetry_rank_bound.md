# A bond-Schmidt bound for microscopic current sectors

Status: proved finite-dimensional symmetry statement. The central-charge
corollary assumes the stated CFT current correspondence and Sugawara
completeness. Neutral or emergent sectors remain uncontrolled.

## 1. Canonical connected bond

Remove the scalar and one-site terms from a Hermitian bond and write

    H=sum_i b_i + sum_i v_(i,i+1),
    tr_left(v)=tr_right(v)=0,
    v=sum_(a=1)^r s_a A_a tensor B_a,  s_a>0.              (1)

The A_a and B_a are orthonormal traceless Hermitian operators for the
Hilbert-Schmidt inner product. This is a REAL singular value decomposition
of the Hermitian coefficient matrix. The positive s_a are the connected
operator Schmidt singular values; r is the connected Schmidt rank r(h).
Let n_1,...,n_t be the multiplicities of its distinct singular values.

Suppose a compact connected on-site group G, acting by U(g) on C^d,
preserves b and v under simultaneous conjugation. This holds
automatically if the uniform on-site action preserves the full ring
Hamiltonian for every length at least three. After centering as in (1),
one-site and distinct two-site traceless supports are orthogonal in
Hilbert-Schmidt inner product, so their changes under conjugation
cannot cancel. Let g be its Lie algebra after removing scalar generators.

There is an ideal k of g such that every generator in k commutes with
b and with every A_a and B_a individually. Moreover

    g/k embeds in direct_sum_j so(n_j).                  (2)

In particular

    dim(g/k) <= sum_j n_j(n_j-1)/2 <= r(r-1)/2.            (3)

These are statements about the interaction channels, not about the full
on-site symmetry dimension. The kernel can be large.

## 2. Proof

The real left Schmidt space span_R{A_a} is intrinsic. It is the image
of partial contractions of v on the right. Since simultaneous on-site
conjugation preserves v, conjugation preserves this space, and likewise
the right space. Unitary conjugation preserves the Hilbert-Schmidt
inner product, so its matrices O_L(g),O_R(g) are orthogonal. They lie
in SO(r) for connected G. Invariance of v is

    O_L S O_R^T=S,   S=diag(s_1,...,s_r).                 (4)

Multiplying (4) by its transpose gives

    O_L S^2 O_L^T=S^2.

Thus O_L preserves each eigenspace of S^2, and its infinitesimal action
lies in the direct sum of so(n_j). Since S is positive, O_L commutes
with S itself. Equation (4) then gives O_R=O_L in the paired Schmidt
bases. In particular their infinitesimal kernels coincide.

Let k be this common kernel. A generator in k commutes with every
Schmidt operator on each side. It also commutes with b by hypothesis.
At any site i its local generator Q_i therefore satisfies

    [H,Q_i]=0                                            (5)

on every ring. The kernel is an ideal, and the quotient representation
is faithful by its definition. This proves (2) and (3).

The use of singular-value multiplicities is material. A connected bond
with pairwise distinct singular values has no nontrivial continuous
action on its interaction space, even if there are local conserved
on-site observables. Degeneracies permit the continuous rotations.

## 3. Conditional central-charge and torus bounds

Assume the following infrared input.

(i) The relevant microscopic charge densities have scaling limits
whose nonzero dynamical components supply the independent affine
currents of a unitary CFT, with nonzero finite velocities.

(ii) These currents exhaust the stress tensor, possibly as independent
factors with different velocities. There is no positive-charge coset.

A generator in k cannot supply a nonzero dynamical affine current under
(i). Equation (5) says its density is stationary at each individual
lattice site, so its dynamical spectral measure is supported at omega=0
at every spatial momentum. A finite-velocity affine current has
nonzero-momentum spectral weight at omega=v|k|. This weight cannot be
the scaling limit of that identically stationary density. Removing
these inactive generators, the joint current algebra is an image of
g/k and also an image of the microscopic algebra in pu(d).

For a compact reductive current algebra l, the supplied Sugawara input is

    c_Sug=rank(abelian part)
           +sum_simple k_a dim(l_a)/(k_a+h_a^vee)
          <=dim(l).                                     (6)

The inequality is strict on any nonzero simple factor. Assumption (ii)
sets c=c_Sug. Therefore

    c <= min[d^2-1, sum_j n_j(n_j-1)/2].                  (7)

No affine-level bound is used. For r<=3 this proves c<=r within this
current-complete class. For r=1 it gives c=0, so any critical rank-one
chain with positive c must have neutral or emergent sectors outside
these hypotheses. This is consistent with the conditional rank-one
construction in RESULT.md Section 6, whose propagating U(1)s are
emergent rather than exact microscopic charge-transfer symmetries.

There is a stronger conclusion if the infrared is Gaussian and all its
bosons are resolved by commuting microscopic on-site charges. The
abelian image in (2) has dimension at most the rank of its compact
ambient Lie algebra. Commuting real skew-symmetric matrices can be
simultaneously put into independent 2-plane rotation blocks, so

    c <= min[d-1, sum_j floor(n_j/2)] <= floor(r/2).       (8)

The identification c=number of independent resolved currents is the
same Gaussian-current completeness hypothesis used in C5/C7. Equation
(8) controls these modes by the actual bond rank, and holds regardless
of perturbative order or cell size when applied at the microscopic
site level. A global on-site symmetry alone is insufficient, because
it need not resolve every infrared sector.

## 4. Examples and limits

For an XX qubit bond, its two equal singular values give n_1=2 and
c<=1 under current completeness. For XXZ with three nonzero channels
and unequal transverse and longitudinal strengths in magnitude, the
multiplicities are 2 and 1, again giving c<=1. This includes arbitrary
levels in principle; the active algebra is abelian.

For an isotropic Heisenberg bond, all three singular values coincide,
so (7) gives c<3 for an SU(2) current-complete theory. It does not derive
the sharper spin-1/2 bound c<=1 or an affine-level ceiling.

For two decoupled XX chains in one d=4 site with unequal strengths,
r=4 and the multiplicities are 2,2. Both (7) and (8) give c<=2, matching
the two charge modes. At equal strengths the dimension bound becomes
6, while the torus bound remains 2.

An Ising critical chain is not a contradiction at rank one. Its c=1/2
sector has no microscopic continuous affine currents. A neutral coset
or emergent currents remain the missing part of any general rank bound.
No conclusion about those sectors is drawn here.

`schmidt_symmetry_checks.py` solves the full on-site infinitesimal
symmetry equations for explicit bonds, computes their induced action
on the actual Schmidt space, and checks dimensions and singular-value
multiplicities. These are algebraic checks, not measurements of c.
