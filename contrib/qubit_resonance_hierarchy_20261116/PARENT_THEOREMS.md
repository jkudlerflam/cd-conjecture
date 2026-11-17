# Exact resonances, their ground spaces, and a quadratic gap upper bound

These statements concern a fixed finite Hermitian nearest-neighbour qubit
bond, repeated with one-site translation symmetry. The parameter r is
fixed before taking the thermodynamic limit. Throughout q>0 and r>0;
the endpoint r=0 is the singular blockade. The central charge under
discussion would be the total usual nonchiral Virasoro charge of the
complete unitary conformal bulk theory. No such c is assigned to these
parents.

## 1. An exact polynomial microscopic curve [proved]

Use `b†=|1><0|`, `n=b†b`, `D1=sum nn`, `N=sum n`,
`T1=sum b†_(j+1)b_j`, `A=sum b†_j b†_(j+1)` and
`F=sum(n_(j-1)-n_(j+1))(b†_j+b_j)`. Let

```
psi(r,q)=q r²|00>-r|01>+r|10>+|11>,
v=1-2r²+q²r⁴,  g=r+q r³,
t=-r²,         p=q r²,       u=2r²-2q²r⁴.                 (1)
```

The physical periodic Hamiltonian is

```
H=v D1+g F+t(T1+T1†)+p(A+A†)+u N.                         (2)
```

Split the onsite term equally between adjacent bonds. With
`G=q r³ X`, its literal bond satisfies

```
h+q²r⁴ I+G tensor I-I tensor G=|psi><psi|.                  (3)
```

This is an elementary entry-by-entry identity. Therefore the periodic
chain, shifted by `Lq²r⁴`, is a positive rank-one parent. All terms in
(1)--(3) are nearest-neighbour or onsite and independent of L. On an
open interval the G difference is an actual endpoint term; the checks
there use the literal projector parent, retaining those terms.

The positive-normalized v=1 chart is obtained by dividing every
coefficient by v, for sufficiently small r. In the earlier coordinate
`epsilon=g/v`,

```
epsilon=r+(q+2)r³+O(r⁵),
p/v=q epsilon²-2q(q+1)epsilon⁴+O(epsilon⁶).                (4)
```

Thus simply holding `p=q epsilon²` misses the exact resonance at the next
order. Formula (1), rather than its leading Taylor terms, is the exact
curve used below.

## 2. Root-of-unity condition and its completeness in the parent family

For (1) the transfer matrix is

```
T=[[-r,1],[-q r²,-r]],
lambda_+=r(-1+i sqrt(q)), lambda_-=r(-1-i sqrt(q)).         (5)
```

Write theta=arctan(sqrt(q)). In this convention
`lambda_+/lambda_-=exp(-2i theta)`. Therefore

```
q_(m,k)=tan²(pi k/m),  gcd(k,m)=1, 0<k<m/2,               (6)
```

has exact projective period m for every nonzero r. The requested branch
has k=1, so q3=3, q4=1, q5=5-2sqrt(5), q6=1/3 and q_m tends to zero.
Interchanging eigenvalues complex-conjugates the ratio and does not
change its order.

This also classifies the nontrivial root-of-unity loci in the earlier
parent `psi=(c,-a,a,b)`, b real positive and a nonzero. Its transfer
eigenvalues are `-a* +/- sqrt(-bc*)`. A ratio exp(-2i theta), different
from -1, is equivalent to

```
bc*/(a*)²=tan²(theta)>0.                                  (7)
```

A common onsite z rotation makes a real, after which c is real positive.
Positive normalization sets b=1, giving (1). A ratio -1 would require
a=0 and lies outside these nonzero-facilitation arcs. No independent
site-dependent change of microscopic translation symmetry is used.

## 3. The exact transfer recurrence [proved]

Define real polynomials

```
s0=0, s1=1, s_(d+2)=-2 s_(d+1)-(1+q)s_d.                 (8)
```

Cayley--Hamilton gives

```
(T^d)_(01)=r^(d-1) s_d(q),
s_d(q)=Im[(-1+i sqrt(q))^d]/sqrt(q).                      (9)
```

For the primitive k=1 resonance,

```
s_d(q_m)=(-1)^(d+1) sin(d pi/m)
          /[sin(pi/m) cos(pi/m)^(d-1)].                   (10)
```

In particular `T^d` is scalar exactly when m divides d. The first
polynomials are

```
s2=-2, s3=3-q, s4=4(q-1),
s5=q²-10q+5, s6=-2(3q-1)(q-3),
s7=7-35q+21q²-q³, s8=8(q-1)(q²-6q+1).                  (11)
```

This explains the roots that first appear in the distance-three/four
penalties and locates the new distance-five/six resonances. It is an
all-order **transfer** identity. It is not, by itself, a formula for all
effective energy coefficients or for the order at which hopping occurs.

There is a precise boundary-data interpretation. Put
`S_r=diag(1,r)` and `M=[[-1,1],[-q,-1]]`. Then

```
T=r S_r M S_r^-1.                                        (12)
```

The two covectors `e0^T M^j` and `e0^T M^(j+d)` have the same
projective zero exactly when s_d=0. Their determinant is a nonzero
factor times `(M^d)_(01)`. Thus two sites can share a pole of the same
product-ground-state parameter exactly at these separations. This is
the recurrence shortcut used here; it avoids guessing new microscopic
counterterms from a few numerical energy coefficients.

## 4. Exact compatible-ring ground space [proved using the parent theorem]

Bravyi--Gosset Propositions 2--3 give the open-chain ground space as
`(tensor_j T^j) Sym^L(C²)` and its periodic dimension as L+1 when T^L is
scalar, two otherwise. These statements apply because `det T=r²(1+q)`
is nonzero. Their Section 2 supplies the complete ground spaces, not
only product-state examples.

Consequently, for (6):

* open intervals have exactly L+1 ground states;
* periodic rings with m dividing L have exactly L+1;
* incompatible rings have exactly two, spanned by the two uniform
  eigenvector product states of T.

For compatible L=mM, (12) rewrites the ground space as

```
G_L(r)=S_r^(tensor L) [tensor_(j=0)^(L-1) M^j] Sym^L(C²). (13)
```

The common `S_r^-1` factors preserve the symmetric subspace and have been
absorbed. Since M^m is scalar, one may use j modulo m in (13). All
vectors in G_L(r) are invariant under translation by m sites, and their
translation momenta lie in the m-point grid `2pi Z/m`.

Equation (13) describes one transformed symmetric multiplet of dimension
L+1, not a tensor product of m independent multiplets. The transformation
is nonunitary and becomes singular at r=0. It preserves exact nullspaces
and Schmidt ranks; it does not establish a uniform spectral equivalence
to the ferromagnetic Heisenberg chain.

Every ground state has Schmidt rank at most ell+1 across an ell-site
interval cut, because a symmetric state does and all the maps in (13)
are products of invertible onsite maps. This bound is independent of m.
It is a statement about the exact ground space, not nearby Hamiltonians.

Source: [Bravyi--Gosset, Sections 2.1--2.2](https://arxiv.org/html/1503.04035).

## 5. Its singular ground-space limit is a sum, not a tensor product

**Theorem [proved].** Fix m and a compatible L=mM. As r tends to zero,
the orthogonal projector onto (13) tends to the projector onto

```
G_L(0+) = span {vacuum;
                 Dicke_(a,N): a=0,...,m-1, N=1,...,M}.     (14)
```

Here `Dicke_(a,N)` is the uniform N-particle state supported entirely
on sites congruent to a modulo m. These vectors are orthogonal and
their number is `1+mM=L+1`.

Proof. The symmetric space is spanned by `v(z)^(tensor L)`, with v(z)
a projective one-qubit parameter. The resulting product ground states
have onsite vectors `S_r M^j v(z)`. Among the m distinct covectors
`e0^T M^j`, choose the unique zero z_a belonging to residue a. Their
zeros are distinct by (9)--(10). Set `z=z_a+r w` in a regular local
projective chart. On residue a, both onsite components are of order r
and their normalized ratio can be chosen arbitrarily by w. On all
other residues the first component is nonzero and the second tends to
zero. After normalization the product state therefore tends to any
coherent symmetric state on residue a, with vacuum on the others.
Varying w spans its M+1 dimensional symmetric space. Their union spans
(14), of dimension L+1. Since G_L(r) has this same fixed dimension,
these independent limiting vectors exhaust the limit. This also proves
convergence of its orthogonal projector. No thermodynamic limit is
exchanged with r tending to zero.

Thus the enlarged ground space supplies m possible residue labels for
one collective ordered sector. It does not supply an extensive tensor
product of independently gapless label degrees of freedom. This does
not exclude extra label dynamics among excited states after detuning.

## 6. A nonconformal variational obstruction at the exact parent

**Theorem [proved].** On a compatible ring L=mM with M>=2, the smallest
positive eigenvalue of the unnormalized projector parent obeys

```
gap_L <= 4r²(1+q_m) sin²(pi/L).                           (15)
```

Proof. Put `rho=r sqrt(q)`, and take orthonormal onsite vectors

```
phi=(1,i rho)/sqrt(1+rho²), w=(i rho,1)/sqrt(1+rho²).
```

The uniform product phi is an exact ground state. On one bond,

```
<psi|w phi>=r(1+i sqrt(q))=x,
<psi|phi w>=r(-1+i sqrt(q))=y.
```

The normalized one-flip trial state
`L^-1/2 sum_j exp(ikj)|w_j phi_elsewhere>` has Rayleigh quotient

```
|x+exp(ik)y|²=r²(1+q)|exp(ik)-exp(2pi i/m)|².             (16)
```

Choose `k=2pi/m+2pi/L`. It is a valid ring momentum but is outside the
ground-space momentum grid in Section 4. The state is therefore
orthogonal to the *entire* ground space. The variational principle gives
(15). For the normalized local projector, divide the bound by
`||psi||²=1+2r²+q²r⁴`.

The trial subspace need not be invariant. Thus (15) is a quadratic upper
bound, not a proof that the exact dynamical exponent equals two. It is
enough to exclude a complete compact unitary conformal finite-size
spectrum with strictly positive velocities at this exact parent. The
growing L+1 ground degeneracy is an additional warning against assigning
it an ordinary positive-c conformal interpretation.

Bravyi--Gosset's Theorem 1 separately gives an O(1/L) upper bound for
the normalized **open** chain here, since the transfer eigenvalues have
equal nonzero moduli. It does not identify the universality class or
control detunings. The stronger compatible-ring bound (15) is derived
above rather than attributed to that theorem.
