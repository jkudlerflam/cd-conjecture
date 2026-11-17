# Exact reductions and an obstruction to the usual third charge

All spins are Pauli matrices, the lattice spacing is one, g is real, and
`H(g)=sum_j h_j`, `h_j=2 X_j Y_(j+1)+g Z_j Z_(j+1)`.
Statements about periodic local charges use L>=6; bulk identities also
hold on the infinite chain. There is one coupling ratio after fixing a
common positive energy scale. The literal monomial spreads at g=1/2,1,2
are respectively 4,2,1. No constants or telescoping terms enter these counts.

## [Proved] Useful maps and symmetry sectors

A common rotation Rx(pi/2) sends Y to Z and Z to -Y, giving

```
H(g) -> 2 sum XZ + g sum YY = 2 H_C144(alpha=0,beta=g/2).
```

Thus the earlier C144 data at alpha=0,beta=1/2 are exactly g=1,
with energies multiplied by two and entropies unchanged. They are not
data on the tuned slow-hopping cancellation curve.

The four-site product of onsite Pauli operators `(I,Y,Z,X)`, repeated,
keeps every XY bond and reverses every ZZ bond. It maps g to -g exactly
on open chains and periodic rings of length divisible by four. For
L=2 mod 4 the closing XY bond acquires a minus sign; neglecting that
twist would give a false finite-size spectral identity. The bulk phase
diagram is even in g; staggered Z order is exchanged with uniform Z order.

Rotating only odd sites by Rz(pi/2) produces the real Hamiltonian

```
H' = -2 sum_(j even) XX + 2 sum_(j odd) YY + g sum_j ZZ.       (1)
```

For even periodic L, Pi=product Z and F=product X commute with (1).
Their four blocks each have dimension 2^(L-2). Physical one-site
translation satisfies `T F T^-1 = Pi F`; hence it pairs the F=+/-
blocks at Pi=-1. This exact pairing of excitation sectors does not
establish spontaneous symmetry breaking. Pi=+1 blocks retain physical
one-site momentum; Pi=-1 blocks can use two-site momentum, with the
other block restoring the k versus k+pi pair. The numerical code retains
all four physical sectors, solving the duplicate only in validation.

The all-sector physical matrices at L=4,6,8, the translation intertwiners,
and the sign duality are checked independently in `check.py`.

## [Proved] Majorana representation, with the boundary qualification

Let `alpha_j=(product_(l<j) Z_l) X_j`,
`beta_j=(product_(l<j) Z_l) Y_j`, and define
`e_A,j=i beta_j beta_(j+1)`, `e_B,j=i alpha_j alpha_(j+1)`.
The interior-bond identities are

```
XY=-e_A,  YX=e_B,  ZZ=-e_A e_B,
H=-2 sum e_A - g sum e_A e_B.                              (2)
```

Open chains satisfy these identities without a seam. Periodic chains
require the usual Pi-dependent Jordan-Wigner boundary twist, retained
by the spin implementation. For g!=0, (2) is an interacting ladder.
Counting its two lattice species is not a bound on interacting Virasoro c.

The earlier onsite-Jordan-Wigner obstruction also specializes here:
the bond coupling matrix has rank two, with left null axis Y and right
null axis X. Site j would have to use both axes as its local fermion-parity
axis to eliminate all longitudinal interactions on its two incident
bonds. They are orthogonal, so no choice of onsite axes makes the whole
chain a nearest-neighbour quadratic Jordan-Wigner chain for g!=0.

## [Proved] No faithful term-by-term free-Majorana map

The anticommutation graph contains a claw with center `(XY)_0` and
leaves `(XY)_(-1),(XY)_1,(ZZ)_(-1)`: the center anticommutes with each
leaf, and the three leaves commute pairwise. Distinct Majorana bilinears
anticommute precisely when their two endpoints share exactly one
Majorana. An edge has only two endpoints, so three pairwise disjoint
incident edges cannot exist. This is a contradiction.

This rules out faithful maps taking every nonzero Pauli summand to one
Majorana bilinear, including faithful Clifford/Jordan-Wigner/Kramers-Wannier
term maps with ordinary boundary parity bookkeeping. It does not rule
out maps to sums, restricted gauge representations, arbitrary nonlocal
diagonalization, or interacting integrability. This is the specific
claw argument already identified in the repository's
[direct-attack alternatives](../direct_attack_20261028/alternatives/RESULT.md).

## [Proved] The boost third charge cannot be repaired at range two

Direct Pauli multiplication gives

```
[h12+h23,[h12,h23]]
 = (32+8g²)(IXY-XYI) + 16g(IZZ-ZZI)
   + 32g(YIY-XIX).                                        (3)
```

The first two terms telescope. The last term has nonidentity operators
at sites one and three, so it cannot equal `I tensor K-K tensor I`.
Consequently the canonical bond fails the regular difference-form
Reshetikhin condition whenever g!=0. See the derivation in
[Mansson, section 2, Eqs. (2.3)--(2.7)](https://arxiv.org/pdf/hep-th/0703150).

There is a stronger, density-gauge-independent statement. Set
`Q3=i sum_j [h_j,h_(j+1)]`. Jacobi cancellation of the four-site terms
and telescoping of the first two terms in (3) give

```
[H,Q3] = 32 i g sum_j [Y_j Y_(j+2)-X_j X_(j+2)].            (4)
```

For any operator R that is a sum of operators supported within adjacent
pairs, `[H,R]` has **no** Pauli word `XIX` or `YIY`. Indeed:

* disjoint supports commute;
* equal or contained supports give range at most two;
* on pairs overlapping at one site, a nonzero Pauli commutator requires
  different nonidentity Pauli matrices on the shared site, and their
  product is a nonidentity Pauli matrix there.

Thus no R of range at most two cancels (4) for g!=0. In particular,
`Q3+R` is never conserved. Replacing the density by
`h_j+u_j-u_(j+1)` changes Q3 only by such an R; a telescoping density
gauge therefore cannot evade the obstruction.

This excludes the standard homogeneous regular difference-form
Yang-Baxter/boost construction, and any construction requiring this
specific leading third charge with a range-two correction. It is not a
classification of all possible conserved quantities or all definitions
of integrability; staggered larger-cell constructions and other charges
have not been exhausted. It implies no bound on c.

`check.py` verifies (3) using exact integer Pauli algebra at four g
values, which certifies the degree-at-most-three polynomial identity.
It checks the repair obstruction on a complete two-site Pauli basis
and (4) separately with literal periodic spin matrices.
