# Common-rotation invariants and an honest simplicity count

Use the canonical bond with field `h(ZI+IZ)/2`. Its exchange matrix is

```
J=[[a,b,0],[-b,-a,0],[0,0,c]],
S=(J+J^T)/2=diag(a,-a,c),
A=(J-J^T)/2,   D=b e_z,   B=h e_z.
```

A common proper rotation sends `(S,D,B)` to `(R S R^T,R D,R B)`.
The following invariants therefore hold before any infrared analysis:

```
tr S=c,
[tr(S²)-(tr S)²]/2=a²,
||A||_F²/2=b²,
||B||²=h².                                               (1)
```

Their Jacobian with respect to `(c,a,b,h)` has rank four whenever
`abh!=0`. After a positive scale is removed, there remain three
continuous ratios. A discrete sign choice cannot reduce that dimension.
The signs of a can be reversed by a common pi/2 rotation about z;
a common pi rotation about x reverses b and h together. In particular
`D.B=bh` is another invariant. No such rotation can make any nonzero
coefficient vanish while retaining a member of the stated family.
The signed trace c, not just its magnitude, is fixed before rescaling.

The singular values of J are

```
{|a+b|,|a-b|,|c|}.                                      (2)
```

Singular values alone would hide a small antisymmetric part. Keeping
`||S||_F`, `||A||_F`, `|B|`, and the relative orientations prevents that.
Adding a constant or a telescoping bond term changes none of J or the
total uniform B. Allowing those equivalences does not evade (1).

## Four monomials, still four independent structures

Let `U=exp(-i pi Z/8)`. Direct conjugation gives

```
U X U†=(X+Y)/sqrt(2), U Y U†=(-X+Y)/sqrt(2),
(U tensor U) h_bond (U tensor U)†
  =(a+b)XY+(a-b)YX+cZZ+h(ZI+IZ)/2.                       (3)
```

Thus a generic member has a four-monomial representation. The new
coefficients recover the old ones by half-sum and half-difference, so
this is no loss of an adjustable structure. A special relation `a=b`
does eliminate YX; together with `h=0` it gives the genuine two-term
subfamily `2a XY+c ZZ`. That relation must be imposed physically; a
common rotation does not impose it on the generic reference.

For the earlier reference from
[the comparable-coupling search](../d2_order_one_20261110/README.md),

```
a=1/2, b=1/160, c=10, h≈-19.9760472506,
rotated transverse coefficients = 81/160,79/160.
```

The original literal coefficient spread is about 3196.17 and the
rotated spread is about 40.46. But the basis-independent ratios
`|a/b|=80`, `|c/b|=1600`, and
`||A||_F/||S||_F = sqrt(2)/[160 sqrt(100.5)]≈0.000882`
are unchanged. The last decimal field is quoted as a diagnostic from
the existing numerical reference, not as a newly specified exact bond.

## Why the zero-field face is not secretly XXZ or free by onsite rotations

For an additive charge `Q=sum n_j.sigma_j`, bondwise Pauli independence
requires `Omega(n_j)J=J Omega(n_(j+1))`, where Omega is the cross-product
matrix. Hence a nonzero n_j must commute with both `JJ^T` and `J^TJ`.
If all singular values are distinct this is impossible. When `ab!=0`,
the two planar singular values in (2) are different. A collision of
`c²` with one of them leaves a unique distinct left singular axis and
a unique distinct right singular axis that are orthogonal planar
directions. They cannot both be the axis of the same n_j. Three equal
singular values require `ab=0`. At `a=±b,c=0`, the rank-one matrix has
noncollinear left and right axes and also has no additive charge.
These cases exhaust the possibilities: **if ab!=0, no nonzero additive
onsite charge exists**, even with h=0 and arbitrary site-dependent axes.
This specializes the [earlier exhaustive charge theorem](../d2_gap_reduction/additive_charge.md).

For a standard onsite-rotated Jordan–Wigner quadratic representation,
the local parity axes must obey both `J n=0` and `J^T n=0` at every
site of a uniform chain. If `abc!=0` and `a²!=b²`, J is invertible.
When `a=±b!=0,c!=0`, its left and right null axes are orthogonal, so
there is still no common null axis. The simple residual `2a XY+cZZ`
therefore does not become free by a different onsite Jordan–Wigner
axis. This does not exclude a nonlocal duality or unknown integrability.
