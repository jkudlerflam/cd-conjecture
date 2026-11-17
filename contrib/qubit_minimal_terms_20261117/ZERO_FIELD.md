# A conditional local bound on the field-free three-coupling face

## Result and hypotheses

**[Derived at the perturbative continuum standard, not a global theorem.]**
Fix b!=0 and a compact subinterval of

```
-1/sqrt(2) < Delta=c/|b| < 1,   Delta!=0.                (1)
```

Starting from a=0, a sufficiently small nonzero |a/b| preserves the
single compact-boson phase with total C=1. The hypotheses are:

1. The exact XXZ low-energy matching at a=0, including its operator
   normalization and compactification, is the reference continuum theory.
2. Perturbative local renormalization about that fixed line controls
   sufficiently small a; generated local operators respect the exact
   microscopic symmetries and the velocity remains positive.
3. The statement concerns this neighbourhood. It does not assume that
   an arbitrary strongly interacting point on the same face has the
   same set of low-energy fields.

No numerical value of the allowed |a/b| is obtained. In particular an
order-one rational point is not certified by calling it "nearby".
The exclusion of Delta=0 is for the strict irrelevance argument; c=0
is classified exactly in REDUCTIONS. Compact subintervals avoid
vanishing irrelevance margins and the critical endpoints.

## Derivation: the operator lattice, not just the number of species

Use the spiral from REDUCTIONS with b>0. The zero-field chain becomes

```
sum_j [b(XX+YY)+cZZ+(-1)^j a(XX-YY)].                    (2)
```

The physical one-site translation is now translation followed by a
quarter-turn about z. It is an exact symmetry of (2); ordinary
one-site translation alone is not. The XXZ reference has

```
K=pi/[2(pi-arccos Delta)],
H_LL=(v/2pi) integral [K (d theta)²+K^-1 (d phi)²].       (3)
```

Use conventions in which `exp(2i phi)` is the staggered density vertex
and `exp(2i theta)` is a pair vertex. The transformed physical translation
shifts both phi and theta by pi/2, modulo their compact periods. Fermion
parity forbids single-spin creation vertices. A local parity-even
vertex can be written

```
V_(p,q)=exp[2i(p phi+q theta)],  p,q integers,
dimension x_(p,q)=p² K+q²/K,
conformal spin s_(p,q)=2pq.                             (4)
```

The translation selection rule is **p+q even**. In particular, the
dimension-K or dimension-1/K mass vertices (1,0) and (0,1) are forbidden.
The leading non-Gaussian allowed operators are precisely

| Charges | Dimension | Role |
|---|---:|---|
| `(1,1),(1,-1)` and conjugates | `K+1/K` | Staggered pairing / velocity anisotropy; conformal spin +/-2 |
| `(2,0)` and conjugate | `4K` | Density umklapp, already compatible with XXZ |
| `(0,2)` and conjugate | `4/K` | Pair umklapp, allowed once a is nonzero |

The last two are also produced by products of the first operators;
ignoring generated operators would miss the boundaries of the argument.
Every other nonzero allowed integer pair has dimension at least one
of these three values. Descendants only increase it. The remaining
dimension-two Gaussian derivatives shift K and v. Symmetries at h=0
exclude a density bias; a total derivative is not a new bulk mode.

For `1/2<K<2` and K!=1, all three non-Gaussian dimensions are strictly
above two. Thus a small perturbation and its generated terms are
irrelevant within the stated continuum regime, leaving C=1. Equation
(3) converts `1/2<K<2` exactly into (1). At K=1 the first perturbation
is marginal: the exact free model gives two unequal Majorana velocities,
still total C=1, until one vanishes. It is inappropriate to infer a
uniform perturbative radius from the limit K->1.

As an explicit normalization check, Delta=1/2 has K=3/4 and dimensions
`25/12,3,16/3`. A nonzero a is therefore irrelevant for sufficiently
small a/b, but this calculation supplies no guarantee for a/b=1/2 or 1.
At Delta=-1/2, K=3/2 gives dimensions `13/6,6,8/3`.

The four-site duality extends the conclusion to a neighbourhood of
b=0 with

```
-1 < c/|a| < 1/sqrt(2),  c!=0,  |b/a| sufficiently small. (5)
```

These are useful controlled exclusions on a three-coupling face, not
proof that the whole face remains connected to them without another
critical sector or intervening phase.

## The remaining simple region

At a≈b the bare B-Majorana speed `4|b-a|` becomes small. The expansion
in c divided by that bandwidth is not uniform, even when all displayed
spin coefficients are comparable. At a=b,c=0 the B band is exactly
flat and extensively degenerate. Consequently the regular-free-CFT
argument cannot be extended across that point by continuity alone.

The earlier [C144](../../codex/work_C144/RESULT.md) and
[C157](../../codex/work_C157/RESULT.md) calculations analyze virtual
interactions near this flat-band point. Their second-order attraction
and dimerization mechanism constrain that asymptotic regime. They do
not give an all-order, order-one-c theorem; their own caveats preserve
that distinction. A two-Majorana lattice representation is not an
upper bound on the number of interacting continuum sectors.

Thus `H=2XY+ZZ` and `H=2XY+(1/2)ZZ`, equivalent to rational points
on a=b,h=0, are particularly simple *residual tests*. Their microscopic
monomial spreads are 2 and 4. In the original grouped coordinates
(a,b,c), the spreads are 1 and 2; these are different stated bases,
not contradictory metrics. We report no new evidence of C>1 there.
The right next diagnostic would resolve the complete low spectrum and
convergence of neutral modes, including a possible very small velocity.

## Why this does not settle the full term-necessity question

The a=0 and c=0 exclusions show that pairing and interaction are
necessary in this family. The b=0 staggered interacting strip and
the h=0 comparable-coupling region remain honest residual problems.
Neither the exact Ashkin–Teller name, a two-colour Ising interpretation,
nor a lack of earlier numerical excess proves that b and h are both
necessary. The general finite bound B(2) requires further input.
