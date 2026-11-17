# Literal model, simplicity and infrared cautions

## Exact source convention

Use Pauli matrices, not spin operators divided by two, and a down spin
as a particle: n=(1-Z)/2. The physical open-chain Hamiltonian is

```
H = sum_(j=0)^(L-2) [(XX-YY)/2 + D(XY-YX) + (V/4)ZZ]_(j,j+1)
    + hz sum_(j=0)^(L-1) Z_j,
D=kappa/(4V),  hz=-V/2+a/V,  Delta=1.
```

The literal C190/C195 symmetric production point has
V=40,kappa=1,hz=-19.976047250613814, hence
a=0.9581099754474565. The `hz0=-19.97625` in C190/models.py is an
untuned reference. STATUS's rounded field is not the production input.
The old C177 points have separately fitted D and hz, recorded in
`source_audit.json`; their kappa is not exactly one.

Writing the onsite term as hz(ZI+IZ)/2 in each bond requires adding
hz Z/2 at **each endpoint** of an open chain. Our MPO puts hz Z on every
site directly. The ED validation detects the missing-half-endpoint
error with operator norm |hz|. The computational energy shift
-(V/4+hz) per site is removed from physical energies and never counted
as a coupling.

The C190 fermion parameters are w=2iD, Delta=1, interaction V and
mu=V+2hz=2a/V. The leading blockade scale and dimensionless chemical
potential are

```
|t2|=1/(2V),   |w|/|t2|=kappa,   mu_eff=4a-6.
```

Thus fixing a=24/25 fixes mu_eff=-2.16 across V and L. It is close to
the production point's -2.1675601, without inheriting a long decimal
field. It is not an exact density constraint: density changes with V,
L, parity and virtual-pair dressing. In occupation notation the full
open-chain field also produces an endpoint potential V(n_0+n_(L-1))/2
relative to the uniform bulk chemical potential. That large endpoint
cost helps explain why small-chain densities require care. It is part
of the stated uniform spin Hamiltonian, not a boundary modification.

## Two counts and two explicitly stated bases

Before fixing scale the family has four independent real parameters:
Delta,D,Jz,hz. After Delta=1 it has three continuous ratios, represented
by V,kappa,a. A chosen point has fixed values, not fewer independent
structures merely because its coefficients are rational.

In the common original Pauli basis there are six nonzero monomials:
XX,YY,XY,YX,ZZ,Z. The physical field is counted once. At every point,
the same common onsite U=exp(-i pi Z/8) gives the exactly equivalent form

```
U^tensorL H U^dagger tensorL
  = sum_j [(1/2+D) X_jY_(j+1)+(1/2-D)Y_jX_(j+1)+(V/4)Z_jZ_(j+1)]
    + hz sum_j Z_j.
```

This has four monomials and the same four independent coefficients
before scale. The identity follows directly from
`UXU^dagger=(X+Y)/sqrt(2)`, `UYU^dagger=(Y-X)/sqrt(2)` and
`UZU^dagger=Z`: `XX-YY` becomes `XY+YX`, while `XY-YX` and `ZZ`
are invariant. The two transverse coefficients have independent sum
Delta and difference 2D, so the parameter count does not fall.
`simplicity.py` verifies the literal conjugation to
2.3e-16. The sibling
[normal-form analysis](../qubit_minimal_terms_20261117/NORMAL_FORM.md)
independently supplies the rotation invariants and cautions against
confusing this shorter expression with the removal of a structure.

For kappa=1,a=24/25 and V>=4, the original literal spread is
R=2V^2-4a; the rotated spread is |hz|/(1/2-1/(4V)). At V=8 these
are 124.16 and 8.27733; at the old V=40 production point they are
3196.17 and 40.45782. Constants, telescoping bond terms and a common
energy scale are excluded. The smaller rotated spread encodes a
small difference between transverse exchanges: their difference is
2D=1/(2V). Invariant ratios (Delta/2)/D=2V and Jz/D=V^2 still distinguish
the hierarchy. The second-order |t2| is a generated kinetic scale;
neither |t2| nor a mode velocity is used as a substitute for the
microscopic coefficient-spread calculation.
No mode velocity is measured in this campaign. The tabulated `t2`
is a leading perturbative hopping coefficient, not a measured velocity
or a lower bound on the smallest infrared scale.

In particular, the rational V=8 point is exactly

```
sum_j [(XX-YY)/2+(XY-YX)/32+2ZZ]_(j,j+1) - (97/25) sum_j Z_j,
```

or, after the stated common rotation,

```
sum_j [(17/32)XY+(15/32)YX+2ZZ]_(j,j+1) - (97/25) sum_j Z_j.
```

At V=12 the rotated exchanges are 25/48,23/48, Jz=3 and hz=-148/25.
The older low-hierarchy V=4,a=1 point has rotated exchanges 9/16,7/16,
Jz=1 and hz=-7/4. That point's unresolved anomaly is documented in
[the prior comparable-coupling search](../d2_order_one_20261110/README.md).

## What no finite-size entropy test here settles

The symmetric family has the same potentially dangerous resonant
two-pair operator discussed in C179 and later REPORT updates. The
new parameter perturbations remain on Jx+Jy=0; they test only this
four-coupling family, not stability to the extra XX+YY coupling or
parity-breaking fields. A neighbourhood of finite-size gas signatures
inside that slice is not an open c>1 region of general qubit bonds.
The
old effective-model pseudospin evidence does not automatically remove
charge locking in the full qubit chain. At V=8 the older infinite-chain
reading near 1.4 is therefore retained as unresolved tension with a
clean two-mode interpretation, not dismissed because a new finite-chain
slope is higher. It uses
different fitted coefficients and an older solver, so it is also not
a direct verdict on the new rational point.

Balanced valley weights do not distinguish a two-component liquid
from a paired state. A parity gap on a few open lengths may have shell
effects, and a small local solver residual does not bound the global
infrared gap. Connected pair correlations are recorded without fitting
an exponent from a few short distances. The one- and two-leg XX controls
diagnose solver bias and finite-size estimator effects; they do not
turn every matching entropy slope into a CFT identification.

The term-removal sibling proves neither that DM nor the field is
necessary in every possible c>1 qubit model. Its zero-field residual
Ashkin--Teller/slow-Majorana region remains open. Those simpler residual
tests currently have less direct c>1 evidence and are not promoted to
counterexamples by their short expressions.

## Observable definitions and a quadrature check

The V=8 tolerance probes change a by +/-1/200 or kappa by +/-1/16,
one parameter at a time. These mean physical changes
`delta hz=+/-1/1600` (0.0161% of |hz|) or `delta D=+/-1/512`
(6.25% of D). In the rotated expression the latter changes each
transverse coefficient by only about 0.37--0.42%, in opposite directions.
Finite-L response to these probes must not be advertised as a broad
verified critical phase, whichever basis is used to display the bond.

The source convention uses the central-window fermion matrix with its
Jordan--Wigner string and a 256-point momentum grid. It partitions
weights by cos(k)>0 and cos(k)<0. Floating-point cos(+/-pi/2) is positive,
so both mathematical boundary bins enter the first group. This produces
an apparent imbalance about 0.0016 at the symmetric candidate.
`analysis.json` preserves the source result as `valley` and also reports
`valley_half_boundary`, assigning half of those two bins to each group.
For the completed V=8,L=128 states the corrected imbalance is below
3e-16. This is a quadrature correction, not evidence for two gapless modes
or against a symmetry-breaking thermodynamic phase hidden by a finite
ground-state superposition.

The additional connected pair diagnostic is
<b_i^dagger b_(i+2)^dagger b_j b_(j+2)> minus the product of the two
pair expectations, where b=|up><down| is the physical hard-core-boson
annihilator. It is a spin-pair observable without an added fermion
string. Both real and imaginary parts are retained, at explicitly
listed separations; no short-distance exponent is assigned to it.
