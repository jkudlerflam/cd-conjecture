# The fifth and sixth resonances, and what detuning remains open

All microscopic paths in this note use the exact polynomial parent (1)
of [PARENT_THEOREMS.md](PARENT_THEOREMS.md), with r the small parameter.
They are not the old path with constant `p/epsilon²` in the v=1 chart.
The numerical implementation retains the full spin basis before each
elimination and every sector of each subsequent constrained space.

## 1. Correct lower constraints and their transfer zeros

At order r² the constraint is `2D2` inside the no-11 space. At order
r⁴, on P2, it is

```
(2/3)(q-3)²D3+(q-1)²D4.                                  (17)
```

The exact parent has no residual chemical or hopping term at these
orders. For real q different from 1,3 this selects P4. The sixth-order
local calculation gives

```
K6=P4 [w5 D5+w6 D6+w7 D7+w8 D8] P4,
w5=(q²-10q+5)²/3,
w6=(3q-1)²(q²-6q+25)/12,
w7=(q³-21q²+35q-7)²/[4(5q²-18q+21)],
w8=(q²-6q+1)²/8.                                        (18)
```

**Status:** derived local Schur identities, verified as complete symbolic
matrices on an open nine-site interval; larger full-sector matrices check
the two resonances with spectators. Equation (18) is not obtained by
fitting just a one-particle band. All four weights are nonnegative for
real q: the two nonsquare factors in the denominators/numerators are
strictly positive. The chart excludes q=1,3 even though some of the
rational expressions extend continuously there; those points require a
different fourth-order projection.

The vanishing factors in (18) are those of the transfer polynomials
s5,...,s8, with factors associated with the earlier resonances removed.
In particular:

* q5=5-2sqrt(5): w5=0 and w6,w7,w8>0;
* q6=1/3: w6=0 and w5,w7,w8>0.

Thus the correct sixth-order zero space forbids **all separations from
one through eight except m**. It is not simply the hard-rod space with
minimum separation m. The other forbidden separations matter in the
subsequent denominators. The code checks the complete occupation-word
description against the actual zero-diagonal space at each size.

The odd coefficients through order seven vanish in the complete tested
spaces. This statement includes potential changes of particle number.
No fixed-N approximation is used.

## 2. Exact interval coefficients and larger checks

The first nonzero retained dynamics for m=5,6 in these calculations is
at order r⁸. On the minimal open interval of m+1 sites the exact endpoint
hop coefficient is

```
-t5,  t5=20992-(28160/3)sqrt(5)
        =2048/[3(123+55sqrt(5))] ~=2.77525120197405,
-t6,  t6=832/2187 ~=0.380429812528578.                      (19)
```

`exact_local.py` derives these in exact algebraic/rational arithmetic,
and reproduces the earlier m=3,4 sixth-order coefficients as controls.
An endpoint calculation alone would not establish a uniform bulk
coupling. The separate larger checks are therefore part of the result:

| m | Open intervals | Compatible periodic rings |
|---|---|---|
| 5 | L=10,13,15 | L=10,15 |
| 6 | L=12,14,16 | L=12 |

Every retained matrix entry is compared with a configuration graph whose
off-diagonal element is exactly -t_m for each admissible hop by m,
including multiplicities. This tests interior hops, endpoint hops,
spectator particles, and the three-particle sectors at the larger open
sizes. There is no extra correlated hop or number-changing entry in the
tested eighth-order matrices. The error is below 1e-12.

On these complete finite spaces the eighth coefficient has the form

```
K8=t_m L_graph+V_diag,    V_diag>=0,                       (20)
```

to numerical precision. Every single-residue block has V_diag=0 and is
exactly the corresponding coarse ferromagnetic SWAP Laplacian. All tested
full kernels have dimension L+1, agreeing with the exact parent limit.
The single-residue interaction is therefore at `V=-2|t_m|`, not an
independently repulsive gas. All potential and off-diagonal tests are
saved in `checks.json`.

**Scope:** (19) is exact on the stated small intervals, while its uniform
bulk interpretation is supported by the larger interval/spectator and
compatible-ring comparisons. This contribution does not claim an
all-length linked-cluster proof of every eighth-order diagonal term,
nor an all-m hopping-order theorem. The exact all-m parent and
ground-space theorems do not depend on that extrapolation.

Incompatible short rings were also inspected during development. Their
winding processes can generate additional eighth-order one-particle
entries. For example, m=5,L=11 has nonzero apparent hops of lengths three
and four, whereas the larger open intervals and compatible rings have
only the five-step graph. Those short-ring entries are not reported as
new bulk hopping terms. This is why boundary and compatibility checks
are essential here.

## 3. What is and is not an all-order conclusion

The recurrence (8)--(10) is exact for every distance, and supplies
arbitrarily large exact periods at fixed local dimension two. The
ground-state pole coincidence at distance m is also exact. Therefore
finite microscopic parameter dimension cannot exclude arbitrarily large
resonance periods.

The calculations establish first motion at order six for m=3,4 and at
order eight for m=5,6. They do **not** establish an all-m rule such as
`2 ceil(log2 m)+2`, even though growing constraint radii might suggest
testing such a rule. In particular the transfer recurrence does not
determine the energy denominators, all allowed correlated moves, or
which lower constraints must be eliminated at every later order.

There is, however, an all-m obstruction to interpreting the *exact*
parent's enlarged kernel as m independent conformal sectors: its exact
ground-space degeneration is the sum (14), its rank is L+1 rather than
an extensive label-code dimension, and its positive spectral gap has
the quadratic upper bound (15). A nonunitary ground-space mapping is
insufficient to transfer an excitation spectrum to a Heisenberg model.

## 4. Explicit normal detunings and their limits

At order r⁸, independent nearest-neighbour counterterms have direct
compression

```
delta K8=P [mu N+tau T1+tau* T1†+h B+h* B†] P,             (21)
```

where P is the correct sixth-order zero space. Adjacent pairing and
facilitated terms have zero direct P compression; their effects at
other orders must still be included if they are inserted earlier.
Equation (21) is a set of admissible physical normal directions,
not an exhaustion of detunings invisible at earlier scales.

The chemical direction alone does not turn the ferromagnetic threshold
into an intermediate-density liquid at the leading order. On a
compatible ring the maximum allowed N is L/m, and the vacuum and each
full single-residue crystal are zero modes of the positive parent
coefficient. Thus `K8+mu N` selects the vacuum for mu>0 and the saturated
crystal sector for mu<0, by positivity and the elementary bound on N.
This is a leading-order statement; no uniform microscopic stability
estimate is silently assumed.

The tau term moves a particle by one and changes its residue. It can
couple the ordered-label sectors, and the h term changes particle
number. Their interacting phase problem is not decided by the m minima
of the unperturbed one-particle compression. Earlier invisible jets,
such as a change in the pair amplitude of order r³, must be reduced
together with their induced number-changing companions. They cannot be
treated as an independent density repulsion merely from a diagonal
matrix element. The sibling resonant-liquids investigation found just
this issue for m=3,4; it is not extrapolated here as an all-m theorem.

The fixed finite detunings that might retain several positive-velocity
sectors are therefore still open. Nothing here proves their existence
for arbitrarily large m, or proves that every detuning reduces to a
single collective mode. The exact collective ground-space result does
not control the excited label sectors once a detuning is applied before
the thermodynamic limit.

## 5. The local-boundedness lemma still required

To obtain a finite bound near the blockade one would need an open
normalized-bond neighborhood U and a finite B such that every fixed bond
in U with a complete unitary conformal bulk limit has total c<=B.
The bound must hold uniformly over m, allowed detunings, shrinking
velocities and crossover scales, and must include neutral/composite
sectors. A theorem that all normal detunings of these resonances admit
a bounded complete continuum description would address part of that
requirement. The parent ground-space identities and variational gap
bound do not supply it.
