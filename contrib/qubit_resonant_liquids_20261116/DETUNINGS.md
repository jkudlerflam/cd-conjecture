# Admissible hopping detunings and a fixed microscopic test point

Date: 2026-11-16. Status: **derived leading dispersion; numerical
few-particle evidence; unvalidated microscopic test point**. The results
below do not assign a Virasoro central charge.

## 1. Smallest tested change that mixes residue labels

Retain alpha=beta=f6=0 and add tau to the sixth microscopic hopping
coefficient. The effective operator is

```
P[K_m + tau T1 + tau* T1† - mu6 N]P.
```

This is an order-six, nearest-neighbour microscopic parameter change that
leaves the order-two and order-four constraints intact. `candidate.py`
checks the complete sixth matrix for complex tau at q=1 and real tau at
q=3, with a chemical shift, on L=10,12,13. The maximum error is
5.69e-14. This does not assume particle-number conservation in the
microscopic elimination; N conservation holds only in this leading slice.

Write tau=t gamma exp(i phi). In units of t the one-particle dispersion is

```
e(k)/t = -2 cos(m k)+2 gamma cos(k-phi).
```

A generic phase selects one lowest valley. Two useful leading slices are
q=1, phi=pi/4, and q=3, phi=0, both with gamma>0. They have two equal lowest
minima in the tested gamma range 0.02,0.1,0.3. At q=3 the minima obey

```
k = +/- arccos(-sqrt((3+gamma)/12)),
e_min/t = -2(3+gamma)^(3/2)/(3 sqrt(3)).
```

At q=1 their sum is -3pi/2 modulo 2pi; at gamma=0.1 they are
-3.1371924942 and -1.5751964862, with e_min=-2.1417324872.
The two equal minima at q=1 follow from the leading-band reflection
k -> pi/2-k. We have not identified an exact microscopic symmetry that
protects that equality. At q=3, real microscopic coefficients preserve
spatial inversion combined with product Z, which sends k to -k. The
facilitated F term changes sign under each factor separately.

## 2. Two- and three-particle tests

`relative_detuning.py` reduces the infinite-line two-body problem at fixed
total physical momentum K to a relative-distance half-line. In units of t,
the off-diagonal relative hops are

```
r -> r+m: -2 cos(m K/2),
r -> r+1: 2 gamma cos(K/2-phi).
```

The constraints and short-range potentials are exactly those of K_m.
For q=3, r=4 is excluded and G8 becomes a potential 8/t at r=8 in the
two-particle sector. There is no crossing term: an m-site hop cannot pass
another particle at the allowed minimum separation m. We impose a hard
wall at rmax=512,1024,2048. The continuum edge is minimized separately
over relative momentum. The cutoff is a relative-distance box, not a
microscopic system size.

For both q and all three gamma values, at K=2k_a and K=k_a+k_b, the lowest
box energy lies **above** the corresponding continuum. The differences
decrease approximately as rmax^-2; short-distance probability decreases
and mean separation grows approximately as rmax/2. At gamma=0.1:

| q | Channel | E_box-E_cont, rmax=512 | rmax=1024 | rmax=2048 |
|---|---|---:|---:|---:|
| 1 | same valley | 0.00120376 | 0.000301719 | 0.0000755259 |
| 1 | opposite valleys | 0.00118155 | 0.000299549 | 0.0000751557 |
| 3 | same valley | 0.000677276 | 0.000169816 | 0.0000425158 |
| 3 | opposite valleys | 0.000661880 | 0.000167931 | 0.0000423084 |

All energies in the table are divided by t. These are evidence against
a localized pair at the tested thresholds, **not a proof of absence of
arbitrarily shallow binding or of binding at every other total momentum**.

Independent finite-ring occupation-space diagonalization in
`scattering.py` finds E2-2E1>0 for L=48,72,96 at all three gamma values,
and E3-3E1>0 for L=24,36 at gamma=0.1. Here E1 is the actual finite-ring
minimum, not the infinite-line edge. At L=96,gamma=0.1 the two-particle
differences are 0.02964169 and 0.01475762 in units of t for q=1,3.
These few-particle tests neither determine the relative-valley stiffness
nor rule out valley polarization, a finite-density cluster phase, or an
infrared gap. No entropy or Casimir fit is made.

## 3. Pairing distinguishes the two slices

At q=3, real tau gives k_a+k_b=0. Creating an opposite-valley pair is
translation invariant and even under inversion times product Z. Thus the
microscopic symmetries do not forbid it. Its coefficient and scaling
dimension must be checked before interpreting two leading minima as two
surviving modes. This symmetry observation is not a calculation of a
nonzero pairing coefficient or a proof that it is relevant.

At q=1, phi=pi/4 the opposite-valley pair instead carries momentum pi/2
modulo 2pi. However other number-changing processes have small nonzero
momentum mismatches. At gamma=0.1 the smallest centre-momentum mismatch
among net creation vectors (a,b) with 1<=|a|+|b|<=4 and a+b!=0 is
|k_a+2k_b| mod 2pi = 0.0044001594. Same-valley pair creation at k_a has
mismatch 0.0088003188. `candidate.json` records the complete enumeration,
keeping only a+b>0 to omit Hermitian-conjugate duplicates.

These are **bare valley-centre kinematics**. They are not an interacting
bosonization operator table: density harmonics, dressed momenta, scaling
dimensions and coupling strengths remain necessary. In particular,
nonzero mismatch alone does not ensure that a potentially relevant
coupling is detuned before it reaches strong coupling.

## 4. Explicit size-independent microscopic test point

The following is a concrete point for further analysis, not a verified
counterexample. It favors the q=1 slice because its opposite-valley pair
has a substantial momentum mismatch. Define exactly

```
epsilon=1/1000, gamma=1/10, phi=pi/4, nu=1/20,
tau=gamma exp(i phi),
e_min=min_k[-2 cos(4k)+2 gamma cos(k-phi)],
mu6=e_min+nu,

H = D1+epsilon F+u N+[h T1+p A+h.c.],
u=2 epsilon²-10 epsilon⁴+(50-mu6)epsilon⁶,
h=-epsilon²+4 epsilon⁴+(-24+tau)epsilon⁶,
p=epsilon².
```

The convention includes h.c. only for h T1+p A. This is a Hermitian,
one-site translation-invariant, nearest-neighbour qubit chain, with no
exact U(1) assumption. Every coefficient is fixed independently of L.
The scalar counterterms can be omitted. `candidate.py` exports its
literal 4x4 bond, including real and imaginary parts, and verifies
Hermiticity. Numerically mu6=-2.091732487236675. The exact definitions
above, rather than rounded decimal coefficients, specify the point.

The chemical coefficient lies nu above the leading band minimum. It is
a trial filling choice; an interacting density or even a two-valley
ground state has not been established. The small epsilon reduces the
formal relative size of the displayed seventh-order terms. No uniform
remainder estimate proves that this particular epsilon is small enough.
The sixth kinetic scale is 1e-18, so direct microscopic energy differences
in ordinary double precision are not a trustworthy validation method;
use a controlled effective calculation or adequate precision.

## 5. Decisive remaining tests

1. Determine the leading finite-density equation of state and the
   relative-valley stiffness at this tau. A positive two-particle binding
   difference alone does not show a stable mixture. Extract both
   intra-valley and inter-valley scattering data or an equivalent
   many-body susceptibility, and test phase separation/polarization.
2. Compute the energy-dependent elimination beyond order seven, including
   its metric, then the coefficients of the allowed number-changing
   operators. The seventh static residue-changing witness in
   `DERIVATION.md` is not itself a new conformal velocity; its projection
   on the undeformed zero manifold vanishes on the tested rings.
3. At a fixed microscopic point, test the competing running coupling and
   momentum-mismatch scales using dressed low-energy parameters. If two
   positive velocities survive, obtain converged independent central-
   charge estimates with the same microscopic coefficients as L grows.

Until those tests pass, the deliverable is an exact leading obstruction,
a correction to the alpha-only detuning argument, and a reproducible
candidate that passes the stated few-particle checks. It is neither a
microscopic c>1 result nor a universal no-go theorem.
