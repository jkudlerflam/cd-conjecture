# Which terms are necessary in the simplest symmetric qubit family?

2026-11-17. OpenAI Codex, GPT-6 family. No independent review is claimed.
The central charge below is the usual total nonchiral Virasoro charge
`C` of the complete unitary compact conformal bulk theory; the exchange
coefficient is called `c`, to avoid confusing the two.

We consider exactly

```
H(a,b,c,h)=sum_j [a(XX-YY)+b(XY-YX)+c ZZ]_(j,j+1)+h sum_j Z_j,
```

with real fixed coefficients, physical qubits, nearest-neighbour range
and one-site translation symmetry. Constants and telescoping onsite
bond terms are removed; only a positive common energy scale is quotiented.
Site-dependent rotations below are solution methods with boundary seams
specified, not changes to the physical problem.

## The classification table

| Removed structure | Exact transformed model or structure | Exact charge / free structure | Infrared conclusion and its standard | Residual region |
|---|---|---|---|---|
| `a=0` | `|b|(XX+YY)+c ZZ+h Z`, with a spiral boundary twist | Uniform `sum Z`; ordinary XXZ in a uniform longitudinal field | **[Established exact-solution physics]** Every ordinary finite-velocity conformal regime has `C=1`; gapped and nonconformal saturation/ferromagnetic endpoints are not extra CFTs | No residual `C>1` region under this standard XXZ classification; `b=0` is classical |
| `b=0` | Alternate `X` rotations give `a(XX+YY)-c ZZ+(-1)^j h Z` | Staggered `sum (-1)^j Z` becomes uniform charge on even rings | **[Proved gaps]** for `h!=0`: `c<=-|a|` has gap at least `2|h|`; `c>sqrt(a²+h²/4)` has exact even-ring gap `4(c-sqrt(a²+h²/4))`. Boundary is nonconformal for `a!=0`. **[Derived / numerical]** `C=1` in the controlled liquid regions | The whole interacting strip `-|a|<c<sqrt(a²+h²/4)`, `c!=0`, is not an unconditional theorem. Prior small-field and near-boundary arguments, large-field assumptions and neutral-sector limitations remain explicit |
| `c=0` | Uniform quadratic Jordan–Wigner chain with imaginary hopping and real pairing | Free Bogoliubov fermions, physical fermion parity | **[Proved free-band classification]** `C<=1`; in this restricted family the ordinary gapless cases have `C=1`. Flat-band and quadratic band-edge limits are not additional conformal sectors | None in the free class; all endpoint qualifications are in REDUCTIONS |
| `h=0` | Two Majorana chains with hoppings `b+a,b-a`, coupled by bond-energy products: asymmetric self-dual Ashkin–Teller | Separate Majorana-species parities for even length, with boundary constraints; generally no additive onsite continuous charge | **[Proved algebra]** exact Ashkin–Teller mapping and `(a,b,c)->(b,a,-c)` duality. **[Derived, with hypotheses]** a stable `C=1` neighbourhood of the XXZ line when `-1/sqrt(2)<c/|b|<1`, `c!=0`, and `|a/b|` is sufficiently small; also its dual neighbourhood | Comparable couplings, especially the slow-Majorana region `|a|≈|b|`, are not excluded by the two-Majorana continuum argument. No unconditional all-face `C<=1` theorem is asserted |

The gap results and the distinction between uniform and staggered fields
come from the project's [exact staggered reduction](../d2_gap_reduction/staggered_slice.md).
The later [C153 analysis](../../codex/work_C153/RESULT.md) explicitly
labels its intermediate/large-field evidence as numerical/heuristic.
Its absence of a second observed mode is not promoted here to a theorem.
The zero-field mapping and duality were already recorded by
[C152, section 8](../../codex/work_C152/PROOF.md); we supply a local
four-site unitary, a boundary audit and a narrower infrared statement.

## Coupling count is not Pauli-string count

A common rotation by pi/4 about z gives

```
H -> sum_j [(a+b)XY+(a-b)YX+c ZZ]+h sum_j Z.
```

| Family | Independent real coefficients | After positive scale | Pauli monomials as written originally | After this common rotation |
|---|---:|---:|---:|---:|
| Generic four-term family | 4 | 3 | 6 | 4 |
| `a=0` | 3 | 2 | 4 | 4 |
| `b=0` | 3 | 2 | 4 | 4 |
| `c=0` | 3 | 2 | 5 | 3 |
| `h=0` | 3 | 2 | 5 | 3 |

The total field is counted once, not as two half-fields in a bond.
The coefficient map `(a,b)->(a+b,a-b)` is invertible. Four independent
rotation invariants recover `c,a²,b²,h²`; no generic member reaches any
coordinate face by a common rotation. This is an algebraic obstruction
to removing a structure, not an infrared bound.

For the earlier `r=40` reference `a=1/2,b=1/160,c=10`, the two rotated
cross couplings are `81/160,79/160`. The small DM coefficient has become
their small difference. Invariant ratios `|a/b|=80` and `|c/b|=1600`
remain; improved literal spread in another basis does not remove that
hierarchy. [NORMAL_FORM.md](NORMAL_FORM.md) gives the invariants and
the diagnostic comparison, including the field.

## What this tells the simpler-model search

Within this family, **nonzero pairing `a` and interaction `c` are necessary
for `C>1`**, using the established XXZ infrared input for the first claim
and the exact free classification for the second. Necessity of `b` or
`h` is **not proved**. The useful remaining three-coupling faces are the
staggered-charge face `b=0` and, especially, the field-free interacting
Ashkin–Teller face `h=0`.

The substantive analytic calculation here treats the latter near its
XXZ limits. It derives the symmetry-allowed vertex lattice, including
operators generated at second order, and specifies the resulting
`C=1` stability region and its uncontrolled boundary. See
[ZERO_FIELD.md](ZERO_FIELD.md). It uses perturbative continuum matching,
not a universal one-fluid assumption or a count of lattice Majoranas.

For a particularly simple **residual test**, impose `a=b` and `h=0`.
The common rotation then gives just `2a XY+c ZZ`, with two independent
coefficients before scale. Points `(a,b,c)=(1,1,1/2)` and `(1,1,1)` are
moderate rational tests of a still-unresolved strongly interacting
region, **not counterexample candidates backed by new CFT evidence**.
Their bare slow-Majorana hopping vanishes, so comparable microscopic
coefficients alone do not guarantee comparable low-energy velocities.
The duality identifies the two signs of c on this line.

## Reproduction and provenance

[REDUCTIONS.md](REDUCTIONS.md) gives explicit algebra and boundary conditions.
[SOURCES.md](SOURCES.md) records primary sources opened and the limited
claims used from each. `check.py` supplies exact symbolic bond identities
and independent complete finite-spin matrix checks. It tests open and
periodic geometries, the six-site twist rather than ignoring it, the
free many-body spectrum, and the normal-form invariants. No new DMRG,
central-charge fit, or independent scientific review is claimed.

```
OPENBLAS_NUM_THREADS=1 python3 contrib/qubit_minimal_terms_20261117/check.py
```

Dependencies are in `requirements.txt`; the small `checks.json` records
parameters, residuals and wall time. The general finite bound B(2), the
resonance hierarchy, and all uncontrolled neutral/composite sectors stay
open. This is a term-necessity classification with a conditional local
infrared result, not a new global bound.

The recorded check takes about half a second with one BLAS thread.
Exact symbolic identities pass. Complete finite-spin identities through
L=8 have errors below 6e-15; complete open free spectra through L=6
agree within 2.4e-14. Omitting the required six-site duality twist gives
an order-one error, so the boundary check is a genuine negative control.
