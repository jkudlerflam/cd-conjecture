# Attractive-route continuation

OpenAI Codex, GPT-6 family, 2026-11-10.

The preferred-ratio entropy flag reproduces at a nearby exact rational point,
but its reference drifts toward an ordered crossover at length 192. This
branch does not establish a new unitary CFT candidate. The older exploratory
reference retains entropy slopes above one through 192, with substantial
ordering and finite-size-energy caveats. These are numerical observations,
not exclusions of these models at larger scales.

## Models and fixed coefficients

All models are uniform nearest-neighbor physical-qubit chains
H=sum[Jx XX + Jy YY + Dz(XY-YX) + Jz ZZ], with exactly absent on-site field.
Each has five expanded Pauli strings and four grouped exchange terms.
No blocking, projection, or enlarged on-site Hilbert space is used.

| Name | (Jx,Jy,Dz,Jz) | Literal spread | Exchange singular values |
|---|---|---:|---|
| `attr_ref` | ((cos(1)+.7)/2,(cos(1)-.7)/2,sin(1)/2,-.6) | 7.76656 | .85,.6,.15 |
| `pref_ref` | ((cos(1.2)+.8)/2,(cos(1.2)-.8)/2,sin(1.2)/2,-.5) | 2.65595 | .9,.5,.1 |
| `pref_rat` | (14,-5,11,-12)/24 | 2.8 | .891036,.5,.0993693 |
| `pref_zrat` | (7,-3,6,-9)/12 | 3 | .943713,.75,.110380 |

The singular spreads are respectively 5.667, 9, 8.967, and 8.550. Thus the
preferred literal ratios coexist with moderate transverse cancellation;
that invariant scale difference must not be hidden. The rational points
have fixed small-integer ratios with denominators 24 and 12, respectively,
not numerical parameter fine tuning. Both contain four distinct nonzero
coefficient magnitudes.

## Completed entropy comparisons

The table selects the lower variational-energy parity state only after
both parity optimizations finish. Every selected state is even. These are
lowest states found, not rigorous global-ground-state certificates.
The estimates are 6 deltaS/log2, using the central three cuts and the
trapezoidal fractional window [3L/8,5L/8].

| Point | Sizes | Central-three slope | Fractional slope |
|---|---:|---:|---:|
| `attr_ref` | 48 to 96 | 1.10133 | 1.13328 |
| `attr_ref` | 96 to 192 | 1.29007 | 1.31081 |
| `pref_ref` | 48 to 96 | 1.36483 | 1.42149 |
| `pref_ref` | 96 to 192 | 1.09192 | 1.14245 |
| `pref_rat` | 48 to 96 | 1.37077 | 1.42759 |
| `pref_zrat` | 48 to 96 | 1.20608 | 1.26586 |

The preferred reference's single-size profile estimates fall from
1.31588/1.52652 at 96 to .82331/.70862 at 192. Its parity splitting falls
.0426376, .0113346, .00114894 at lengths 48, 96, 192. XX(L/2) grows
.156824, .180551, .221077, while normalized central-window staggered-X
order decreases .326403, .282950, .267526, with a slowing decline.
Density remains .5 and cannot by itself reject a symmetry-restored cat.
These trends favor an ordered crossover over a stable above-one CFT.

Independent saved-MPS diagnostics decompose its 96-to-192 fractional
slope into .319231 binary Schmidt-parity entropy and .823221 conditional
entropy. The binary contribution approaches its bound log2: at 192 it is
.686223, so its contribution on another size doubling is at most about .060.
Conditional entropy is not a separate physical central charge. The parent's
independent periodic rings 14/16/20 also give entropy estimates near/below
one; see ../preferred_rings.jsonl.

The R=3 reserve instead has strong Z-order evidence: ZZ(L/2) grows
.244171 to .305398, while central-window uniform-Z order decreases only
.406341 to .386568. Its 96-site profile values are 1.04120/1.05025.
Its parity gap is .137188 at 96; the higher odd sector's entropy must not
be substituted for the even state's entropy. Z ordering preserves the
optimized parity, so binary-parity subtraction does not remove this caveat.

The exploratory reference retains 96-to-192 slopes 1.29007/1.31081 and
profile estimates 1.17874/1.19955. Its whole-window staggered-X order falls
.244717 to .214373, despite XX(L/2) rising .121892 to .140484. The
fractional slope decomposes as .349998 binary plus .960807 conditional.
A post-hoc energy diagnostic also shows substantial finite-size corrections:
D(L)=E(4L)-3E(2L)+2E(L) cancels bulk and boundary constants. Its effective
correction power from 24/48/96/192 is .3899, versus .9497 for exact XX on
the same sizes. This is not a dynamical exponent measurement or an exclusion
of eventual conformal scaling; see energy_scaling.json.

## Numerical audit and deferred scope

The unchanged, previously validated C190 solver is called through the
full-microscopic wrapper in ../../d2_order_one_20261110/dmrg_screen.py
(the actual import and source hash are in run.py and run_inventory.json).
Sixteen new runs cover eight model/size configurations, each with both
parities. All registered sweep gates, parity checks, canonical-versus-sweep
entropy checks, solver hashes, and checkpoint/JSON metadata checks pass.
The new records retain 4311.43 seconds of optimization sweeps. Previous
reference 48/96 files are read-only and separately identified as baselines.

All new lengths 48/96 use chi128; the two 192 references use chi192.
Selected-state entropy changes from chi64 to128 at preferred-reference96
and rational96 are 3.11e-5 and 2.80e-5. Changes from chi128 to192 at
preferred-reference192 and exploratory-reference192 are 5.09e-6 and
1.36e-5. These are much smaller than the observed size trends. Completed
schedules and local residuals alone do not certify bond-dimension convergence,
global optimality, or a CFT. This branch did not independently compute
variances; the parent's full-MPO measurements, with record/checkpoint hashes
and their limitations, are retained in ../variances.json.

Following the new 192 reference evidence, the unstarted `pref_rat` 192 pair
and later preferred96 chi192 refinements were explicitly deferred. Both
waiting 192 children were interrupted before any variational sweep, raw
record, or checkpoint existed. The earlier exploratory `attr_rat` comparison
was deferred before output when preferred-ratio points took priority.
No reserve192 promotion or alternate boundary-pinning solver was started.
PROTOCOL.md preserves the registrations and reasons; sentinel files prevent
accidental relaunch. No numerical jobs remain active, and no commits or
pushes were made by this branch.

Run analyze.py, audit.py, then plot.py to reproduce summaries, the audited
inventory, and diagnostics.png/SVG. Independent order/decomposition data
are in ../diagnostics/attractive_diagnostics.json and include input hashes.
The figure uses completed parity-selected data and hash-matched diagnostics.

The mechanism caveat is consistent with the primary literature:
[Peotta et al., section 3](https://arxiv.org/html/1403.4568v3) find first-order
Luttinger-liquid/ferromagnet boundaries and leave their meeting region
unresolved. That motivates coexistence and size controls; it does not assign
these particular numerical points to a phase by itself.
