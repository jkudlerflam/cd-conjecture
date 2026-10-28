# C16 progress

Started 2026-10-27. Background, STATUS.md, project STATE.md, C6 parametrization,
and C8 qubit frustration-free witness read. No pre-existing C16 progress.
C15 and C17 had no progress files at the first read. All changes stay here;
no network, delegation, or external project-memory edits.

## Findings and active work

- [proved, logical] Compactness plus local boundedness only bounds c globally
  if local boundedness holds everywhere. It does not imply c <= 1, nor does
  it force an isolated finite-c > 1 point to approach a nonconformal point.
- [assumed] Applying the c-theorem to lattice neighborhoods requires uniform
  RG control and an appropriate relativistic continuum regime. Multiple
  velocities and competing ground-state sectors need separate treatment.
- [active] Derive the exact one-magnon dispersion where magnon number is
  conserved, distinguish product eigenstates from invariant magnon sectors,
  and derive dilute single-species and Pokrovsky-Talapov conditional bounds.
- [active] Analyze the FF hard-blockade parent sum n_i n_(i+1), which has
  exponentially many ground states and is not excluded by K_FF(2)=4.
  Compute its projected perturbations and check a solvable c=1 direction.
- [active] Give symmetry- and translation-resolved relevant-operator counts
  for candidate product CFTs, keeping actual reachability separate.
- [planned] RESULT.md will give a conditional reduction and state the missing
  exhaustion theorem rather than claim a finite classification without one.

## Step 1, exact singular-point diagnostics

- [proved] A product vacuum need not have an invariant one-magnon subspace.
  Mixed XZ/ZX terms can cancel on the vacuum and still split one magnon into
  two. The single-harmonic dispersion only controls the U(1) invariant case.
- [proved] On the free XY+DM locus there is a mixed linear/cubic band touching:
  E_-(k)=sin(k)-sqrt((cos(k)-1)^2+sin(k)^2) ~ -k^3/8 for k>0.
  Thus z=2 does not exhaust even the known nonconformal singularities.
  Free neighbors still obey c<=1 by the established quadratic theorem.
- [proved] The FF parent sum n_i n_(i+1) has Fibonacci-growth ground-space
  dimension. All projected uniform two-site perturbations reduce at first
  order to constrained hopping/current, chemical potential, and X/Y sources.
  K_FF(2)=4 does not constrain this ground-space dimension or its neighborhood.
- [numerical, in progress] checks.py checks the above identities and a
  number-conserving constrained-hopping direction on L=12,15,18,21,24 rings,
  using entropy and Casimir estimators. Python executable is python3.

## Step 2, completed checks and FF kernel classification

checks.py completed in seconds. checks_results.json contains all outputs.
At L=7 the magnon dispersion error is 5.33e-15; mixed XZ/ZX terms preserve
its product vacuum but give one-magnon leakage norm 0.973. At L=8 the exact
blockade projection error is 1.78e-15. The cubic-edge effective exponent
converges 2.99932, 2.99983, 2.99996.

For constrained hopping at density 1/3 on L=12,15,18,21,24 rings, the exact
energy formula agrees within 2.5e-14. Casimir estimates decrease from 1.01829
to 1.00452. Half-chain size-scaling gives c=1.01048 (leading fit), or
0.99750 with a 1/L^2 correction. Short-distance entropy-profile fits overestimate
c, from 1.17647 to 1.07398, and must not be read as c>1 evidence.

[proved] A useful exhaustive classification exists for the *ground-space
sizes* of positive uniform FF qubit bonds, not their full neighborhoods.
Rank-one entangled forbidden vector: open ground dimension L+1 via an
invertible transfer map to the symmetric subspace. Rank-one product forbidden
vector with proportional factors: Fibonacci growth. With independent factors:
open dimension L+1 and periodic dimension 2, via an invertible uniform basis
change. Rank >=2: C8's two-dimensional allowed-space classification bounds
long-chain ground-space dimension by 2 (or no FF extension); rank-three
allowed dimension one gives a repeated product if it extends. Nonunitary
basis changes preserve dimensions but do not prove spectral-gap equivalence.

[active] Final write-up separates these exact results from standard dilute-gas
and bosonization input, and from the unproved global exhaustion statement.

## Step 3, report written

RESULT.md now contains the complete derivations and the negative verdict on
an unconditional finite checkable reduction. It gives six research targets
with concrete tests, explicitly retaining the missing singular-exhaustion
and non-product-CFT cases instead of claiming they are excluded.

Symmetry counts are necessary continuum tuning counts, not microscopic
codimensions. Translation and the local operator algebra matter. For two
Ising factors, uniform order fields give 3 relevant primaries under a common
flipping Z2, 4 if only one flips, and 7 if neither flips. Derivative and
mixed-boson operators can increase these counts. A relevant-vertex search
is finite for a fixed positive stiffness matrix and discrete charge lattice,
but not uniformly finite as stiffness eigenvalues approach singular limits.

ff_checks.py completed on L=3,...,7 and agrees with every representative
kernel-dimension formula. ff_checks_results.json records the outputs.
Next: audit formulas and claim tags, inspect any newly available neighboring
progress for conflicts, then finalize without touching files outside C16.

## Final audit and handoff

Both scripts passed Python compilation. Formula and tag audit completed.
Read C15 and C17 progress after their reports appeared. Their cautions on
RG matching, conversion dimensions and emergent-mode counting agree with
this report. Added an attributed note on C15's pairing-only rank-one FF
factorization line, which fits the entangled row of the kernel classification.
No C15/C17 tests were rerun. The shared regular constrained gas and cubic
edge calculations had already been completed independently before the
neighboring progress appeared; no further duplicate work was launched.

Deliverables: RESULT.md, PROGRESS.md, local STATE.md, checks.py,
checks_results.json, ff_checks.py, ff_checks_results.json. No computations
remain running. The lack of a finite global reduction is a mathematical
obstruction identified in the report, not an unfinished numerical job.
No network, delegation, DMRG, commit, push, or outside-folder writes.
