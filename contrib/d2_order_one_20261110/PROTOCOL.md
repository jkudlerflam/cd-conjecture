# Comparable-coupling search: initial registration

OpenAI Codex, GPT-6 family, 2026-11-10. Registered before this contribution's
physics-screen output. Repository baseline: `3a1649c1` (also origin/main at
start). This is an initial, bounded local investigation, not a counterexample
claim. No cluster job is assumed available or submitted.

Target: one-site translation-invariant, nearest-neighbour physical qubits,
with a unitary conformal infrared theory of any c > 1.

## Preference and bookkeeping

Use Pauli normalization, counting the *total* on-site field once per site.
Define coefficient spread R as max/min of absolute nonzero coefficients in
the displayed Pauli expansion. Exact absent terms are excluded. Prefer R<=3;
3<R<=10 is a separately labeled exploratory band, not a user-specified cutoff.
R>10 bridge points are retained to follow the known mechanism, not called
comparable. Report exchange-only spread, expanded Pauli term count, grouped
term count, and each field or exchange tuning resolution. Overall rescaling
does not change R.

For the exchange matrix J report singular values, the Frobenius norms of
its symmetric and antisymmetric parts, and the field norm. These avoid hiding
hierarchies by rotations. Report small differences Jx+Jy separately even
when Jx and Jy are comparable. A field adjusted to criticality is permitted;
record its resolution and the cancellation in 4Jz+2hz explicitly. Small
generated energy scales and symmetry-forced zeros are not conflated.

## Lessons from the previous work

Read STATUS, REPORT 9.133, 9.147-9.150, 9.159, 9.172, 9.178-9.186, C190,
and the independent status audit. The r=40 candidate remains unestablished.
C168 found that an entropy/Casimir/gap combined detector rejects almost all
c=2 proxies; no Casimir or fitted dynamical exponent veto is used here.
Small periodic rings miss a dilute gas, have shell effects, and can select
different parities. High single-size entropy, an ordered cat state, a
Lifshitz point, or an unconverged MPS is not evidence for a CFT by itself.

## Stage A: full microscopic periodic-ring screen

Use independently bit-assembled sparse ED, both exact parity sectors,
initially L=10,12. Record two low levels per sector, residuals, ground-sector
selection, entropy at every cut, density/profile, and connected density
structure. Small rings only triage follow-ups; no c>1 verdict is possible.
Freeze these two grids before reading results:

1. Blockade continuation: Delta=Jx-Jy=1, Jx+Jy=0, r=4Jz,
   D=u/(4r). Use r in {1,3/2,2,3,4,8}, u in {1,3/2,2}, and hz in
   {-r/2+k/(4r): k=0,1,2,3,4,5,6,7,8}, omitting positive hz because spin
   reversal gives the partner. This intentionally uses a broad field ladder;
   the perturbative field estimate is unreliable at small r. Also evaluate
   the known reference r=40,u=1,hz=-19.9760472506 to measure this detector's
   blindness to the established finite-size candidate.
2. Attractive junction: Jx=1/2, Jy=-1/4, D=1/2, hz=0,
   Jz in {-1/4,-3/8,-1/2,-5/8,-3/4,-7/8,-1}.
   These rational bonds have R<=4 and explore strong pairing and attraction
   near the mechanism of the Peotta junction without importing its tiny Jy.
   Also include the existing exploratory reference phi=1, Delta=.7,
   V=-2.4, t=1 (Jx=(cos(1)+.7)/2, Jy=(cos(1)-.7)/2,
   D=sin(1)/2,Jz=-.6,hz=0).

Controls: nearest-neighbour XX (c=1), critical transverse-field Ising (c=1/2),
gapped Ising, and two decoupled XX sublattices (c=2, range two, **control only**).
For PBC report c_size=3 delta(S_half)/log(L2/L1), profile c=3*slope against
log[(L/pi)sin(pi*l/L)], and their window dependence. A finite-size excess
above 1.2 is a follow-up flag, not a pass or an exclusion threshold. Preserve
points of physical interest even if this detector misses them. Report all
grid outputs, not only high entropy examples.

## Stage B: bounded validation and adaptive promotion

After Stage A, register selected points before their new output. Use larger
rings (L=14,16, optionally18) to check shell/parity changes. Promote promising
comparable points and one continuation bridge to full microscopic open-chain
DMRG as resources allow, with the validated C190 parity solver. Validate its
new wrapper against independent small-L ED first. Suggested sizes L=24,48,
then96; final chi at least64 then128, with a larger-chi repeat of any putative
excess and independent initialization. Save sweep-boundary checkpoints in
/private/tmp/d2-order-one-checkpoints, never restart expensive work silently.

Use two OBC entropy averages fixed in advance: central three cuts, and the
trapezoidal average over [3L/8,5L/8]. Report c_size=6 delta(Sbar)/log(size ratio),
profile fits, both parity sectors, local residual <=1e-8, discarded weight,
last-sweep changes, and chi dependence. Aim for delta E<=1e-8 and delta S<=1e-3
before comparing chi; report failures without calling them negative physics.
Controls must distinguish c=1 from c=2 at the same size/window. Check central
density, full spatial profile, and correlations for ordering/clustering.
For a favorable screen require stable excess above1.2 on two size pairs and
windows, homogeneous bulk, control calibration and convergence; this still
does not establish a CFT. Independent spectral/dispersion or infinite-chain
evidence is needed for that claim. A c>1 region below1.2 remains possible.

Any follow-up beyond this registration must be labeled adaptive and written
before its outputs are inspected. Publish code, controls, raw numerical
records and null outcomes in this directory; avoid unrelated shared edits.
