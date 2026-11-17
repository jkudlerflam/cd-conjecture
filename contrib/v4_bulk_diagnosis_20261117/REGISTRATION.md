# Bounded V=4 bulk-phase test

Written 2026-11-17 before new production in this contribution. Reusing
existing profiles and checkpoints is explicitly allowed. Author: OpenAI
Codex, GPT-6 family. The thermodynamic central charge is UNKNOWN.

## Fixed model and question

The original common-basis bond is `(XX-YY)/2+(XY-YX)/16+ZZ`, with
the uniform onsite field `-7/4 Z`, including both endpoints. The same
onsite `exp(-i pi Z/8)` rotation gives `(9/16)XY+(7/16)YX+ZZ` and
the unchanged field. Four independent real couplings become three
ratios after scale fixing; six original monomials become four in this
common rotated basis. Spreads are 28 and 4 respectively. Pairing/DM=8,
ZZ/DM=16 and the generated leading t2=1/8 remain distinct measures;
no mode velocity is known. All bulk coefficients stay fixed across sizes.

The question is whether the L=96 density range near 0.301 reflects
boundary oscillations, density/bond order (including a symmetry-restored
state), or an unresolved crossover. No one-size entropy fit will assign c.

## Executed set and resource cap

1. Reuse both saved L=96,chi=256 states from the previous contribution.
   Check checkpoint metadata and archived energies before measuring.
2. Validate the new connected density, bond and pair contractions against
   literal small-chain ED in both parities before larger measurements.
   Calibrate density metrics on exact open one-/two-leg free controls
   at density near 0.296, and on ordered product/cat examples. These are
   diagnostic controls, not new nearest-neighbour counterexamples.
3. Benchmark fixed work at L=144,chi=256 with the actual V=4 bond before
   new production. Proceed only below 8 GB per process and a projected
   10 minutes per sweep, using the prior measured full-sweep times as
   the primary calibration. If not met, report the resource limitation.
4. Main sequence: L=48,96,144, both parity sectors. New L=48 schedules
   end at chi=128,192. New L=144 schedules end at chi=128,192,256.
   Save each unmixed stage's checkpoint, and measure the last two chi
   levels. Existing L=96,chi=256 states supply the middle size.
5. One genuinely independent seed-7 L=96 even run through chi=192,256
   checks the state whose modulation motivates this test. This is in
   addition to the previously agreeing independent chi=192 calculation.
6. Only if L=144 fails the last chi comparison, permit one continuation
   to chi=320 in each affected sector, after its fixed-work cost benchmark.
   No chi above 320, no size above 144, no field scan in this round.
7. Small periodic rings L=8,12,16 in both parity sectors, at most four
   low states per sector, provide a boundary-independent comparison of
   density/bond structure and low gaps. Translation-averaged observables
   must be identified as symmetry-restored averages. Small-ring gaps
   cannot establish a positive thermodynamic velocity or a CFT.
8. A boundary probe is optional only if the original-boundary modulation
   amplitude at 144 exceeds 85% of that at 96 while normalized connected
   peak weight falls by at least 15%. In that case, allow at most the two
   L=96 parity states with HALF endpoint fields, explicitly labeled as
   a different boundary Hamiltonian, using chi=192,256. Do not mix them
   into the original-boundary size sequence. No other boundary scan.

At most TWO active computational workers (production, measurement or ED),
one BLAS thread each; no overlapping extra validation job. Checkpoints are
ignored and outputs are small JSON records. Commit completed checkpoints
of the research, not incomplete state files.

## Numerical gates and stopping rules

Retain C190-style local gates: final residual <=1e-8; last energy change
<=1e-8, max entropy-window change <=1e-3 and density change <=5e-4;
discarded weight <=1e-7; last chi-step Sbar6 change <=0.01. The stricter
within-stage stopping rule is delta E<2e-9, delta S<5e-4, delta n<2e-4,
residual<1e-8, after at least three unmixed sweeps, at most six. Preserve
limit hits. Mixer stages are initialization, not converged comparisons.

Also compare the final two chi levels' bulk profile RMS and connected
structure-factor peak divided by window length. Require change <=
max(5% of the larger magnitude,1e-4) for each density and bond diagnostic
to call that observable locally stable. A failed entropy gate stays
visible even when these observables stabilize. Independent-start checks
use |delta E|<=1e-8 and the same observable tolerance; failure leaves
possible metastability unresolved. These are local numerical gates,
not statistical confidence intervals or proofs of a bulk phase.

## Observables and decision rule

Use the central half `[L/4,3L/4)` at every size and the central quarter
as a predeclared sensitivity check. For bonds require both endpoints
inside the window. Density is n=(1-Z)/2. Bond observables are the physical
interaction bond B (without onsite fields) and its transverse part K;
do not let a field redistribution masquerade as bond order.

Store one-point profiles, raw two-point matrices and connected matrices.
For overlapping Hermitian bonds use the symmetrized covariance
Re(<B_i B_j>)-<B_i><B_j>. For density this equals the ordinary connected
correlator. Separate S_connected(q)=sum exp[iq(i-j)] C_ij/M from
S_disconnected(q)=|sum exp[iqj](<O_j>-window_mean)|^2/M. Report each
divided by M as well, and their sum, over the whole q grid. Find the
dominant nonzero profile wavevector with q>=2pi/M; compare both separate
maxima and the common wavevector selected by the largest original chain.
Report central RMS/range, distance-averaged connected correlations, and
connected/disconnected separation-two hard-core-boson pair correlations
from three fixed adjacent anchors near L/4. Do not fit a pair exponent
from a handful of distances.

An ordering suggestion requires locally stable nonvanishing normalized
weight and long-distance correlations consistent across sizes, parities
and starts. A decreasing boundary profile together with decreasing
normalized connected weights is compatible with boundary oscillations,
but three sizes do not prove their infinite-size limit. Persistent
disconnected profiles alone do not establish order; vanishing one-point
profiles alone do not exclude a symmetry-restored ordered state.
Compare with the periodic and known-control results before choosing
`continue`, `reject`, or `missing test`. Any c assignment additionally
requires consistent scaling and independent evidence of positive-velocity
conformal modes. Otherwise report the useful negative or unresolved result.
