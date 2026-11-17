# Microscopic qubit continuation, frozen 2026-11-16 21:02 EST

Author: OpenAI Codex (GPT-6 family). Written after auditing the old saved
records, before any new solver validation or production output in this
contribution. This file is immutable; amendments must be separate dated files.
Baseline repository: bbeb1403a00564f1316d42b0afc0549e9115db8e.

## Model and question

The bulk is a fixed finite Hermitian qubit nearest-neighbour bond with one-site
translation symmetry and size-independent couplings. In Pauli normalization,
Jx=1/2, Jy=-1/2, D=1/160, Jz=10. The continuation uses the **full-precision
production field hz=-19.976047250613814**, on every site including both ends.
The question's -19.976047 is its rounded display value, not numerically equal:
the difference is 2.50613814e-7 per site. No result for the rounded field will
be claimed without a separate run. The operator-norm bound L|delta hz| is
1.28314e-4 at L=512, much larger than the observed parity splitting, so that
bound cannot justify transferring a branch selection to the rounded field.

Use `work_C195/models.py:candidate('sym40',L,hz,gauge='stag')` unchanged.
The staggered onsite unitary is a computational gauge, not a change of physical
translation symmetry. Each onsite field is included once. Equivalently the
bulk bond includes hz(ZI+IZ)/2, and the finite open chain includes additional
hz Z/2 at each end. Omitting those endpoint halves changes this calculation.
Retain the original energy shift and reconstruct E_phys with E_off. There is
exact Z2 parity, **not exact microscopic particle-number U(1)**.

Here c means the total nonchiral Virasoro central charge of the complete
unitary conformal bulk theory in the repository's entropy convention (one XX
liquid has c=1; no extra factor for left/right movers). The entropy estimator
does not establish existence of that bulk CFT.

## Old evidence, not new output

The raw ground-branch values are c_two(128,256;6)=2.02533855 and
c_two(256,512;6)=1.98273733, using chi=512,1536,2048 respectively, odd parity
at all three sizes. At L=512 the step 1536->2048 in Sbar6 is 0.01037536,
outside C4=0.01. Only odd was continued to 2048; the latest even is 1536.
At common chi=1536, E_even-E_odd=1.5682581e-6. At unmatched highest chi it
is 2.0338057e-6. These are variational numerical comparisons, not rigorous
sector-energy bounds. The matched-sector check remains part of the new work.

## New stage and stop rules

Choose **chi=2560**, the next economical stage above 2048. Evaluate L=512
in both parities. Resume only independently verified checkpoints, read only:
odd from its 2048 state, even from its 1536 state through an intermediate
2048 stage. Run each new stage with alpha=0, residual target 1e-9,
maxit=300, mmax=32, minimum 3 and maximum 8 sweeps; stop after the minimum
when |dE|<4e-9 and |dSbar6|<5e-4 between consecutive sweeps. Reaching the
sweep cap is not automatically convergence. Do not loosen a failed criterion.
Use float64 in the real staggered gauge and the unchanged C195 offload path.
Separate output/checkpoint paths; never overwrite historical runs.

L=256 reference: the odd chi=1536 run, whose last chi increment is 0.00307639.
Its even competitor is at chi=1024 and lower accuracy. If compute exists,
continue even to 1536 with minimum 3/max 8 sweeps, tolE=2e-9 and the same
other settings, to verify its order at matched chi. Further chi=3072 is not
silently substituted: a dated amendment is required, motivated by resource
checks or a failed 2560 convergence test and identifying outputs already seen.

Before a full allocation, benchmark the real sym40 local update at chi=2560,
including canonicalization, environments, solve, split and transfers, with
bounded memory. The old p153 complex benchmark is only an estimate. Reject a
launch if available memory cannot cover tensors, environments, solver scratch,
checkpoint replacement and final measurements with 25% headroom. Use one run
at a time unless an existing allocation explicitly has capacity for both.
No changes to other jobs; no assumption of remote access or old checkpoints.

## Pass/fail/undecided, fixed before outputs

1. State provenance gate: validate L, parity, model, field, gauge, checkpoint
   sweep count and final bond dimensions. Recompute energy, all entropies and
   density against the saved run (tolerances 1e-8 total energy, 1e-7 entropy,
   1e-7 density), check finite tensors and exact parity blocks, and hash files.
   File names alone do not authenticate a checkpoint. CPU/GPU and resume
   validation must meet 1e-10 energy and 1e-8 entropy on small chains.
2. T1 for **both L=512 sectors**: residual <=1e-8, last-sweep |dE|<=2e-8,
   max over W=1,6,12 of |dSbar_W|<=1e-3, |dn_mid|<=5e-4, maximum discarded
   weight<=1e-7. Require the stop rule as well. Publish every final-stage
   sweep, not just the final energy. These residuals are local, not a global
   eigenstate certificate. Apply the analogous 1e-8 energy threshold at 256.
3. C4: |Sbar6(2560)-Sbar6(2048)|<=0.01 for both sectors, comparing the last
   sweeps of individually stationary stages. Report W=1,12 increments too.
4. Numerical ground-branch check: both sectors reach the same final chi and
   satisfy 1-3, have the same ordering at the last two matched chi stages,
   and their energy splitting exceeds ten times the sum of their absolute
   last-chi energy changes and last-sweep energy changes. This added guard is
   numerical, not a rigorous bound; failure leaves the branch undecided.
5. Same formulas: Sbar_W averages cuts |l-L/2|<=W, c_two=6 delta Sbar/ln2,
   n_mid averages sites [L/4,3L/4). Final canonical entropies decide c_two;
   last-sweep values, all four parity pairings, matched-chi values and c_prof
   are retained. Require |n_mid(512)-n_mid(256)|<=0.006 on the selected branch.
6. Controls: verify C190 one-XX and two-XX raw data against exact correlation
   matrices using these windows and the same solver family; run local controls
   with the actual C195 driver and verify interruption/resume. For a fully
   matched **new production verdict**, run/verify both controls at L=256,512
   with the C195 solver and final-chi stop/convergence tests, with mu at the
   midpoint of free levels for N nearest 0.1957L (balanced legs for xx2).
   Require each Sbar_W within 0.001 of its exact finite-size value and
   |c_two_DM RG-c_two_exact|<=0.03 at W=6,12. Report deviations of the exact
   finite-size estimator from asymptotic c=1 or 2 rather than tuning mu to c.
   C198's c2 data use a different solver and particle sector and cannot alone
   satisfy this new matched-solver production gate.
7. PASS only if 1-6 hold and c_two>=1.6 at W=6,12. FAIL (entropy criterion)
   only if 1-6 hold and c_two<=1.3 at both. Otherwise UNDECIDED, identifying
   unavailable compute separately from actual failed convergence checks.
   C4 failure explicitly blocks this new pass, unlike historical C195 T2.

An increasing chi trend is an observation, not a lower bound on entropy or c.
Raw values are primary. Optional linear-in-discarded-weight and power-in-chi
fits are labelled model-dependent extrapolations with their input points and
spread. No fitted value replaces a failed convergence gate. A pass supports
finite-size c>1 evidence; it does not prove an asymptotic counterexample.

## Correlations, if checkpoints and resources are available

Use C193 `meas193.py`, its checkpoint checks, `frozen.json` windows and
`ana193.py states` unchanged for both sectors at L=256,512. Publish eta_b,
pair exponent (d=3 primary, d=2 diagnostic), both K_s estimates, uncertainties
and 128->256 vs 256->512 drifts. Apply the original B2/B3 inequalities shifted
up one size doubling, retaining the rule B2 must hold for both estimators,
B3 if either triggers, and the same KT drift/error formulas. Test calibration
at these sizes first. Keep the original C193 B3 result and its known small-L
structure-factor failure visible. The boundary-corrected two-exponent fit is
secondary, labelled separately; do not change windows after inspecting outputs.
Saved C195 gwin=48 correlators do not reach the registered fractional windows
at L=512, and their CP is not the registered Dpair observable. Without the
checkpoint/full correlators this test is NOT EVALUATED, not gapless.

C196 and C196/b were read: neither infinite-chain route converged at the
candidate, including liquid starts. No new infinite-chain run is authorized
by this registration, and no exact U(1) constraint will be imposed.

## Local fallback fixed now

This host has 48 GiB RAM, no CUDA GPU or local C190/C195 checkpoint; the old
Typhon SSH connection timed out. If these conditions persist, complete the
saved-data arithmetic, exact finite-size controls, small-chain ED and
offload/resume gates, authenticated continuation/config generation and a
resource estimate. The new production verdict remains UNDECIDED. Local
small-chain tests must never be presented as chi=2560 production evidence.
