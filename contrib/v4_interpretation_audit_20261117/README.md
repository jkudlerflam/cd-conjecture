# V4: saved-state entropy, symmetry and momentum audit

## Claim and status

**[numerical] The asymptotic central charge remains undetermined.** The
microscopic chain is Hermitian, translation invariant in its bulk, nearest
neighbour, and has physical on-site dimension **d=2**:

\[
H=\sum_{j=1}^{L-1}\left[\frac9{16}X_jY_{j+1}
+\frac7{16}Y_jX_{j+1}+Z_jZ_{j+1}\right]
-\frac74\sum_{j=1}^LZ_j.
\]

Open chains retain the full endpoint fields. The target is the Virasoro c
of a unitary compact infrared CFT, conditional on that limit existing.
Here **c_ent means a finite-size entropy slope**, not c_eff. No finite-member
c>2 certificate or arbitrary-c family follows from this work.

The latest preserved L=96 to 144 estimate has the exact decomposition

\[
2.5668852487=1.9950941061+0.4591286499+0.1126624926.
\]

The first term uses completed chi256 **odd-parity** states at both sizes;
the second replaces the excited odd L=96 state by the lower-energy even
state; the third replaces completed L=144 odd chi256 by its **interrupted**
chi320 continuation. This is saved-record arithmetic, not an additive
physical correction formula. Selecting the excited parity to obtain a
number near 2 would not establish the ground-state c=2.

No DMRG was restarted. Six bounded, read-only checkpoint measurements took
**27.49 seconds total**, one worker and one BLAS thread, maximum RSS 1.60 GB.
The preceding JSON archive is unchanged. Work begun November 17 was
finalized November 18 following the battery warning.

## Entropy, parity and accuracy

For open critical chains the leading term is
S(l)=c/6 log[(2L/pi) sin(pi l/L)]+constant. We average entropy and log chord
over identical cuts. The factor 2 cancels in size differences. Applying this
formula presupposes a CFT scaling regime; a finite fit does not verify it.
[Calabrese and Cardy, opened](https://arxiv.org/abs/hep-th/0405152).

W6 denotes the 13 cuts with |l-L/2|<=6. Quarter and half denote fixed central
fractions; audit.py specifies inclusive-cut conventions.

| Sequence | Sizes | W6, log L | W6, chord | Quarter, chord | Half, chord |
|---|---:|---:|---:|---:|---:|
| Lower energy, completed | 48 to 96 | 2.81843 | 2.72787 | 2.74594 | 2.77089 |
| Lower energy, completed | 96 to 144 | 2.45422 | 2.42911 | 2.46681 | 2.53351 |
| Lower energy, interrupted 144 | 96 to 144 | 2.56689 | 2.54062 | 2.57427 | 2.62605 |
| Fixed odd, completed | 96 to 144 | 1.99509 | 1.97468 | 1.92762 | 2.11197 |
| Fixed odd, interrupted 144 | 96 to 144 | 2.10776 | 2.08619 | 2.03509 | 2.20451 |
| Fixed even, completed | 96 to 144 | 3.12829 | 3.09628 | 3.06840 | 3.06519 |

The lowest computed energy has parity **odd, even, odd** at 48, 96, 144.
Odd-minus-even W6 entropy differences are -0.020094, +0.031027, -0.045552,
despite energy splittings of only -3.886e-4, +1.917e-4, -1.054e-4.
L times the absolute splitting is 0.01865, 0.01840, 0.01518. Three small,
shell-sensitive gaps cannot discriminate a slow mode from a small eventual gap.

Central-half density changes by -0.003967 in the completed ground-state
96-to-144 comparison, versus -0.000481 at fixed odd parity. The near-2
branch also matches density better. But microscopic N is not conserved;
no measured dS/dn justifies subtracting a density correction. Parity, density
and entropy changes are correlated, not independent statistical errors.

Completed L=144 chi256 states failed the registered inter-chi entropy gate:
changes from chi192 were about 0.0184--0.0189. Four saved odd chi320 sweeps
increase W6 entropy another 0.00761, shifting c_ent by +0.11266. The last
two sweeps differ in W6 entropy by only 7.69e-6, but the final energy change
1.15e-8 exceeds the 1e-8 gate (the strict stop was 2e-9). Small last-sweep
drift is not a chi-convergence bound. The independent L=96 restart agrees
with its baseline, validating that accuracy rather than infinite chi.
Full histories, every window, and matched-parity comparisons are in audit.json.

## Boundary oscillations: strong evidence, no validated subtraction

After removing a fitted constant and log chord from each profile, we compare
entropy with cut density, physical bond energy and transverse bond energy.
The motivation is the observed entropy/bond alternation in open XXZ chains,
not a universal identity assumed for V4.
[Laflorencie et al., opened](https://arxiv.org/abs/cond-mat/0512475).

| State | Entropy residual q/pi | Density correlation | Bond correlation | Transverse-bond correlation |
|---|---:|---:|---:|---:|
| 96 even, chi256 | 0.88315 | +0.9594 | -0.9624 | -0.9855 |
| 96 odd, chi256 | 0.88291 | +0.8836 | -0.8831 | -0.8660 |
| 144 even, chi256 | 0.88151 | +0.8355 | -0.8324 | -0.7872 |
| 144 odd, chi256 | 0.88151 | +0.9802 | -0.9816 | -0.9929 |
| 144 odd, interrupted chi320 | 0.88151 | +0.9803 | unmeasured | unmeasured |

These are central-half residual correlations. At 144 odd chi256 the
transverse-bond regressor reduces residual RMS from 0.02232 to 0.00265,
but changes the single-profile chord coefficient from 2.3840 to 2.4070,
not to 2. At 96 even the coefficient changes from 2.7041 to 2.7582.
Oscillations explain much cut-to-cut scatter, without fixing the remaining
smooth slope. Middle-three-quarter fits and left/right holdouts are also
saved; reflected halves are not independent samples or uncertainty estimates.

Exact fixed-parameter XX controls give W6 log-L slopes 1.0061, 0.9950 for
one leg (c=1), and 1.9779, 1.9540 for two interleaved legs (total c=2), for
48-to-96 and 96-to-144 respectively. Two-leg midpoint-only slopes are
0.9004 and 2.9389. The same residual diagnostics count all leg bonds crossing
the cut in the kinetic regressor. These controls calibrate the procedure;
their errors are not subtracted from the interacting model.

## Momenta and mode counting

One-body matrices use the pre-common-rotation physical frame and a central-
quarter window M=L/4. We define Gij=<b_i dagger b_j> and
n(k)=sum exp[ik(i-j)]Gij/M, whose Brillouin-zone average is window density.
Both bosonic and Jordan-Wigner fermionic matrices are saved. Bosonic peaks:

| State | k_A/pi | k_B/pi | Folded separation/pi | Density-profile peak/pi |
|---|---:|---:|---:|---:|
| 48 odd, chi192 | -0.04663 | -0.95337 | 0.90674 | 0.88721 |
| 96 even, chi256 | -0.05762 | -0.94238 | 0.88477 | 0.88428 |
| 144 odd, chi256 | -0.05884 | -0.94116 | 0.88232 | 0.88013 |
| 144 odd, interrupted chi320 | -0.05884 | -0.94116 | 0.88232 | not remeasured |

Window resolution 2/M in units of pi is 0.167, 0.083, 0.056; dense Fourier
sampling does not improve it. Density and entropy oscillations are
consistent with interference between these **same two valleys**, also in
connected density correlations. No extra mode is needed for this wavevector.
Two peaks and equal weights do not prove two independently gapless sectors.

Raw n near 0.292 would give 2n near 0.584 for the fundamental harmonic of
an ordinary conserved-charge-one liquid, unlike the observed 0.880. However,
the exact U(1) assumption fails here: a conserved-charge density-momentum
theorem cannot be applied directly to raw n.
[Yamanaka, Oshikawa and Affleck, opened](https://arxiv.org/abs/cond-mat/9701141).
Higher-harmonic coincidences, including 3n near 0.876 for a paired description,
are not independent mode-count evidence.

[PHYSICS.md](PHYSICS.md) details the exact antiunitary symmetry, periodic-ring
momentum partners, full leading projected interactions, and number-breaking
vertices. The effective model includes a blockade and strong attractions;
it is not a free two-valley band, and its V4 truncation lacks an error bound.

**Conclusion:** c=2 is a plausible two-component description on these scales.
A flow to c=1 through locking one sector, or ultimately to a gap, remains
possible. An extra gapless sector requires evidence beyond c_ent>2, static
peaks and symmetry partners. No measured velocities, converged long-distance
exponents or independent low-energy towers select among these alternatives.

## Cheapest continuation when power is available

The single cheapest useful continuation is **at most two additional sweeps
of the saved L=144 odd chi320 state**, unchanged settings, one worker.
Recorded sweep times 261--312 seconds imply about **9--11 minutes**; cap
wall time at **12 minutes**, preserve the last complete checkpoint, and never
increase chi or L automatically. Stop early only after two consecutive
energy changes below 2e-9 and W6 changes below 1e-5; retain truncation and
residual diagnostics. Otherwise report stationarity unresolved. Passing
settles remaining same-chi drift at the 1e-4 entropy scale; it does **not**
pass the inter-chi gate or identify c. This continuation was not executed.

There is **no demonstrated cheap decisive continuation for c=1 versus c=2
versus an extra sector**. Accuracy-controlled larger scales and independent
low-energy mode evidence are required. Without measured slow velocities or
crossover lengths, a guaranteed decisive size or runtime would be invented.

## Checks, attribution and reproduction

- Author and analyst: OpenAI Codex, GPT-6 family. No independent human or
  second-model review is claimed.
- The unchanged C190 contraction code credits Claude Opus 5.5. Production
  states and checkpoints were reused; no solver changes or new eigenstates.
- onebody.py exactly reproduces archived L=96 even density and both matrices
  before allowing further measurements. Matrices pass Hermiticity and
  positivity; fermion eigenvalues lie in [0,1].
- algebra.py checks parity, antiunitary and rotation identities on eight
  sites and the full second-order projected Hamiltonian by virtual-state
  elimination. Errors are below 2.3e-14; runtime 0.007 seconds. This is not
  an infrared proof.
- audit.py checks the decomposition, Fourier normalization and unchanged
  archive hashes; runtime 0.52 seconds. Outputs: audit.json, algebra.json,
  measurement_*.json and source_manifest.json.

From this folder with numpy, scipy and threadpoolctl installed:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B audit.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B algebra.py
```

Committed onebody/*.json files suffice for these reductions. Optional
`python -B onebody.py` requires ignored local checkpoints named in the
records; its total measurement budget is 180 seconds, per-state alarm 60
seconds. No checkpoints are committed. Used numpy 2.5.5 and scipy 1.18.1.
AUDIT_PLAN.md records post-hoc scope, not a blind hypothesis test.

## Relation to existing work

This extends [the stopped bulk diagnosis](../v4_bulk_diagnosis_20261117/README.md)
and [the candidate work](../simple_qubit_candidates_20261117/README.md).
Recommended shared-status wording: **V4 remains an unresolved d=2 candidate;
latest interrupted ground-state slope 2.5669, large parity and accuracy
sensitivity, no c>2 certificate.** Top-level edits belong to the concurrent
coordination task and are not made here.
