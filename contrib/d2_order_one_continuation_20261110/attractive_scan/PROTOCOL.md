# Attractive-junction continuation: registration

OpenAI Codex, GPT-6 family, 2026-11-10. Written before new output in this
directory. The user asks to continue until a promising comparable-coupling
candidate is found. This is one independent route, not a counterexample claim.

The full microscopic model is the uniform nearest-neighbour qubit chain
Jx XX + Jy YY + D(XY-YX) + Jz ZZ, with zero on-site field. Physical sites
are never blocked or projected. Both exact total-parity sectors are optimized.
The unchanged C190 solver is called through the earlier independently
validated wrapper. Earlier contribution files and checkpoints remain read-only.

The existing Peotta reference has
Jx=(cos(1)+0.7)/2, Jy=(cos(1)-0.7)/2, D=sin(1)/2, Jz=-3/5.
Its previous C190 central-three/fractional slopes are 0.954/1.007 at 24->48
and 1.101/1.133 at 48->96. The older C184 work reads about 1.31 at 96->192,
near a first-order boundary. This is a reason for an independent extension,
not evidence that the infrared central charge is above one.

## First comparisons

1. `attr_ref`: extend the reference to L=192 at chi=192, both parities.
   The existing L48/96 raw records supply the smaller sizes without being
   overwritten. The new chain starts independently with seed17.
2. `attr_rat`: nearest simple rational comparison
   (Jx,Jy,D,Jz)=(5/8,-1/12,5/12,-3/5). Run L=48,96 at chi128, both
   parities. If its finite-size behavior remains useful, register L192
   before that output. Its literal coefficient spread is 7.5, in the
   exploratory band, not the preferred band <=3. Both models have five
   expanded Pauli strings, four grouped terms, and exactly absent field.

Use chi stages16,32,64,128 for L<=96, adding192 at L192; every stage runs
three to six zero-mixer sweeps to the existing stopping conditions.
Final local residual target is1e-10; report the registered physics gate1e-8,
energy step1e-8 and entropy step1e-3, truncation, final canonical/sweep
agreement, density profiles, and XX/YY/ZZ correlations through L/2.
Save every sweep boundary in /private/tmp/d2-order-one-attractive-checkpoints.
At most two single-BLAS-thread CPU processes run here; no cluster submission.

Entropy estimators and finite-size controls are unchanged from the initial
search: central three cuts and trapezoidal [3L/8,5L/8], c_size=6 deltaS/log2.
Exact XX and two-XX controls are available through L192. The two-XX model
is range two and is only an estimator control. Profile fits, parity-sector
comparisons and window dependence must accompany every entropy slope.

No entropy threshold excludes c only slightly above one. A promising point
must retain excess under increased size and chi, have a credible profile,
and not owe it to an ordered cat, phase separation, a boundary crossing or
unconverged optimization. A full CFT claim needs independent infrared tests.

Record raw and invariant interaction scales, coefficient spread, term count
and tuning resolution. The rational comparison has exact fixed coefficients;
no fine tuning is claimed. Any later field/exchange grid will be registered
before output and separated from this initial comparison.

## Preferred-ratio lead, registered before new output

The raw C184 summary (not just its REPORT headline) contains an unextended
phi=1.2, Delta=.8, V=-2, hz=0 point with old four-cut c(48,96)=1.429.
It has not been tested at192 in that grid. Its Pauli spread is about2.656,
while exchange singular values are .9,.5,.1 (spread9): the transverse
cancellation t-Delta=.2 is moderate and must be disclosed.

3. `pref_ref`: Jx=(cos(1.2)+.8)/2, Jy=(cos(1.2)-.8)/2,
   Dz=sin(1.2)/2, Jz=-1/2,hz=0. Run48 and96, both parities,chi128.
4. `pref_rat`: exact coefficients (Jx,Jy,Dz,Jz)=(7/12,-5/24,11/24,-1/2),
   hz=0. Equivalently integer ratios(14,-5,11,-12)/24; literal spread2.8.
   Run48 and96, both parities,chi128.

Both are five expanded Pauli strings/four groups. Neither has parameter
precision below1/24 in the rational comparison. A positive size slope
excess in both registered windows, in lower-energy parity states without
convergence failure, triggers192/chi192 for that point and chi192 at96.
The trigger is a screen, not a conformality conclusion; ordered-cat and
first-order diagnostics remain required. The original running reference192
jobs retain checkpoints and finish. Unstarted `attr_rat` work is deferred.
Parent authorized a temporary maximum of three single-BLAS-thread workers
while those two old192 jobs finish; otherwise maximum two. No work above
L192 or chi256 is authorized here without parent confirmation.

Scheduling refinement after preferred reference48/96 completed: parent asked
that pref_ref192 start as legacy192 slots free, without waiting for rational96.
A priority driver starts P0 after one legacy completion and P1 after both.
New launches use per-job locks and two global CPU locks; thus no duplicate
optimization and at most two numerical workers after legacy workers finish.
The pre-existing legacy workers plus these two slots obey the authorized
brief maximum of three. Rational192 extensions wait for pref_ref192 both
parities; this changes scheduling, not the registered comparison or gate.

## Second preferred-ratio reserve, registered before output

Parent's independent raw audit found phi1.2,Delta.8,V=-2.8 with old
fractional c(48,96)=1.335, profile values1.22--1.369 and raw spread3.199.
The following simple rational near-neighbor moves into the preferred band:

5. `pref_zrat`: H=sum[7XX-3YY+6(XY-YX)-9ZZ]/12, hz=0.
   Pauli spread3 exactly; five strings/four groups. Fixed integer ratios,
   four distinct nonzero coefficient magnitudes; no decimal fine tuning.
   Fermion t=sqrt(10)/3, phi=atan(3), Delta/t=sqrt(10)/4,
   V/t=-9/sqrt(10). Run48/96, both parity sectors, chi128, after current
   preferred48/96 jobs finish. This reserve takes priority over an expanded
   Jz offset grid; additional larger sizes require a recorded promotion.

Its motivation is another entropy flag near the other attractive boundary,
not evidence against a ferromagnetic-cat or first-order explanation. The
same entropy, convergence, parity, density, and long-correlation checks apply.
One sequential worker shares the global two-slot lock; pre-existing legacy
workers remain subject to the temporary three-worker allowance.

## Explicit deferral after new reference-size evidence

Parent directed stopping the unstarted preferred-rational192 pair and the
subsequent96 chi128-to192 refinements after the preferred reference192
preliminary result dropped to c_fractional about1.14 and its parity gap
fell rapidly, with persistent order/cat warnings. The parent concurrently
has a different R8.75 candidate retaining about1.4 through192.

The extension orchestration process and both pref_rat192 children were
interrupted while waiting for the completed pref_ref192 gate. Their logs
show KeyboardInterrupt at that wait; no pref_rat192 record/checkpoint or
new variational sweep exists. The original48/96 outputs remain intact.
The two already-active pref_ref192 parity calculations finish normally;
pref_zrat96/P1 finishes to complete its already measured48/96 pair.
The later preferred96 chi192 refinements and any further reserve promotion
are deferred, not silently counted as null results or completed checks.
`defer_preferred_extensions` prevents accidental relaunch. No alternate
boundary-pinning solver was started.
