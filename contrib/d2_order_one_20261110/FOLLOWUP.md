# Adaptive follow-up registration, after Stage A

2026-11-10, before any Stage B physics output. Stage A is complete: 145 points
at L=10,12. XX and Ising calibrate but the c=2 two-leg control reads 3.07 at
this incompatible odd/even leg size pair. The r=40 reference is missed entirely
(density .0022, c_size .053). No Stage A point can pass a c>1 criterion.

Selected from the frozen grid:

| Label | Purpose |
|---|---|
| bridge_r2_u1.5_k4 | Preferred-ratio physical seed, even though its small-ring size slope is negative |
| bridge_r1.5_u2_k1 | Preferred-ratio entropy flag; explicitly momentum-degenerate at L10 |
| bridge_r1.5_u1.5_k3 | Preferred-ratio entropy flag with uniform momentum eigenstates at both initial sizes |
| bridge_r4_u1_k4 | Wider-ratio continuation bridge, R=28, density near .16 |
| attractive_Jz-0.5 | Best Stage A rational attractive point, R=2 |
| reference_peotta_exploratory | Existing unresolved attractive-junction reference, R=7.77 |

Larger-ring follow-up: L=14,16, optionally18 for these points; controls L=16
and20 (and existing12), so the two-XX legs have compatible even lengths.
Use at least4 low eigenvectors and record degeneracy/translation diagnostics.
If a low multiplet is found, diagonalize translation within it and report raw
and momentum-selected entropies separately. The original Stage A results
remain untouched. Single-vector Lanczos with k=2 can miss degenerate partners;
its returned gap is not a certification of nondegeneracy.

Open-chain follow-up: L=24,48, both parity sectors for the three preferred
bridge points, r4 bridge, attractive rational point and Peotta reference,
as local runtime permits. Controls XX and twoXX at L=24,48 in their known
even ground parity. Start chi=64 at24 and96 at48; extend a surviving excess
to L=96, chi128, and compare chi128/192 at48 or96. Register any reduced
selection before larger-size results. The negative small-ring slope at the
r2 seed must not veto its OBC check.

Wrapper gate: C190 solver matches independent bit-assembled ED at L8 for
symmetric and attractive bonds, XX and twoXX; max energy error5e-14 and
entropy error8e-15 for nondegenerate sectors. Degenerate XX odd-sector
states are compared by energy only; arbitrary ground vectors have different
entropies. A point-level distance plumbing bug in the draft wrapper was fixed
before any physics/control OBC run. Checkpoints are written every sweep.

## Control continuation, registered after the first OBC controls

The twoXX control completed its schedule but did not meet the energy-step
tolerance at chi64/96. Continue its L24 and48 states to chi192. Exact
correlation-matrix calculations now give c_fractional(24,48)=1.71271 and
c_central3=2.05876; the fractional window has a sizeable finite-size effect
even without truncation. At48,96 the exact predictions are1.93810 and2.03332.
Use these actual finite-size control values rather than forcing them to2.
The control's original stage histories remain in its continued record.

## Third size and chi checks, registered after initial bridge output

At24->48 the r4,u1,hz=-7/4 bridge reads c_central3=2.763 and
c_fractional=3.008. Its density drifts .3344->.3013 and its chi32->96
entropy step is .0875, so it is not a converged c>1 result. Continue both
L48 sectors to chi192 and run both L96 sectors atchi192. Add an independent
L48 seed23 atchi192 if the excess persists after the continuation.

The preferred r2,u1.5,hz=-1/2 seed reads .198 and1.002 in the two windows;
run both L96 sectors atchi128 to test whether this disagreement settles.
This is a negative-control continuation of the physical route, not promotion
by its entropy. These selections precede any L96 output.

Retain the Peotta reference for an L96 check atchi128 in both parities,
regardless of its24->48 reading, because its existing48->192 excess was
the independent rationale for this route. Do not run the other preferred
flags beyond48 in this initial bounded investigation: their larger-ring
slopes are already belowone and their OBC slopes are at most.78.

## Two nearer bridge points, after the r4 chi continuation

The r4 L48 excess survives chi96->192 (entropy step .00181 in the odd ground
sector), before the L96 result is available. Add exactly two nearer bridge
points atL24,48 in both parities: bridge_r3_u1_k4 and bridge_r3_u1.5_k4,
with the same chi64/96 schedule. Their coefficient spreads are14 and28/3;
the second enters the exploratory band. Both fix hz=-7/6 and Jz=3/4,
with D=1/12 or1/8. Stage A was already seen: the first had density~.22 and
c_size1.23; the second had density~.41 and a large degenerate-shell artifact.
The decision is motivated by continuing the r4 anomaly, not interpreting
that artifact as a CFT. This closes the initial local exchange grid; no
additional unregistered broad sampling is planned for this milestone.

## One same-chi convergence extension

The r2 seed's odd L96 sector exhausted six final sweeps atchi128 and fails
the sweep test; its canonical entropy also disagrees with the sweep profile.
The even sector is lower and converged, but ground-sector comparisons should
not silently treat the odd state as settled. Continue the odd checkpoint for
up to six further final sweeps at the samechi128. This is a convergence fix,
not a change of Hamiltonian or a new high-entropy selection.

## Density anchors, before finishing the initial milestone

The selected comparable points have raw densities .36-.45, so they are not
a density-matched continuation of the original dilute gas. To measure this
limitation directly, add two already-screened lower-field points atL24,48
in both parities, chi64/96: bridge_r2_u1.5_k1 (hz=-7/8,R=14/3,
ring density~.195) and bridge_r3_u1.5_k3 (hz=-5/4,R=10,
ring density~.173). These are fixed coarse density anchors, not an optimized
constant-density path; raw density includes virtual pairs. They are selected
for that physical comparison despite low ring entropy slopes. No new
exchange grid is added. This registration precedes their OBC output.
