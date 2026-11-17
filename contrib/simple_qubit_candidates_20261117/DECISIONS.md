# Sequential choices and conventions

1. Source reuse selected V=8 and V=12 before new diagnostics. C190's
   actual symmetric production field is -19.976047250613814, corresponding
   to a=0.9581099754474565, not a=0.95 or the rounded STATUS field.
   C177's old low-V values are a=0.9669800331,kappa=0.9285827243 at V=8
   and a=0.9610913978,kappa=0.9658725946 at V=12. Those old datasets do
   not test exactly kappa=1. Their source hashes and raw small outputs
   are in `source_audit.json`.
2. The initial a=1,L=32 pilots had middle densities about 0.35--0.43,
   depending on V and parity, far above the older candidate's density.
   The next rational choices were a=31/32 at V=8 and a=24/25 at V=12.
   The V=8,31/32 pilot remained at 0.23--0.26. A further bounded move to
   a=24/25 yielded 0.196--0.228 there. We therefore select the common
   a=24/25,kappa=1 for the first size comparisons at both V, and retain
   the other pilots as density/tolerance evidence, not as discarded data.
3. Validation at L=8 in both parities compares the literal Pauli matrix,
   full endpoint fields, energy, every entropy cut, site density and a
   pair observable. All pass. The largest missing-half-endpoint operator
   norm is exactly |hz|; omitting endpoint halves changes the model.
4. The fixed-work L=64,chi=128 benchmark measured 0.047 GB of persistent
   arrays and peak RSS 0.231 GB; three central 100-matvec updates take
   about 0.13 seconds each. This uses the audited C195 numpy kernel;
   the actual campaign uses the corresponding C190 CPU solver and logs
   every complete sweep, whose measured times decide extensions.
5. Counts use the physical onsite field once, not two half-fields in a
   bond representation. There are four independent real coefficients
   (Delta,D,Jz,hz), three after fixing Delta=1. In the original common
   X,Y,Z basis the six monomials are XX,YY,XY,YX,ZZ,Z. The same onsite
   U=exp(-i pi Z/8) at every site sends the Hamiltonian exactly to
   (1/2+D)XY+(1/2-D)YX+(V/4)ZZ+hz Z: four monomials. No identity or
   telescoping term is counted. The two bases are held fixed across all
   points; `simplicity.py` verifies this conjugation.
6. Both spreads must be visible. At V=8,a=24/25,kappa=1 they are 124.16
   in the original basis and 8.27733 in the rotated basis. The small
   difference between the rotated exchanges is still 2D=1/16, equal
   to the second-order t2 here. The rotation shortens the expression,
   but does not remove the small dynamical scale or field cancellation.
7. At chi=128 the ground branches of both V=8 and 12 pass the four
   local numerical gates on L=32,64. Their W=6 size readings are 2.421
   and 2.369, while the exact two-leg control reads 2.389. That agreement
   justifies a size extension; it does not establish a central charge.
   W=1 and W=12 show substantial finite-size effects, retained in full.
   The odd V=12 branch has a 0.011 entropy change from chi=96 to 128,
   and the L=64 two-leg control still fails the energy-stationarity gate
   despite an entropy error below 0.0005 against its exact value. These
   specifically justify chi=192 continuations and matched L=128 controls.
   The L=128,chi=192 fixed-work benchmark uses peak RSS 0.673 GB and
   0.42--0.51 seconds per central 100-matvec update. The first size
   extension is capped at L=128,chi=192, with any further chi increase
   requiring the measured convergence results.
8. The sibling `qubit_minimal_terms_20261117` independently obtained the
   same four-monomial rotation and invariant hierarchy warning. It does
   not prove the necessity of the DM term or field in every possible
   c>1 model. Its zero-field two-monomial residual tests have no new
   c>1 evidence and do not displace the supported blockade family here.
9. Reading that sibling's sources revealed the earlier
   `d2_order_one_20261110` search. Its point D is exactly our family at
   V=4,kappa=1,a=1: original spread 28, common-rotation spread 4.
   It already has L=24,48,96 full-chain data. The last size gives slopes
   2.51--3.00 across parity/window choices, a central density range about
   0.302, and no completed chi comparison there. The old checkpoints
   are unavailable locally. This is useful low-hierarchy evidence but
   not a more reliable CFT identification than the V=40 baseline.
   Preserve it in the ranked comparison and, if resources permit,
   perform a bounded new rational-field check rather than a new grid.
10. For the legacy V=4,a=1 point, the useful missing check is chi
    convergence at L=96. A fresh two-parity calculation through
    chi=128,192,256 directly addresses it despite the absent old
    checkpoints. It is more informative than rerunning tiny-ring
    entropy slopes. Keep this as a third, explicitly exploratory
    comparison; its large old slope does not make it the preferred CFT
    candidate. Its global variance will also be measured. No claim of
    a phase boundary follows from a handful of such points.
11. At V=8,L=128 both chi=192 runs pass residual, sweep-stationarity
    and discarded-weight gates, but fail chi convergence: Sbar6 rises
    by 0.0447 (even) and 0.0417 (odd) from chi=128. Extend only this
    preferred lower-hierarchy point to chi=256 for the next check.
    The chi=256,L=128 fixed-work benchmark stayed below 0.86 GB RSS.
    Two production DMRG processes remain the concurrency limit; tiny
    L=8 validation jobs and fixed-work benchmarks have also run while
    production proceeds, without changing numerical tolerances.
12. The matched L=128 controls are complete. Exact one-/two-leg W=6
    readings at 64-to-128 are 1.004026 and 1.858274. The corresponding
    two-leg DMRG value is 1.831440, with final energy drift 1.22e-6 and
    entropy error -0.003144 at L=128; it stays flagged. We use its exact
    reference, not an assumed asymptotic value 2, and do not spend this
    bounded campaign repeatedly fitting the solvable control. V=12's
    chi=128-to-192 entropy increments 0.0474--0.0514 fail as V=8's do;
    its flagged W=6 size reading is 2.109011. It receives no further
    extension because it offers no simplicity advantage over V=8.
13. The final queue contains only two V=8,L=128 chi=256 continuations,
    the two legacy V=4,L=96 calculations, eight L=32 tolerance pilots
    (a=0.955,0.965 at kappa=1 and kappa=15/16,17/16 at a=0.96), and one
    independent-seed V=8,L=64 even run. Two production workers run these
    sequentially. A further failure at chi=256 will be reported as
    unresolved rather than triggering an unbounded resource escalation.
14. The C190-style local diagnostic gates use delta E<=1e-8,
    delta S<=1e-3,delta n<=5e-4,residual<=1e-8, discarded weight<=1e-7,
    and last chi-step delta Sbar6<=0.01. They are looser than the
    registered within-stage stopping rule (2e-9,5e-4,2e-4,1e-8).
    The analysis explicitly reports both, with `schedule_limit_hit`
    preserving every exhausted final stage. A passed local gate does
    not erase a stopping-rule failure or certify the infrared limit.
15. The raw C190 valley imbalance includes both +/-pi/2 grid bins in
    one half because floating-point cosine is slightly positive there.
    Retain that value and also split those two bins equally; corrected
    imbalance is at numerical roundoff on the symmetric candidates.
    Neither convention tests charge locking or proves two gapless modes.
16. A fresh-directory L=32 one-leg XX reproduction test exactly matches
    the archived energy, density and all three entropy-window values.
    This tiny control validation overlapped the two production workers;
    it is not an additional parameter scan. `reproduce.py` reconstructs
    warm-start dependencies sequentially in a separate directory, so
    missing ignored checkpoints do not prevent fresh reproduction or
    silently reuse the committed JSON records as new calculations.
17. V=8,L=128,chi=256 completes in both parities. The strict stopping
    rule is met after six sweeps, but chi increments 0.013937/0.013014
    still fail the 0.01 gate. Its 64-to-128 lower-energy-branch W=6
    estimate moves from 2.059737 at chi=192 to 2.172839 at chi=256;
    all four parity pairings span 2.172839--2.480305. Retain the flags
    and stop this extension at its stated cap. The two continuations
    took about 16.4--17.0 minutes each and stayed below 1.95 GB RSS.
    They are followed by the preselected V=4 and tolerance/start tests,
    with no further V=8 size or chi escalation in this campaign.
18. The fresh V=4,L=96,chi=256 runs pass all four local gates and the
    strict stopping rule in both parities. Entropy increments from 192
    are 0.005120/0.005618 and variance per site is 1.52e-7/1.68e-7.
    This changes the simplicity-first follow-up ranking: V=4 precedes
    V=8, with V=12 third. It does not establish a CFT: density range
    remains 0.30147 in the lower-energy state, profile/window estimates
    disagree, and only the old L=96 chi gap was addressed here.
19. All eight V=8 tolerance pilots pass local gates. Their lower-energy
    densities span 0.15540--0.26437 and sectors switch. The physical
    field steps are +/-1/1600 and DM steps +/-1/512. Do not call this
    a robust critical window or a thermodynamic boundary. The L=64
    independent even seed agrees to 3.6e-12 in energy and 2.1e-8 in Sbar6.
20. All 46 scheduled runs are complete; no production job remains.
    Eight earlier unused config alternatives are retained and listed
    as unrun in the inventory, not counted as numerical results. The
    final rank concerns next tests, not an established c>1 phase.
