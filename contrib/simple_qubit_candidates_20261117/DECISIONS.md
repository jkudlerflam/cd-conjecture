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
