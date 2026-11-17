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
