# C28 progress

[proved, scope] Work is confined to work_C28. The project STATE and the requested C21, C22, C25, C26 and assembly sources have been read. No existing C28 progress file was present. Source results are inherited only where explicitly identified. No network or agents are used.

[proved, field obstruction] A fixed endpoint two-site diagonal readout on alphabet d has a real operator space of dimension d^2. Compression is linear, so it cannot provide arbitrarily many independently adjustable logical Z fields. C8's colored readout increases its port and color alphabet with the number of readouts. It cannot be imported unchanged as a fixed-alphabet solution. This does not prohibit a prescribed nonresonant field sequence.

[proved, candidate scope] C21 scalar Peierls phases on its cursor path are gauge removable even with its operator-valued on-site potentials. They leave endpoint return spectral measures unchanged. They cannot break its reversal restriction. A loop or additional noncommuting transport changes this hypothesis.

[ongoing] Examine C22 with chemical potentials before declaring its zero-field gap an obstruction. Its exact U(1) forbids pairing, but quartic Cooper transfer remains allowed. Compute the actual order-g^4 interaction and its weak-coupling continuum projection for an n=2 doped reference. Test complete two/three-cell legal-sector Hamiltonians at n=1,2 with local fields. No successful all-n construction or all-design no-go is claimed.

## New positive core and doped reference

[proved] `ff_bus.py` replaces C22's vacuum adjacency matrix by the uniform path Laplacian. The literal 17-letter positive repeated bond plus positive endpoint terms has exactly the uniform vacuum history times n unconstrained data bits as its kernel. Excited rails are unchanged and have diagonal d=4. Its source visibility is 1/sqrt(n+1). Endpoint density terms preserve the code exactly. They are retained order-lambda^2 perturbations, not positive terms required to annihilate the code.

[numerical, n=1,2] Literal physical dimensions 17^2=289 and 17^3=4913 agree with complete legal restrictions of dimensions 16 and 48 with zero stored-arithmetic residual. The local positive bond has minimum numerical eigenvalue -3.5e-15, consistent with its analytic positivity proof.

[proved] Unequal positive/negative conversion strengths only multiply the order-g^2 hopping by their product. They do not break C22's common reflection center. The earlier exploratory idea failed.

[numerical, n=2] Choose rail phases (pi/2+0.005,0,0), link phase pi/2 and field mu=0.00032039083881467 in leading g^2 units. The four roots are +/-0.673661697190 and +/-2.416591538066. The two positive velocities are 2.46936224408e-6 and 2.62593151981e-6. Exact U(1) forbids anomalous pairing at every order.

[heuristic, one-loop continuum input] The generated order-g^4 pair amplitude is -6.45074367714e-7+8.04111022848e-8 i in the recorded spinor gauge. The dimensionless BKT coordinates satisfy x/g^2=-0.09435238 and |y|/g^2=0.08124030, so x<-|y|. Pair tunnelling is irrelevant at this order. Other quartic resonances must be detuned before using this flow; the accidental backscattering resonance of the pi/2 reference has mismatch 0.10267884 after the flux change.

[ongoing] Full legal-sector ED with fields is being completed at all four requested (n,cells) pairs. The open-chain lambda->-lambda symmetry forces a fourth-order, not a cubic, leading remainder. The task's all-n independent field construction and all-n stability remain open. A finite number of diagonal spatial readout profiles is not a proof that no differently encoded cell can work.

## Completed handoff

[proved] RESULT.md now contains the complete positive 17-letter bond and endpoint table, all-sector code proof, gap bound, exact total U(1), density intertwining, the rank-four diagonal-field obstruction and the two restricted phase-repair obstructions. The general fixed-code linear-control dimension bound is stated separately from the all-n route-A problem.

[proved] `certificate.py` gives exact rational scalar resolvent weights at n=2 and certifies four simple leading Fermi crossings at the actual flux shift 0.005. This goes beyond a momentum-grid scan. It does not certify the continuum interpretation.

[numerical, n=1,2; two and three cells] All eight full legal-sector ED jobs at lambda=.08,.04 completed. Maximum dimension 110592, maximum individual runtime about 120 seconds. Full-error halving exponents are 3.99935 to 4.00114. Maximum eigenpair residual 7.14e-13. An initial literal n=1 ARPACK spectrum missed the known zero space and was discarded; dense diagonalization corrected that check. Assertions now verify every expected zero level. No failed output is used in RESULT.

[heuristic] The actual n=2 leading one-loop matching is on the pair-irrelevant side, with a strict margin 0.01311209 g^2. The other leading distinct-pair mismatches are nonzero, minimum 0.10267884. This is a small-coupling continuum inference, not a measured c=2 and not a finite-lambda phase certificate for the ED parameters.

[proved, files] `logical_chain_n2.npz` is the direct dimension-four numerics handoff. It contains the leading bond, order-g^4 bond and a full-in-g second-order bond at g=.02. `final_algebra.json`, `checks_results.json`, `certificate_results.json`, `handoff_results.json` and the eight ED JSON records contain the final checks. `algebra_results.json` is the earlier unsuccessful imbalance exploration.

[assumed/open, final scope] This run is partial. The core exists at every n, but independent diagonal fields stop at four controls and n conducting pairs plus their many-pair stability are not established. No all-design no-go was proved. The architecture counts 68 frozen and 1680 mobile are conditional counts, not demonstrated d0 values for unbounded c. The project STATE remains untouched because outside writes were forbidden. Local STATE.md records the same handoff.
