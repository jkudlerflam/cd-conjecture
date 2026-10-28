# C44 progress

2026-10-28. Active audit, no final verdict yet.

Read the project STATE, STATUS, C39/C42 results, the C39 referee, C32 result and proofs, C25 result, and the four requested C8 notes. Also read the length_filter and arrow_wrapper dependencies. All writes stay in work_C44. The project STATE is outside the permitted write scope.

Established from the definitions so far. C39 supplies positive repeated bulk and positive endpoint terms, a complete initialized code, exact U(1), and vacuum versus bus superselection. Its biased vacuum does not obstruct the wrapper. The filter is an exact invariant ground factor in port-source calculations, even when its individual bond coefficients depend on hole positions. Endpoint compensation does not remove the adjacency change of the excited core. The full scalar second-order term must be kept.

Remaining work. Reconstruct the literal positive core terms locally, evaluate all four wrapped source bonds including helix excitations, diagonalize small complete legal cell sectors, check transport and physical ring gluing, and write explicit all-n error and exclusion budgets. A full 4464-state many-site tensor-product ED is not feasible; numerical restrictions will be stated explicitly and full-space exclusion will be proved analytically.

Numerics completed so far. All four length pairs were contracted with every charge and reachable helix excitation at n=1,2, two vacancy ratios and two core/charge-gap ratios. Scalar length mismatch is nonzero and much larger than the tiny diagnostic hopping. No relative payload-accuracy claim follows from those diagnostic points. Exact finite norm bounds, rather than these points, will close the all-n schedule.

The n=1 complete legal initialized core, charge and occupation-filter blocks have dimensions 1728 and 4320. They give the predicted two-dimensional zero code and positive next level. The first subspace check failed because degenerate ARPACK vectors were not orthogonal; QR of their span fixes the projector calculation. The original code residual was already below 4e-16. This was a numerical projector error, not a failed code test. Endpoint helix spectra, direct typed-hole transport and physical cyclic-gluing checks are in progress.

The finite audit is complete. PROOFS.md now gives the full repeated-bond composition, every sector margin, a nonasymptotic positive charge-gap bound, the exact four-length source formula, and explicit equalization and mixed-term budgets tied to J0 and C39's lacunary margin. The order u_mov*lambda vanishes by vacuum/bus conservation. The first potentially nonzero marker-port coefficient is lambda^2*u_mov; its resolvent formula and norm bound are retained. Readout leakage is budgeted separately. Odd periodic cell rings are not incorrectly assigned an even-in-lambda expansion.

Completed ED and transport checks. The projected physical-ring dimensions 104,432,1024 reproduce the joint momentum spectrum to 2.4e-14. Direct typed-hole transport checks use C39's actual biased code at n=1,2. Complete helix matrices at lengths 4,5,6 have unique kernels and the predicted gap bounds. Independent bare-source and direct charge-source inverse checks pass. The final writeup will distinguish these checks from full microscopic many-cell ED, which was not run.

Final result written to RESULT.md, detailed derivations to PROOFS.md, and handoff to STATE.md. The finite wrapper passes for every n. The thermodynamic elimination, auxiliary-sector stability and Gaussian infrared conclusion remain assumed. No full microscopic many-cell ED or central-charge claim was made. All computations completed in seconds. No work outside the folder was performed.
