# C44 local state

2026-10-28. Complete.

[proved] The finite mobile-wrapper audit passes on C39's 46-letter initialized core for every n. The full alphabet is 4464. Two selected physical lengths n+3,n+4 have the same complete 2^n code. Full-space unperturbed exclusion, exact first-order maps, physical momentum gluing, exact second-order source formulas, finite equalization errors and mixed-term bounds are in PROOFS.md. RESULT.md is the concise verdict.

[proved] The wrapper is only approximately equalized at finite parameters. Scalar length-density and traceless length-payload terms are retained. One port insertion with any number of marker moves cannot return to the code. Lambda^2 times marker hopping and marker/readout terms can survive. Odd periodic cell rings can have winding odd-port terms. The parameter schedule includes J0, the lacunary momentum and velocity scales, the actual charge gap, compensator norms, readout leakage and scalar retuning.

[numerical] Checks cover all four source length pairs at n=1,2, complete legal core/charge/filter ED at n=1 in dimensions 1728,4320, endpoint helix ED in dimensions 16,32,64, direct code transport at n=1,2, and full projected physical-ring spectra in dimensions 104,432,1024. No full microscopic many-cell ED or central charge. Diagnostic parameters are not the all-n lacunary schedule. Small nonscalar source differences are not relatively certified.

[assumed] The remaining step is thermodynamic local elimination, massive-sector stability and Gaussian infrared stability at fixed nonzero couplings. The finite audit is not a proof of c=n+1. Finite-cluster perturbation theory is explicitly separated from this input.

[proved, file map] audit_numeric.py computes full source bonds and small legal-cell spectra. transport_ring.py checks typed-hole transport and cyclic gluing. bounds_checks.py checks selector inequalities and exact shortest-path coefficients. consistency_checks.py checks literal local terms and independent resolvents. The local C39/C32 snapshots are core_model.py, core_literal.py and wrapper_tools.py. Their provenance and modifications are documented in RESULT.md.

[proved, write scope] No files outside work_C44 were modified. This local STATE substitutes for the project STATE update forbidden by the task's write restriction.
