# C58 local state

2026-10-28. Completed thermodynamic transfer.

[proved conditional on imported finite and massive results] PROOF.md closes estimate (c) for the original C44 wrapped C39 bond at d0=4464. The positive completion is an intermediate proof device. The original/completed rooted pressures differ by at most T log(1+A_def exp(-delta_def/T)), with delta_def>=U/4 after choosing lambda and a length-source window. The proof uses invariant per-cell selector mismatch counts, S union successors(S), pin deletion on at most two cells per bad cell, and a partial-code-trace lower bound on the common selector-legal partition function.

[proved] Fixed-marker pinning is now covered. For M>=1 on an L-site ring with ||h||<=J, (L/M)Z_pin<=Z<=L exp(4 beta J L/M)Z_pin. The compression isometry is V psi=direct_sum_i p_i psi/sqrt(M). A bounded rooted seam then gives equality of thermodynamic pressures for M/L tending to a positive density. This is not a spectral comparison at scale 1/L.

[proved] A block-cut and finite-type convolution proof gives canonical length duality for the countable all-length spaces at slopes bracketed inside a trace-class interval. C44's reference length curvature plus the finite norm-error budget verifies the required finite-difference bracket for an interval of scheduled densities, without an IR susceptibility assumption. Together with the uniform tilted embedding bound, p_original and the exact matched thermal logical p_matched differ by C T exp(-delta/T). Pressure-zero and canonical-supremum operations preserve the required o(T^2) accuracy.

[proved application] C45/C49 dressing, full heat trace, two-charge separate Jordan-Wigner extension, C48 coefficient norms and massive winding, and the holomorphic unbounded length-source interface are checked. Use the complete gap including charge tunneling. Strengthen the support exponent before fixing SW constants and choosing lambda last. Keep all scalar terms and the exact massive thermal endpoint.

[assumed, remaining] Only the source-uniform all-scale RG/counterterm/source input (a) and the uniform subtracted scalar thermal input (b) remain in the final free-energy theorem. Source-free beta estimates alone are insufficient for canonical density selection or length compressibility. C55-C57 progress was checked; none claimed these full inputs solved. A complete Gaussian scaling identification, if supplied by (a), gives c_R=c_L=n+1 through C51. No spectral exhaustion or unconditional central charge is proved here.

[numerical] checks.py passes 75 quantum pinning tests, 6400 selector-move checks and nine defect comparisons on matrices at most 64 dimensional. Compression residual 3.29e-16. checks.json records the run. No microscopic Form W simulation.

[proved, file map] RESULT.md is the concise result. PROOF.md contains all new proofs and a hypothesis ledger. PROGRESS.md records completion. All writes stay here; the parent project STATE is outside the authorized scope. No subagents, cluster jobs or network commands. No independent review of this note has occurred.
