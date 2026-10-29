# C65 progress

2026-10-28.

[proved work record] Read S2 and S backgrounds and the complete C57, C59 and C62 text, scripts and saved checks. The project memory was read first. All writes stay in work_C65.

[assumed target] The uniform joint comparison in C59 Lemma TP is not established. C62 supplies an exact trace normalization for its Hamiltonian reference but no uniform kernel estimate for the finite fermion windows.

[heuristic under investigation] A sharp spatial fermion window gives boundary terms with algebraic spatial decay. At positive temperature, its static density bubble appears to have a nonzero cusp at zero transfer momentum even after a smooth infrared subtraction. This would refute the weighted absolute kernel estimate for that particular regulated construction, rather than only leave it unproved. The calculation and its scope are being checked.

[assumed next work] Check smooth regulators and the spacetime-smoothed BM reference separately. Removing an ultraviolet edge obstruction does not by itself establish TP or identify the Euclidean scalar with a Hamiltonian trace.

[proved, sharp-window failure] The cusp is verified explicitly. Its coefficient is -n_Lambda(1-n_Lambda)/(2 pi T), and a smooth low-covariance subtraction changes the bubble only by a smooth function. The loop generates a same-chirality four-field term at order g^2. Antisymmetrization leaves the cusp. A first spatial absolute moment would force continuous differentiability, so the TP fifth-moment preparation and thermal-replacement estimates fail for these finite sharp windows. This does not refute taking their trace limit first.

[proved, repair] Constructed a smooth regulator using a fermionic attenuation isometry relative to the filled sea. It multiplies normal-ordered vertex legs by smooth cutoffs, preserves the quartic degree and the relative form stability bound, and converges in the full Gibbs trace at fixed L,beta. No monotonicity of this new trace sequence is asserted. The smooth static loop has a uniform moment repair, and its infrared shells satisfy the standard cutoff-independent covariance and Gram estimates.

[assumed, remaining] The smooth reference retains C62's scalar normalization, but the all-depth joint normalized tree comparison and actual four-history derivatives remain unproved. BM's smooth fermion regulator removes the cusp too; its spacetime interaction alone does not provide a trace identification. PROOF Section 7 records the remaining estimate.

[numerical, complete] checks.py passes on the first run. It verifies the boundary formula and sharp cusp slope, the growing first spatial moment, smooth-bubble differentiability, and the attenuation identity and inherited form bound. Largest matrix dimension 256. RESULT.md, PROOF.md, SOURCE_AUDIT.md and local STATE.md record the result. No files outside work_C65 were changed.
