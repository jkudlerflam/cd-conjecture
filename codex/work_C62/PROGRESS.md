# C62 progress

2026-10-28.

[read] S2 and S backgrounds, project memory, C57 and C59 proofs and source audits, and the collaborator's thermal winding, localization, moment-reserve and four-trajectory notes. C59 Lemma TP is the target. All writes stay in work_C62.

[assumed, target not proved] The task requires a uniform comparison of complete interacting effective actions, not only localized couplings, and a scalar identity with a Hamiltonian trace. The existing thermal mark allocation is conditional on ordinary renormalized-tree routing. It is imported, not rederived.

[in progress] Checking the primary BM reference construction, its ultraviolet regulator and the passage from normalized Schwinger functions to a thermal trace. Investigating whether the mixed scale sums actually follow from the available short-memory estimates and where a shared reference-defect and thermal mark need additional bounds.

[read] The December 2024 Renzi-Giuliani-Mastropietro-Santini announcement remains a relevant prior-work lead. An announcement is not being used as a theorem with verified hypotheses. No author contacted, no subagents, no shell network commands.

[proved partial] PROOF.md now constructs a trace-defined Hamiltonian reference. Its finite fermion-window Gibbs traces converge by min-max and form-core approximation, and the exact finite-regulator Grassmann identity retains its scalar. The reference has the correct matched v,K. This is not a proof that its kernels satisfy TP's common all-depth norm.

[proved] The parameter conversion between the BM and Hamiltonian conventions is explicit. Finite-cutoff currents have boundary occupation terms in their commutator; direct oscillator substitution at finite cutoff is invalid. A double scale sum retains both required powers using only one short-memory factor. A separate analytic-flow example shows why affine bare interpolation between matched critical inputs need not remain controlled.

[unproved, stopping point] The full joint marked GN majorant, actual four-history estimates and compatible vacuum derivatives are not established. PROOF Section 7 states the obstruction. Thermal c=1 for interacting XXZ remains conditional. No counterexample to its physical prediction was found.

[checks complete] All controls in checks.py pass, including independently assembled finite fermionic reference matrices compared with the oscillator trace. Largest dimension 256. RESULT.md, SOURCE_AUDIT.md and local STATE.md record the partial outcome and its limits.
