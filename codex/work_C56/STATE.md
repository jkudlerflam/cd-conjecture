# C56 local state

Updated 2026-10-28.

[proved, task status] The requested C50 reference comparison (7.1) has not been proved, even for n=1. The physical counterterm remains conditional on the unproved retarded-action RG bounds. Do not describe estimate (a) or input S as closed.

[proved, files] PROOF.md gives the derivations, RESULT.md the outcome, PROGRESS.md the work record. checks.py produces checks.json. All writes are in this folder because the task explicitly prohibits writes elsewhere.

[proved, results]

1. A diameter-only all-degree kernel norm cannot bound a polynomial terminal-tree norm after any fixed field-radius loss. The counterexample is a sparse space-time grid of one-species densities. This is a norm counterexample, not an instability of the physical chain.
2. Actual C48 occupied polymers retain enough geometry for the stronger norm. Safe output rates are alpha/8 and L_0/(128 tau), within its existing cluster budget. This adapts the collaborator's connected-network proof. It is conditional on C48's polymer construction.
3. The all-degree one-endpoint quadratic coefficient defines a small real continuous first-order counterterm after compact-patch initialization, under the stated moment and smoothness assumptions. No uniform higher-order remainder is known here.
4. Given the unique RG trajectory of the collaborator's COUNTERTERM_INDUCTION, uniform bounds on a compact root domain and fixed-scale continuity imply continuity of the infinite-scale relevant counterterm a_0(p'). Finite-cutoff a_0 maps converge uniformly. The proof is compactness plus uniqueness and avoids root derivatives. The final physical map is nu=A^(-1)(u-a_0), retaining the exact initialization convention. Brouwer then applies as in C50.

[assumed, unresolved] Derive the normalized marked-tree bounds for actual arbitrary-degree retarded endpoints, including wave-function normalization, on one common history domain. Section 3 of PROOF gives sufficient comparison estimates for (7.1). The collaborator's G1-G2 hypotheses and exact physical initialization must be verified before the full counterterm construction applies. Sources, thermal scalar estimate and physical finite-size transfer are separate.

[proved, provenance] Use C39/C42 finite bands, C44/C45 mixed margins and exact charges, C48 matching, C50 finite integration and field-radius conversion. The collaborator already supplied Gaussian local mode diagonalization, UV integration, degree reparameterization, marked-tree transfer lemmas, the conditional coupled recursion and finite-cutoff Brouwer existence. They are cited in PROOF and were not claimed as new C56 results.

[numerical] Small checks passed. No heavy computation, physical-chain simulation, or independent referee was used.
