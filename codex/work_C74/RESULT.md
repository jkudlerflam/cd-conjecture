# C74 result

[assumed, unresolved target] C60 (17) and corrected Lemma N remain open for general GN trees. I did not identify a contraction pattern forcing unbounded counters under every reorganization. The results below close a restricted contraction family and one counter in the general problem.

[proved, general derivative counter] Factoring the momenta of surviving fields before each Gaussian contraction gives an all-depth induction with at most two localization derivatives on each original field. Tree covariances have at most four, and determinant rows or columns have at most two. The proof retains the determinant and applies R_m(P_a F)=P_a R_(m-a)_+ F after each merger. This does not bound the coefficient moments hidden in F. See [PROOF.md Section 2](PROOF.md).

[proved, contraction formula] A projector moved through C_W F(q)=integral W(c,q)F(c,q) dc produces additional contracted jets. For quartic child and parent the exact identity is

    R_1^q C_W R_1^(c,q) K
      = integral W(c,q)[K(c,q)-K(c,0)] dc
        + integral [W(c,q)-W(c,0)][K(c,0)-K(0,0)] dc.

The second term is absent from a fixed-kernel flag collapse. A legal two-quartic BBF merger has a nonzero regenerated contact, explicitly negative in (15) for a temporal cosine endpoint and the reflection-symmetric free lattice covariance. That merger has a bounded repair, so it is not a counterexample to (17).

[proved, arbitrary-depth ladder theorem] Sections 4-6 establish (17) for quartic pair ladders with a retarded initial endpoint and local forward vertices thereafter. Both terms above have an exact spatial representation. Their path records are respectively 0 to h, and 0 to s_1 together with s_1 to h. Every required scale cut is crossed once. Subsequent contractions multiply in the conserved pair momentum and preserve its zero. The maximum endpoint moment, tree-line moment, covariance derivative and field derivative orders are all one, independently of ladder depth. The bound is

    C^L ||kappa||_exp product_j |lambda_j| 2^h,

with no constant per unit scale span. The one-by-one residual BBF determinants, weighted-line arguments and anchor integrals are retained in the same representation.

[definition] This is the derivative-decorated root estimate with normalized covariances and explicit field-normalization ratio products factored out, as in (17). It is not an untested kernel total-variation bound uniform in an independent root scale.

[assumed, remaining assertion] General branching mergers retain several independent port variables and generate the jet terms of PROOF (8). The general derivative induction does not control their accumulated endpoint or line moments, or assign their gains to all cuts. Source comparison and the bound for D_O(u)-D_O(0) also remain separate. Neither the nonzero contact nor a large intermediate Taylor order proves an unavoidable obstruction.

[numerical, finite exact-arithmetic checks] `python3 -B checks.py` passes 64 contraction-jet identities, 32 child-parent identities, 192 moment-differentiated identities, 405 spatial splits, 826 field-counter cases, 240 repeated-multiplier cases, 1,024 ladder cut records and 40 BBF checks. Multiplier chains were checked through depth 40 and cut records through depth 64. These are finite controls, not verification of general GN routing.

[record] PROGRESS.md, PROOF.md, STATE.md and the check files are local to work_C74. No outside files were modified, and no independent review is claimed.
