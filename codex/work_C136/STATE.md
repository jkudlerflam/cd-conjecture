# C136 handoff

2026-11-07. All files are local to work_C136. Read RESULT.md, then PROOF.md. No independent referee pass on this new proof. C134's inherited proof does have an independent referee report, despite stale predecessor status lines.

[Proved] Task 1 is complete at fixed anisotropy gap delta>0, including every filling and both parities. For Delta in [-1+delta,0], all sectors and all L>=2,

    |F_L(n)-J L e_bulk(Delta,n/L)| <= 10^15 delta^(-25/3) J/L.

PROOF Section 6 gives a smaller explicit K_all(delta). C129's reductions give the physical all-field/all-DM theorem

    L(E_(2L)-2E_L) <= 6 * 10^15 delta^(-25/3) M(B_field).

M includes the field. The two-oriented-bond convention applies at L=2. Odd doubling uses (2L,2n), not C134's (2L,L).

[Proved, new ingredients] The clipped cumulative fixed-point comparison gives sup CDF error B/(2L), independent of filling. Finite-boundary dressing subtracts h(Q), with epsilon_Q(+-Q)=0, and compares at the finite counting boundary Q_L. Its continuum mass can differ from n/L by at most B/L. Bulk curvature is bounded by pi sin gamma. Restricted positive resolvents have an exponential majorant from the beta=1/2 kernel moment. This proves kappa r_infinity<=r_Q<=r_infinity uniformly in Q, and exponentially localized boundary corrections through two derivatives. These imply uniform variation of epsilon_Q'/r_Q. Combining with the cumulative discrepancy, C134 clipping, and explicit tails proves the all-filling theorem. Section 6.5 checks the simplified constant arithmetically.

[Proved] Section 5 has smaller explicit constants below n/L<=1/2-eta, uniform down to zero density. Section 9 also gives a smaller-constant high-field theorem with |B_field|>=(15/16)J(1+Delta), using a convex-combination sector lower bound. Neither is needed to state the all-filling compact theorem.

[Assumed, unresolved] Task 2 and the full requested uniform constant through Delta=-1 remain open. C134's cutoff scale is epsilon^(-1), beyond the epsilon^(-1/2) band. Its large-size coefficient still grows as epsilon^(-3/2), even after its small-size branch is discarded. The new constants also diverge. An explicit uniform bulk-error estimate beyond epsilon L^2>=A, or a direct one-sided defect estimate there, would close the remaining gap. PROOF Section 8 states it and Section 9 gives the conditional assembly. Do not promote the compact theorem to a theorem on all of [-1,0].

[Numerical] 1,032 Bethe pairs pass, all sectors through L=12 and at 16,17,32,33,64,65, six sectors at 128,129,256,257. ED error <=1.25e-14 at every sector through 12. Max Q=2.3431457505076203, the two-site free value. Field, twist, kernel, and arithmetic checks pass. Only 43 of 54 independent continuum diagnostics are resolved. Eleven near-ferromagnetic cases remain underresolved, with explicit errors retained. CHECKS.md details the scope, including the untested interacting large-size clipping branch.

[Assumed, pending] run_typhon.slurm supplies the full 98,304-pair grid through L=256. run_intermediate_typhon.slurm supplies denser continuum diagnostics. Neither was submitted. Typhon is inaccessible. Optional positive-Delta work was not pursued because Tasks 2 and 3 remain incomplete. Project STATE.md was not edited under the outside-write prohibition.
