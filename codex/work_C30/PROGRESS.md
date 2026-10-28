# C30 progress

2026-10-28. Complete for the selector task. Read RESULT.md.

[Proved] A nine-letter arrow selector with one C7 walker on each arm selects exactly P_L L^m M R^m P_R for any m>=1. The only variable coefficient needed in this selector family is the P_L chemical potential. With f(r)=-2t cos(pi/(r+1)), d_r=f(r+1)-f(r), choose e=(d_(m-1)+d_m)/2 and mu=(2m+3)e-2f(m) for m>=2. The all-length excess is at least eta_m(|a-m|+|c-m|), eta_m=(d_m-d_(m-1))/2. For m=1 use mu=-10t,e=-2t,eta=t. K=32t excludes arbitrary illegal ring sectors by at least 4t per violation relative to eL. Frozen carrier dimension is 120. Compatible-ring code and gap are explicit.

[Proved] At mu=(2m+3)d_m-2f(m), all four pairs in {m,m+1} squared coexist. Thus the independent-arm crossing has zero-cost asymmetric cells. Changing the right arm also changes the number of logical modes. This is not C8's fixed-payload length fluid.

[Proved] A positive symmetric core with one neutral helix qubit on each L/R site has local dimension 84 and exactly 2^m code states at k=2m+3 occupied sites. Its initialized arm and illegal-sector gaps are explicit. It satisfies C8's positive-core hypotheses. Applying C8 gives physical lengths k+1,k+2, with the same symmetric occupied core and neutral padding. Dimension is 2028 before endpoint compensators, 8112 including them. The helix core encodes m in its phase and pin, so the one-affine-chemical-potential property is asserted only for the frozen walker family.

[Proved, finite matrices] A Schur-complement bound extends C25's ordered-limit mechanism to this noncommuting carrier core: first w/U small but nonzero, then core scale s below the resulting fast-wrapper gap g. Compensated second-order source bonds converge to the bare source bond with errors O((w/U)^2)+O(s/g), with constants depending on finite core gaps. Exact equality at finite parameters is not asserted. Conservation of non-code-dimer count makes the mixed u lambda term vanish for this implementation. Full wrapped-switch source ED was not performed.

[Proved conditional graph result] A hypothesized fixed-payload step-two exchange fluid with odd b has one length degree of freedom, constraint L=bM+2r, M=L mod 2, and a connected physical marker graph. Its complete free spectrum and modified winding gluing are in RESULT Section 5. A microscopic symmetric pair-transfer parent was not constructed. The actual mobile design uses step-one padding instead.

[Numerical] checks.py and checks_results.json are complete. Full walker spectra for arms 1..8 (dimension <=64) agree to 6.22e-15. Selection inequalities pass for m=1..12 and both arms 1..48. All nine-letter ring sectors through L=5 (59049 states) obey the illegal-sector bound. Thirty arm-helix matrices for m=1..4, r<=9 (dimension <=512) pass gap/rejection tests. All legal initialized core patterns at k=5,7,9,11 have the stated unique structural kernel. Step-two graphs of dimensions 36,58,120 agree with the exact spectrum to 7.55e-15. No central-charge estimate is made.

[Assumed remaining physics] The routeA_compile bare carrier's c=n evidence covers n=1..4, displaced seas only n=2. Its all-m band count, momentum margins, and finite-coupling infrared transfer remain open. C13/C20 require actual nonresonance, positive stiffness, local massive-sector elimination, and gluing. Exact total charge does not remove interflavor pair-transfer resonances. C25's finite-w source mismatch and scalar length interactions remain in the scale budget.

[Proved clarification] The compile README's 5 percent asymmetry changes arm hopping, not arm length. It gives no geometric arm-length tolerance. The present frozen construction selects equal lengths, so no geometric fallback was needed.

Only this directory was modified. No network, cluster, subagents, or external messages were used. The project STATE.md is outside the permitted write scope, so local STATE.md carries the handoff. Sources read include the requested route-A notes, C7, relevant C8 finite-cell and wrapper proofs, C13, C20, C25, and concurrent C21/C22/C23 progress.
