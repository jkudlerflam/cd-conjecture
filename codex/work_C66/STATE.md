# C66 local state

[proved, task record] Completed on 2026-10-28. RESULT.md is the summary and PROOF.md the derivation. Work concerns the n=1 C44/C45 mobile wrapper, with one payload fermion and one length fermion. The microscopic infrared theorem remains unproved.

[proved, explicit specialization] The four physical projected payload bonds Q_ab=lambda^2 K2_ab give A=Q00, B=Q10-Q00, C=Q01-Q00, D=Q11-Q10-Q01+Q00. Static action contains A+m B+m' C+mm' D. V_l=tr(D)/4 is the scalar length interaction. PROOF (13)-(14) give all six first-order forward coefficients, including contractions of six- and eight-field terms. Fully matched coefficients are defined by finite matching, not numerically evaluated.

[proved, scheduled reference] At n=1 Gamma=1/256, t_p=sqrt(1+eta^2)/256, v_p=sqrt(4t_p^2-mu_p^2), j_l=W=lambda^2 J0/s, and v_l=2sin(pi nu) in W units. J0=2g^2 q^4 D_aux^(-6)/(1+q^2+q^4)^2. The four velocities are (v_p,-v_p,v_l,-v_l). Constants must pay inverse v_min, not inverse speed differences. The numerical verifier instead sets length hopping equal to payload hopping and is not a test of this scheduled velocity hierarchy.

[proved, stability] The isolated attractive length chain has linear fixed-density energy for V_l<=-2j_l; require a strict reserve |V_l|<2j_l. The scheduled scalar bound is |V_l|<=4W d_eq/J0. Mixed terms require the full Gaussian stiffness bound, or its block Schur complement, in addition. Gaussian positivity alone does not establish scaling dimensions or microscopic RG.

[proved, conditional comparison] PROOF Sections 6-7 specialize corrected C61 Lemma N and C50 (7.1). The free source difference D_O(0) stays separate. A four-label lift of a species-independent routing bound costs a fixed exponential endpoint factor. This is an implication, not a proof of the missing routing. Comparison uses separate forcing d_for and marginal size eta; explicit relevant-tail and geometric-sum constants are given. C56's root map remains conditional on its RG hypotheses.

[assumed, next work] Prove one simultaneous spatial routing with paired cut coverage, bounded line/endpoint loads and controlled coefficient mass; retain normalized history marks and multichannel signed reference cancellation. Then establish source-uniform critical comparison. Scalar thermal P/TP remains separate. Do not mark input S or microscopic c=2 complete.

[numerical] checks.py passes, with output CHECKS.md. It tests finite algebra, not an interacting trajectory. Maximum error in 18 forward-coefficient polynomial derivative comparisons is 1.36e-14. No full-chain simulation or large computation was run.

[proved, scope] All files are local to work_C66. No subagents, network commands or writes elsewhere. Primary RG papers were checked with the browsing tool. Parent STATE.md was not modified. No independent review is claimed.
