# C71 progress

2026-10-28.

[proved work record] Read the project memory, S2 and S backgrounds, and C57, C59, C62 and C65, including proofs, source audits, scripts and saved checks. Read the collaborator's winding, thermal localization, moment reserve and four-history lemmas. All writes remain in work_C71. The project memory outside this folder will not be modified.

[assumed target] C59 Lemma TP requires a common full-kernel norm, a vanishing zero-temperature reference defect, a mixed thermal replacement with p > 2 + theta, and compatible complete vacuum derivatives. C65's smooth sea-preserving attenuation and trace convergence are retained. Neither a trace limit nor a smooth individual shell supplies uniform ultraviolet preparation.

[heuristic under investigation] Audit the high-frequency and large-momentum preparation of the instantaneous Hamiltonian reference. In particular, check which coefficient cancellations are required before taking absolute moments, and whether smooth cutoff removal is bounded in the precise off-shell norm. Also distinguish finite-regulator infrared matching from matching only the limiting trace.

[proved partial result, write-up in preparation] For the smooth reference covariance, the absolute two-line chiral ultraviolet loop has integral (4 pi A)^(-1) log Lambda + O(1) on a fixed positive-time strip, even after a fixed smooth infrared subtraction. Spatial smoothing of the interaction does not cure this particular absolute-before-integration majorant. This is not a divergence of the integrated effective coefficient.

[proved partial result, write-up in preparation] The integrated density bubble at zero temperature is p/[2 pi (A p - i Omega)] on transfer momenta inside the regulator plateau. After spatial interaction smearing and infrared subtraction it has every fixed absolute spacetime moment. A derivative identity moves its short-distance singularity onto the smooth spatial potential. For the thermal difference, equal-image singular products cancel, leaving only distinct-image products with integrable singularities. This supplies arbitrary fixed thermal powers for this second-order prepared coefficient, uniformly in the smooth regulator, including two speed derivatives. It does not prove the all-order ultraviolet preparation or TP.

[proved, completed write-up] PROOF.md now proves the logarithmic absolute-majorant obstruction, the integrated-loop moment bound and the uniform differentiated thermal coefficient estimate. It separates this method obstruction from C65's actual sharp-regulator norm failure. Section 6 states the next ultraviolet lemma and the remaining full-kernel/history/scalar estimates. No q_h^theta defect theorem for XXZ has been established.

[numerical, complete] checks.py passed on its first run. The logarithm's doubling increment agrees within 3e-17, the largest dynamic Matsubara extrapolation discrepancy is 2.82e-11, and the finite image identity residual is 7.95e-17. The winding quadrature is an explicit high-covariance control, not an interacting flow.

[proved scope] RESULT.md, SOURCE_AUDIT.md and STATE.md record the outcome and limitations. All files remain in this folder. No subagents, shell network commands, large computation or independent review. TP remains unproved; the next concrete work is an all-order preparation that retains the cancellations demonstrated here, followed by the common-domain infrared comparison.
