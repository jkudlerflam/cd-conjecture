# C51 local handoff

2026-10-28.

[proved conditional, main result] PROOF.md establishes c_R=c_L=N=n+1 for a complete compact Gaussian scaling limit with positive 2N-dimensional quadratic form, arbitrary velocity ratios and density-current mixing. c_cd=(c_R+c_L)/2=N. The O(N,N) canonical transformation is explicit. The Sugawara proof works at nonzero zero-mode charge. This is a conditional continuum identification, not a microscopic construction theorem.

[proved, finite-size qualifications] A shifted compact lattice x=B(m+alpha_L) gives full current Fock multiplicities product p(k_i). Irrational filling moves alpha_L and can spoil a naive ground-energy Casimir limit. Correcting the quadratic zero-mode energy restores the Gaussian Casimir term. Subsequence convergence and uniform theta bounds are proved. An o(1/L) spectral matching/exhaustion contract F is needed for microscopic towers. The actual wrapper compactification lattice is not derived. Individual Gaussian charges need not be exact lattice symmetries.

[proved conditional, thermal theorem] The physical coefficient is A=(pi/12)sum_i 1/u_i over all 2N chiral speeds. It reduces to (pi/6)sum_a 1/v_a only with paired equal speeds. A alone is not an unweighted central charge. Contract P requires a ground-subtracted scalar shell comparison |delta_h(T)|<=C q_h^2 eta_h min(1,(T/q_h)^p), p>2, eta_h->0, plus UV control. Its sum is o(T^2). For eta_h=O(q_h^omega), the relative power is min(omega,p-2), with a log at equality. P is a sufficient unproved microscopic estimate, not an inference from normalized correlations.

[proved, prior work and cross-read] C46 already supplies the diagonal Gaussian pressure, common-bond quotient transfer and basic Sugawara bound; cited, not claimed as new. C48/C49/C50/C52 were empty initially but their progress was read after it appeared. C50 PROOF explicitly leaves all-scale comparison and source bounds open, so no microscopic LL theorem can be imported as already proved. Its finite-scale scalar norm does not prove P. C48's rooted/enlarged pressure result retains physical-ensemble qualifications.

[numerical] checks.py and checks.json contain small standalone controls. Maximum canonical error 3.56e-14, Virasoro error 6.22e-15. Fock cutoff dimension 139. All assertions pass without warnings in the final run. No microscopic spectrum, entropy or pressure was computed.

[proved, files and scope] RESULT.md is the concise conclusion; PROOF.md is the derivation; PROGRESS.md records the work. No network, subagents, large computation or edits outside work_C51. Parent project STATE.md was not updated because the task forbids writes outside this folder. The construction remains conditional. Next work is a complete current/field scaling identification, scalar thermal estimates P or full physical finite-size matching F, not another free-boson mode count.
