# C129 progress

2026-11-07.

[Proved] Field-sector reduction, DM gauge transformation, the field-dependent optimized bond scale, an exact twist inequality, and a minimizing flux have been derived. The proofs are being written in PROOF.md.

[Proved] The unrestricted ground energy is extensive for J_z <= -sqrt(J_perp^2+D^2), for every field. This closes the original target on that region without requiring the stronger sectorwise droplet lemma.

[Proved, literature audit] The premise about Kozlowski's half-filled expansion is incorrect. Proposition 3.4 is the all-order result in the finite-Fermi-boundary regime. Theorem 3.5 gives convergence at half filling for -1 < Delta <= 1, without a finite-size rate. Corollary 4.3 already states the energy expansion below half filling. The paper explicitly flags missing uniform energy control near half filling after Lemma 4.1.

[Assumed, unresolved] No constant for the entire U(1) class has been obtained. Parameter-uniform estimates and the ferromagnetic crossover from above remain to be established.

All writes and computations are confined to this directory. No cluster connection or large local computation is planned. Larger runs will be supplied as SLURM-ready scripts.

## Completed work

[Proved] PROOF.md contains Steps 1-3 with the exact trigonometric twist estimate and explicit Taylor remainder, a minimizing flux, both free-fermion parities, the polarized and classical regions, the dilute estimate Q<=2pi^2 n^3/L, and the crossover estimate Q<=(Delta+1)L^2. It also records the saturation-field bound and the literature dependency audit.

[Numerical] All 250 small Bethe/ED comparisons passed. Twist checks through L=12, all 144 requested survey-point/length field checks, and 72 variational checks passed. Results and commands are in CHECKS.md. The largest sampled doubling ratio is 4. No numerical output is used to certify a theorem.

[Assumed, unresolved] The full C_U1 bound and the requested uniform pieces (a) and (b) remain open in this run. Half filling for the whole massless interval is absent from the cited finite-size expansion. Delta-uniform remainder constants and crossover control outside the explicit variational bands have not been proved.

[Work record] RESULT.md and local STATE.md are written. run_typhon.slurm prepares the requested large scan, which was not submitted. All files are local to work_C129. No outside file was changed.
