# C81 progress

2026-10-28. Work is confined to this folder. No network commands or external
writes. The project memory and requested C75, d24, decisive-attack and
direct-attack source notes have been read. C75's n=15 certificate is an
inherited computer-assisted result, not a new computation here.

[Proved, inherited] Exact squeezing and the cyclic charge-sea projection
hold for any occupied bond. The phase-model reference has K_gap=4 at half
occupancy, while the physical particle field has K_particle=1.

[Numerical, complete] Independent 16-state bond matches C75 entrywise.
All three charges and Hermiticity pass. Full 1350-state ring ED verifies
the cyclic sea projection and J^2 finite-size errors. Two-cell mobile
sectors of dimensions 2320 and 6960 retain illegal strings. Scalar phase
fluid and full-cell frozen selector checks also pass. Outputs are in
checks.json and auxiliary_checks.json.

[Proved from inherited C75 bound, complete] Full finite-K cells replace
the d24 conserved-illegal-bond argument. The boundary-energy comparison
proves nonzero-port frozen selection on compatible rings. A full-cell
thermal argument gives exponential suppression of wrong-cell contributions
to frozen pressure. Neither statement establishes mobile defect stability.
Both derivations and their inherited input are explicit in PROOF.md.

[Complete] RESULT.md, PROOF.md and local STATE.md give the bond, schedule,
transfer audit, all numerical checks and the conditional c=16>15 claim.
The finite-K auxiliary charge grading requires a new matching check,
since the d24 microscopic Z2 assignment does not transfer verbatim.

[Assumed, unresolved] Uniform finite-J vacancy matching, a dressed defect
comparison, and the final infrared/CFT identification are not supplied by
the finite band certificate or the existing ultraviolet normal form.

[Numerical, final validation] Both programs were rerun with assertions and
exited successfully. The cyclic projection error is 1.8207657603852567e-14;
the largest two-cell mobile eigenpair residual is 1.4802183424669638e-11.
Python syntax and local Markdown links pass. No prohibited prose markers
were found. The final run transcripts are checks_run.txt and
auxiliary_run.txt. All outputs remain in this folder.
