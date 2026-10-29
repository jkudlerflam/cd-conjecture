# C118 progress

2026-10-29.

[Assumed] Read the project STATE.md, every update of specs/_bound_background.md, REPORT.md sections 9.68 to 9.79, and every file in work_C116. Its settled constants and results are used without rederiving them. All writes and Python execution are confined to work_C118. The parent STATE.md cannot be updated under the task restriction; a local STATE.md will carry the handoff.

[Open] The requested estimate for R has not been proved. The fixed-L zero-temperature limit requires a uniform bound on L(E_(2L)-2E_L)/M, stronger than a bound on its infrared limit.

[In progress] Audit seam averaging and finite-rank insertions, derive the exact relative-entropy identities, and test convexity and thermal-remainder candidates on full spectra. Large XX checks use linear-cost free-fermion recursions. Larger interacting diagonalizations will be provided as SLURM-ready scripts, without attempting Typhon access.

[Proved] Completed a full-trace XX estimate, R <= 2 log 2 + (5 pi/3) beta/L, including parity projections. It extends to independent XX tensor products. Also proved R <= log d for every diagonal classical bond, without assuming a Hermitian spatial transfer matrix.

[Proved] Derived exact seam Jensen and relative-entropy comparisons. Averaged seam norms remain at least 2M on XX, Heisenberg and SU(3). Explicit partial-swap commutator witnesses show the missing continuous symmetry. XX has both signs of R curvature. The continuum scaling form and a single-temperature bound leave slow-sector multiplicity uncontrolled without microscopic input.

[Numerical] Full-trace tests and all information identities completed. XX paired length up to 8192, interacting ED up to Heisenberg length 8 and SU(3) length 6. The first two local runs checked the proof and additional partial-swap tests; the final run also recorded the explicitly unproved trial extension to the interacting models. Full trace discrepancy is 2.14e-14, information identity residual at most 2.60e-12. Failed candidate inequalities are listed explicitly in CHECKS.md.

[Assumed published context] Consulted the primary abstract of Wolf et al., arXiv:0704.3906, through the browsing tool. It states that the quantum thermal mutual-information area prefactor may grow with inverse temperature. The relevant Gibbs variational inequality is derived independently in PROOF.md.

[Work record] Wrote RESULT.md, PROOF.md, CHECKS.md and local STATE.md. A SLURM-ready larger run is supplied but unsubmitted. No Typhon access, large local computation, or write outside this folder. The general ratio estimate remains open; no independent referee review was performed.

[Work record] Final validation passed Python syntax, saved-result counts, prose restrictions, and the SLURM guard. Invoking --cluster without SLURM_JOB_ID exits with the expected message, "Large ED requires a SLURM allocation." No large run was started.
