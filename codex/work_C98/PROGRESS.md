# C98 progress

2026-10-29. All work is confined to work_C98. No agents or cluster jobs are used.

[Assumed inputs] Read the project memory, the requested bound background, and the relevant C14, C35, C69 and survey results. The sharp proposed coefficient is pi^2 M log_2(d)/6 in the energy deficit. C35 already proves the logarithmic comparison for the symmetric SU(N) families and critical XXZ, conditional on their stated continuum inputs. Those proofs will be cited, not repeated.

[Proved, initial reduction] With D(L)=E(2L)-2E(L), the Casimir expansion gives D(L)=pi A/(4L)+o(1/L). A useful gluing theorem must upper-bound D(L), not lower-bound it. Product gluing gives a constant boundary cost and misses a factor of L.

[Completed] Checked exact finite-ring XX formulas, small XXZ/Heisenberg rings, fundamental SU(N) rings and the inherited symmetric-representation formulas. The no-remainder sharp finite-size statement fails at Heisenberg L=4. The requested asymptotic statement is not falsified.

[Open] No proof or counterexample to the universal logarithmic Casimir bound has been obtained.

[Proved] Product gluing has O(M) cost; raw XX closure has a nonzero boundary energy; XX has exponential half-chain Schmidt rank. The smooth partial-copy swap has an extensive cost, calculated explicitly for XX and Heisenberg. Reflection concavity and a bounded boundary deficit are insufficient at the level of energy sequences.

[Proved] Derived the exact two-seam susceptibility formula for the doubling defect using finite-temperature free energies, then taking beta to infinity at fixed L. This avoids differentiability assumptions at level crossings. A quantitative positive-part integrated-response bound is sufficient for the logarithmic Casimir inequality, but remains an additional unproved hypothesis.

[Numerical] checks.py and response_check.py passed with python3. The initial executable name python was absent and was corrected. Ring eigenpair residuals are below 6e-12; the two-seam identity quadrature error is below 1e-15. No large computation was run. CHECKS.md and RESPONSE_CHECK.md contain the full results.

[Completed deliverables] RESULT.md gives the outcome, PROOF.md all derivations and qualifications, and STATE.md the local handoff. No file outside work_C98 was modified. The project-level memory was left untouched as required by the task restriction.
