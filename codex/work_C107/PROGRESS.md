# C107 progress

2026-10-29.

[Proved] Read the project STATE.md, the complete updated bound background, and REPORT.md Sections 9.68 through 9.75. Read the prior C16 compactness discussion and C98 boundary-energy obstruction. The task confines all modifications to work_C107, so the memory handoff will be a local STATE.md.

[Proved] The literal set M(h)=1 is not compact because scalar shifts and coboundaries are unbounded. The quotient sphere is compact. Finite-volume Casimir quantities are continuous with a Lipschitz constant of order L squared. Uniform convergence on a neighborhood containing an ordinary critical point and its massive perturbations is impossible.

[In progress] Deriving an exact XY-chain obstruction to monotonicity, checking length-sector dependence of the proposed limit, and constructing a compact one-parameter energy countermodel satisfying the elementary concavity and boundary estimates but with an unbounded asymptotic coefficient. The energy countermodel is not a spin-chain construction.

[Open] Upper semicontinuity of a physical Casimir extension, local boundedness near all singular bonds, and a finite universal F(d) have not been proved or disproved.


## Completion

[Proved] Completed the quotient compactness argument, finite-L Lipschitz estimate, subadditivity audit, XX and gapped classical Ising sector examples, and explicit failure of uniform convergence near critical Ising.

[Proved] Derived the critical XY correction a_gamma=7pi^2(4gamma^2-3)/(240gamma). Its sign change refutes universal monotonicity. Its divergence at gamma=0 refutes a universal decreasing correction Q_L+C M/L^2.

[Proved logical counterexample] Constructed a compact analytic energy family with finite pointwise Casimir limits and an unbounded coefficient. It preserves concavity and the elementary boundary/additivity estimates but has no asserted local Hamiltonian realization.

[Proved] Identified local boundedness everywhere on the quotient sphere as equivalent to the desired finite bound. The canonical upper semicontinuous majorant is allowed to be infinite, so defining it does not close the argument. Gave conditional finite-cover and finite-scale formulas for making a bound explicit.

[Numerical] checks.py completed successfully. Nine small spin matrices agree with the free-fermion formula to maximum error 1.95e-14. Scalar sums and quadrature check the asymptotic coefficient. No heavy computation was run.

[Work record] An initial Python file-writing command had a quotation syntax error and wrote no PROOF.md. The corrected command succeeded. All requested reports and local memory are now written. Only this folder was modified.

[Open] No general finite F(d), upper semicontinuity theorem for physical A, or actual divergent-A bond sequence was found. The final report states these limits explicitly.
