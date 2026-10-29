# C114 progress

2026-10-29. All writes are confined to work_C114.

[Assumed inputs] Read the project memory, all updates in the bound background, REPORT.md Sections 9.68 through 9.78, and the complete text and Python sources in C98 and C107. Their settled formulas are inherited, not counted as new results. The vacuum-compatible ring expansion is assumed throughout coefficient extractions.

[Proved] An energy-preserving lift of an L-ring ground state to a kL ring would imply E(kL) <= k E(L). This contradicts the positive Casimir defect for A > 0. Convex combinations of larger-ring energy densities require a correction of order A/L^2. This invalidates the exact-lift premise, not the proposed bound.

[Proved] An exact finite-range positive decomposition plus a spatial coboundary wraps without a seam error on every sufficiently large ring. It would force E(L) >= L e0 and hence cannot describe a positive Casimir coefficient. The same obstruction applies to a finite-range Hamiltonian commutator after taking its expectation in a ring ground state. Quantitative approximate decompositions remain under investigation.

[Operational] First Typhon connection used ssh -o LogLevel=ERROR with a 15-second connection timeout. It exited 255 with 'Connection timed out during banner exchange'. No job was submitted and no computation ran on the login node. A later retry is planned. Small local checks and preparation of a SLURM suite can proceed.

[Heuristic candidate] A coherent sum of translations of two concatenated ring ground states is a different trial state from a mixture of seam positions or a product of partial swaps. Its Gram matrix and energy matrix will be checked explicitly. No inverse-length bound for this trial family is assumed.

[Open] No explicit universal C and f(d), and no counterexample to the target, have been obtained.

## Completed derivations and checks

[Proved] Derived the coherent reconnection inequality D_L<=q_L/S_L. The Gram matrix is the interval-purity matrix and Re K=q_L I. For real h and psi with q_L>0, this is the optimum over the whole coherent prefix-swap span. The uniform coherent upper bound itself does not require real h or any internal symmetry.

[Proved] On XX, a direct correlation-matrix estimate gives p_r<=min(r,L-r)^(-2/pi^2), using the number variance and the odd harmonic sum. This proves L times the optimized coherent cost diverges. No CFT entropy formula is needed for this exclusion.

[Proved conditional] A support-r positive certificate with norm residual <=B M f(d)/r^2+o(r^-2) would give A <=(24B/pi) M f(d). No such certificate was constructed. Exact certificates or faster residual decay contradict a positive Casimir coefficient. Averaging seam positions does not reduce a fixed seam expectation.

[Numerical] checks.py passed on XX, both XXZ signs, Heisenberg, and SU(3), with the largest full-ring space of dimension 729. Maximum eigenpair residual is 1.33e-13. It checks the failed zero-cost lift, the coherent variational inequality and optimum, the bulk-deficit sign, the bare norm bound, and exact XX energies. No approximate SOS certificate was numerically tested because none was constructed.

[Numerical] purity_check.py passed, using XX correlation matrices of dimension at most 32 and three arbitrary complex qubit bonds. single_ring_check.py passed through qubit L=8 and SU(3) L=6. All local numerical jobs were small and completed in under a second each. No heavy computation was moved off Typhon.

[Operational] A second SSH attempt after the first small-ring calculations, and a third after the analytic write-up, both used the same 15-second connection timeout and exited 255 with 'Connection timed out during banner exchange'. No job was submitted and nothing ran on the login node. The expanded suite is prepared as run_ed.slurm and checks.py --full; the latter refuses to run without SLURM_JOB_ID. That suite remains unrun.

[Completed deliverables] RESULT.md and PROOF.md are written, together with the generated check reports and local STATE.md. Only this folder was modified. The parent memory and shared reports were not changed because the task prohibits writes outside this folder.

[Open] The general bound and the adapted local-unitary reconnection estimate remain unresolved. The negative results concern the specified exact lift, exact finite-range certificate, and coherent prefix-swap trial family. They are not counterexamples to A <= C M f(d).
