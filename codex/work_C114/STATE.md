# C114 local state

2026-10-29. Investigation written up. The universal bound remains open.

[Open] No explicit universal C and f(d) for A <= C M f(d), and no counterexample to that target. General adapted seam unitaries remain unconstructed and unexcluded. Do not describe C114 as proving the Casimir bound.

[Proved] New coherent-reconnection bound D_L <= q_L/S_L, with q_L=2[tr(h sigma_L tensor sigma_L)-E_L/L] and S_L=sum_(r=0)^(L-1) tr(rho_r^2). Construct prefix swaps U_j on two copies of a ring ground vector. Gram matrix G_ij=tr(rho_|i-j|^2), and Re K=q_L I for the energy-defect matrix K. This upper bound has no internal symmetry assumption. For real h and real translation eigenvector psi with q_L>0, K=q_L I, so q_L/S_L is the optimum over the entire coherent prefix-swap span.

[Proved] Entire coherent prefix-swap span fails to achieve O(1/L) on XX. Gaussian interval purity satisfies p_r <= min(r,L-r)^(-2/pi^2), proved by p_r<=exp(-2 Var N_r) and Var N_r>=log min(r,L-r)/pi^2. Since q_L=2/[L sin(pi/L)]>=2/pi and S_L=O(L^(1-2/pi^2)), L times the minimum trial cost diverges. This does not exclude more general reconnections.

[Proved] Exact energy-preserving lifts contradict E(kL)-kE(L)>0 for A>0. Longer-ring convex comparisons require an A-dependent order-L^-2 error. Exact finite-range positive-plus-coboundary certificates, even with finite Hamiltonian commutators, force E_L>=L e0 and cannot give a positive Casimir term. A support-r certificate with residual <=B M f(d)/r^2+o(r^-2) would imply A <=(24B/pi) M f(d); its construction is open.

[Assumed] The vacuum-compatible ring expansion and inherited Bethe bulk energies are inputs. C98 and C107 were read fully. Their scalar almost-additivity results, gradient exclusions, and partial-swap failures are inherited, not new discoveries here. Standard CFT purity scaling is used only as calibration; the XX obstruction is independent of it.

[Numerical] checks.py, purity_check.py, and single_ring_check.py passed. Full-ring ED uses at most 729 states, comparing qubit rings 4 and 8 and SU(3) rings 3 and 6. Additional single-ring ED reaches qubit L=8 and SU(3) L=6. Maximum eigenpair residual 1.33e-13. XX correlation matrices reach ring L=64 with matrices no larger than 32. Three arbitrary complex bonds test the symmetry-free identities.

[Operational] Three Typhon connections timed out during banner exchange, exit 255. No cluster job submitted. run_ed.slurm and checks.py --full prepare a larger suite; full mode requires SLURM_JOB_ID. That suite is unrun. No heavy numerics ran locally. Only work_C114 was modified, so the parent project STATE.md and shared REPORT.md were left untouched. This file is the local memory handoff.

[Files] RESULT.md is the verdict. PROOF.md contains full derivations, with the new theorem in Section 5 and analytic XX obstruction in Section 6. CHECKS.md, PURITY_CHECKS.md, and SINGLE_RING_CHECKS.md are generated numerical reports. checks.json stores the small-ring values. PROGRESS.md records the work and cluster attempts.
