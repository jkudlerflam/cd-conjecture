# C123 progress

2026-10-29.

[Work record] Read the project memory, all background updates, REPORT sections 9.68 through 9.82, C121 RESULT, PROOF, CHECKS and source in full, and C122 PROGRESS and RESULT. All writes and Python execution stay in this folder. The local STATE.md will carry the handoff because the task prohibits changes to the project memory.

[Proved, inherited consequence] C121 gives R_L-3N log d <= G <= R_L+2N log d. Combining its upper bound with C118 proves the requested bounded-cost lemma for XX for every fixed finite spacing a. There is no finite largest admissible a in this subclass. This consequence does not compute the actual G.

[In progress] Implement exact finite-grid occupation-insertion determinants for multi-pinch XX, including both Jordan-Wigner sectors and singular parity traces at the two zero modes. Calibrate against full small-ring ED before large rings. Implement charge-block full-spectrum ED for interacting models within modest local resource limits.

[Open] A universal interacting reconnection bound is not established. C122 rules out full operator-order comparison at the requested scale and explains why a gap alone does not justify a uniform spectral truncation.

[Proved] Completed the XX theorem with b_0=2 log 2 and b_1=2 log 2+(10pi/3)a, and proved that admissibility at one spacing implies admissibility at every finite spacing with changed constants. Derived exact multi-pinch occupation determinants and the zero-mode bordered determinant. At fixed a,N, their large-L limit gives G=R_L-N log 2+O(L^(-2)).

[Numerical] Completed 51 direct determinant calibrations, 156 XX cases through L=8192, and 192 full-spectrum ED cases through 12 qubit sites and 8 SU(3) sites. The actual G at L=8192,beta=4L,N=1 is 2.372985925; the prior 0.986691567 uses the unpinched denominator. All mathematical assertions pass. SU(3) L=4 is explicitly excluded from vacuum-sequence inference.

[Proved] Derived exact positive Hilbert-Schmidt temporal transfers and a sufficient residual in terms of their spectral-radius ratio and denominator boundary overlap. Twelve independent algebraic checks pass. The spectral-rate estimate is necessary for the all-N target, while the boundary estimate is a sufficient stronger condition.

[Open] No interacting uniform estimate or all-quasi-free channel theorem is proved. The requested larger interacting computations were not run. Network commands are forbidden, so no Typhon submission was attempted. Wrote RESULT.md, PROOF.md, CHECKS.md and task-local STATE.md with these limits explicit.

[Work record] Final validation passed Python syntax, prose restrictions, all saved-result counts, the explicit XX reconnection ceiling on every saved case, all interacting R/G inequalities and history-sum checks, and all transfer bounds. No failed mathematical check remains. The general interacting residual is still open.
