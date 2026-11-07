# C124 handoff

2026-10-30.

[Open] The general noncommuting stoquastic and interacting reconnection lemma is not proved. No violation of it or of A<=C M f(d) was found.

[Proved] Commuting adjacent bonds admit an on-site direct sum of incoming, outgoing and multiplicity factors. In the adapted basis the positive reflected histories obey q_s=B_ss and q_D(s,t)=B_st B_ts with entrywise nonnegative B. Hence G<=0 for every spacing. In an arbitrary basis, the sector transfer matrix gives R<=log r<=log d and G<=log r+2N log d. The class includes entangled bond eigenstates but has A=0 under the stated infrared assumption.

[Proved] The boundary-overlap condition proposed as sufficient in C123 is not uniform. For a diagonal on-site qutrit with levels 0,epsilon,1, epsilon=L^(-2), M=1/2, tau=2aL, its C_L is (2L-1)log 2-2a+O(1/L), whereas G is between -log 3 and zero. Both temporal spectral radii equal one and w_D=1. Tensoring XX with a slow on-site qubit gives the same obstruction in a critical stoquastic d=4 family, while its G remains bounded by the XX value. Preserve the full thermal contraction; a separate top-eigenspace boundary estimate loses a cancellation.

[Open exact remaining issue] At fixed tau, G_N=G_1+sum Gamma_n, where Gamma_n=log(S_(n+1)/S_n)-log(D_(n+1)/D_n). A uniform G_1 bound and a uniform increment bound in a/2<=tau M/L<=a would suffice. Their Rayleigh quotient form retains the evolved endpoint vectors. Neither has been proved in the noncommuting class. Bounds on each increment are sufficient, not necessary.

[Numerical] All new applicable checks pass. Large XX reaches L=8192, doubled length 16384, N=6. Full-spectrum ED reaches doubled length 12 for qubits and 6 for SU(3), N=8. Commuting separator residual <=3.50e-14; free-fermion dense-history calibration <=8.92e-15; direct temporal-transfer G versus ED <=1.87e-13. The noncommuting extension G<=0 fails in all requested families at sampled parameters.

[File map] RESULT.md gives the outcome. PROOF.md contains full derivations, especially Sections 4-7. CHECKS.md gives scope and reproduction. checks.py, extra_checks.py, interval_checks.py, transfer_checks.py perform the checks; xx_kernels.py and ed_kernels.py are locally copied and calibrated C123 kernels. Scalar outputs and run transcripts are listed in CHECKS.md.

[Work record] Only work_C124 was modified. Project-level STATE.md was left untouched because the explicit task rule forbids outside writes. The excluded independent attempt was not read. No subagents, network shell commands, cluster jobs or external writes were used.
