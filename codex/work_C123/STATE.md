# C123 local state

2026-10-29.

[Work record] This is the task-local handoff. The background prohibits writes outside work_C123, so the project-level STATE.md was not modified.

[Proved] XX bounded-cost reconnection holds at every fixed a>0, even L>=4, all beta, with b_0=2 log 2 and b_1=2 log 2+(10pi/3)a. It follows from C118 and C121. Independent XX stacks replace log 2 by log d. The separate general quasi-free Casimir bound is not rederived.

[Proved] Existence of the general lemma at one a implies existence at every finite a with changed constants. No finite largest admissible spacing exists under the stated criterion. Optimize the cost divided by a instead.

[Proved] Exact XX histories use occupation determinants on a reflected 2N-time contour, both physical Jordan-Wigner sectors, and a bordered determinant for the two zero modes. The nonselective sums are sums of determinants and are generally not Gaussian. At fixed a,N, G=R_L-N log 2+O_(a,N)(L^(-2)).

[Numerical] 51 direct Gaussian calibrations pass to 1.96e-14. 156 actual-G XX cases reach L=8192 and doubled length 16384, N<=6. At L=8192,beta=4L,N=1, G=2.372985925. The older 0.986691567 was log(S/Z_(2L)), not G. 192 full-spectrum ED cases reach 12 qubit sites and 8 SU(3) sites, N<=8, maximum block 924. The XX ED/Gaussian discrepancy is at most 6.75e-13. SU(3) L=4 is not vacuum-compatible; the only vacuum-compatible SU(3) pair is L=3 versus 6.

[Proved] General temporal transfers give S_N=<u_S,T_S^(N-1)u_S> and Q_D=<u_D,T_D^(N-1)u_D>. With lambda_j the largest eigenvalue and w_D the full top-eigenspace weight, G<=log(||u_S||^2/w_D)+(N-1)log(lambda_S/lambda_D). Both maps retain the full quantum state. Twelve independent algebraic checks pass to 3.25e-14.

[Open] Need uniform spectral-rate and boundary-overlap estimates in PROOF (19), or a weaker comparison of the full contracted spectral measures. C122's order obstruction and the gap/multiplicity obstruction prevent the proposed immediate locality or ground-projection argument. No general interacting bound was proved.

[File map] RESULT.md is the main handoff. PROOF.md contains all derivations and the precise residual. CHECKS.md records scope and reproduction. checks.py implements XX, ed_checks.py implements full symmetry-block ED, transfer_checks.py tests the temporal identity. JSON files contain the generated results. No network commands, external writes, cluster jobs or subagents were used.
