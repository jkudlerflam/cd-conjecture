# C98 local state

2026-10-29. Investigation complete; universal logarithmic Casimir bound remains open.

Read RESULT.md for the verdict and PROOF.md for derivations. All new files are restricted to this folder. The project-level STATE.md was not modified because the task explicitly prohibits writes outside work_C98.

[Proved] D_L=E(2L)-2E(L)=pi A/(4L)+o(1/L) on a vacuum-compatible sequence. Elementary product gluing gives D_L<=4M and |E(L)-e_0 L|<=2M, not an inverse-length bound. Exact dyadic telescoping is in PROOF.md Section 2.

[Proved failures of stronger claims] Sharp no-remainder finite-L bound fails at Heisenberg L=4. Raw XX closing cost has the nonzero limit -(1/2-1/pi). XX half-ring Schmidt rank is 2^(L/2). A distributed partial-copy swap has a positive bulk cost, 0.704986412078 sin^2(2theta) per doubled Heisenberg bond and 0.447316023917 sin^2(2theta) for XX. Reflection concavity plus bounded deficit does not fix the 1/L coefficient; the abstract sequence e_0 L-B xi/(L+xi) shows the inference gap, without constituting a Hamiltonian counterexample.

[Proved reduction, open estimate] Two distant bond reconnections V_x,V_y define H(s,t). Exact identity D_L=(1/2)lim_beta integral_[0,1]^2 chi_xy^(beta), with chi the connected Kubo susceptibility. A positive-part integrated response bounded asymptotically by pi^2 M log_2(d)/(2L) proves the sharp conjecture. This response bound is unproved and has not been checked in the thermodynamic limit. Its signed version would only restate the target. Average local-marginal extension error <=(pi^2/8)log_2(d)/L^2 is another sufficient, unproved criterion.

[Inherited] C14/C35 supply the symmetric SU(N) Bethe/WZW formulas and the proof of their logarithmic comparison. Do not repeat those derivations. New checks are small-ring calibrations, not a new continuum proof.

[Numerical] Run OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 checks.py and the same prefix with python3 response_check.py. Both passed. CHECKS.md and RESPONSE_CHECK.md are generated reports. The initial command `python` was unavailable; `python3` worked. No other check failed. No large computation or cluster submission.
