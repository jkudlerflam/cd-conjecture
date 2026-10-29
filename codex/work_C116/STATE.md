# C116 handoff

2026-10-29. Work is confined to this folder. The parent project STATE.md was read but not edited because the task prohibits edits outside work_C116.

[Open] The universal crossover-free bound A <= C M f(d) remains unproved. No counterexample to it was constructed.

[Proved] R = 2 log Z(L,beta) - log Z(2L,beta) isolates A through the ordered large-L then large-tau slope, A = (4/pi) lim R/tau, under the full infrared trace assumption and compatible lengths. Exact general bound |R| <= 4 beta M is insufficient. An upper bound K(d)+k M f(d) beta/L would give C = 4k/pi.

[Proved] Constant spatial transfer dimension fails even on XX, without positivity if the exact representation has the correct pressure. Necessary log D_beta >= beta (sqrt(2)-4/pi) - 4 log 2. Same-shape unweighted ratios cancel the Casimir contribution. The corrected 4-to-1 combination retains unbounded reciprocal-velocity terms at fixed d,M,A in a two-XX-sector family. R cancels them.

[Proved under stated infrared assumptions] Massive spectators do not enter the slope. Exact massive Ising check gives 0 <= R_massive <= log 2. Twisted XX shows R need not be nonnegative even at even L; compatible ring lengths matter.

[Numerical] Completed small full spectra for XX/Heisenberg L = 4,6,8 and SU(3) L = 3,6. Full-trace XX recursions agree with ED to 5.7e-14. O(L) checks include paired length 4096 for a two-velocity example. No heavy local work.

[Blocked computation] Two separated SSH attempts to Typhon timed out during banner exchange. No job submitted and no remote writes. run_typhon.slurm is prepared for XX/Heisenberg through L = 12 and SU(3) through L = 9. check_rings.py refuses --cluster outside SLURM.

Files: RESULT.md is the result summary; PROOF.md contains complete derivations and primary references; CHECKS.md states numerical scope and values; PROGRESS.md records work; check_rings.py and check_comparisons.py reproduce checks; checks_small.json and comparisons.json contain the numerical output. No independent referee review was performed.
