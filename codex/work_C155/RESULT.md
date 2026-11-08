This run was done by a Claude Opus agent (2026-11-07), not by Codex.

# C155 result: adversarial search for c > 1 at a nearest-neighbour qubit bond

Status. The estimator is built and calibrated, the full search is packaged and submitted to Typhon (job 17594766, still pending at the end of the session), and a local pilot of every stage ran. No candidate with c > 1 survived refinement in the local runs. The Kitaev [110]-field line is not decided at accessible sizes. The log is PROGRESS.md, the derivations are PROOF.md.

## Estimator and calibration

- [numerical] est.py: momentum-resolved ring ED with precomputed sector tables (exact to 1e-14 against dense diagonalization), about 1.5 s per point at rings 12 and 16. Screening value c_obj = min(c_S, c_cas): c_S from the half-ring entropies at two sizes (velocity free), c_cas from the Casimir term divided by a velocity read from the lowest level at relative momentum one. Set to 0 on flat-band manifolds (three or more levels within 1e-3/L).
- [numerical] Calibration (calib_12_16.out). Ising: c_S 0.50 to 0.51, c_obj 0.42 to 0.45. c = 1 (XX, Heisenberg, XXZ, XX + DM, XX in field): c_S 0.85 to 1.05, c_obj 0.75 to 0.995. Excess controls with range-3 terms (outside B_2): two XX chains (c = 2) read 1.90, cluster-Ising (c = 3/2) reads 1.39, two Ising chains (c = 1) read 0.92. So the screen reads low by about 0.1 at c = 3/2 to 2 and by 0 to 0.25 at c = 1, and a genuine c = 3/2 would show as about 1.4. Gapped points read 0 to 0.6 when xi is well below 16, and up to 0.98 when xi exceeds the ring. FM, saturation and flat-band points read 0.
- [numerical] Criticality indicator crit = |z_eff - 1| from the ring gaps (calib_crit_12_16.out): 0.002 to 0.04 at Ising, XX and two XX chains, 0.41 to 0.90 at gapped points, near 1 or above at z = 2 points. It misses cluster-Ising (0.46) because of a sector change between rings 12 and 16.
- [numerical] Refinement: open-chain DMRG (C144's validated code) with Calabrese–Cardy fits plus parity term, smoothed fits, and two-size midpoints. At L = 32, 64 it gives 0.52 (Ising), 0.96 to 0.98 (XX), 0.97 to 0.99 (XXZ Delta = -0.5), an error of about 0.04 (dmrg_calib_32_64.out).

## What the search found (local runs only)

- [numerical] Largest screening values: 4.0 (the objective cap) at a point near the Kitaev [110]-field line, and 2.40 from an unpenalized pilot restart. Both are artifacts. The first reads c_S = -2.83 at rings 16/20. The second reads 0.51 at rings 16/20, and DMRG drops from 1.65 (L = 32) to 0.85 (L = 64) with an incommensurate entanglement oscillation. After the pilot, all objectives outside the XZ neighbourhood carry a penalty 2 max(0, crit - 0.2).
- [numerical] Largest refined value: DMRG profile 1.13 (smoothed 1.06) at L = 64, at the uniform-sample point v = (0.112, 0.266, -0.029, 0.344, 0.087, 0.310, -0.679, 0.438, -0.205). It falls from 1.18 at L = 32, the two-size estimate is 0.80, and rings 16/20 give c_S 0.95. I read this as c = 1, not converged. It is queued for L = 128.
- [numerical] XZ neighbourhood (coordinator request). Along YX, ZY and YX + ZY at delta = 0.3 and 0.45, c is 0.93 to 1.00 at L = 64 and 128, and the best XZ-seeded c-maximization point refines to 1.00. This confirms C144's c of about 1 slow sector at moderate delta. It is not more than one.
- [numerical] Kitaev [110] field (C156 seed). For B = 0.1 to 0.8 along the bisector, ring gaps are 1e-5 to 6e-3, and the entropies are not in a scaling regime: c_S swings between 0.79 and 7.1, single-ring slopes are 2.5 to 3.0. DMRG at L = 64 did not converge within the local time limit. [heuristic] The c of this line is not decided here. Its 1e-3 energy scale implies crossover lengths of order 10^3 sites, beyond every size used.
- [numerical] Volume fraction: of 320 uniform points on S^8, 5.6 percent read gapless (crit < 0.15). Their c_obj is at most 1.00, apart from the one point above. The three gap-first minima are critical hypersurfaces: ball tests at radius 0.03 and 0.08 give gapless fractions 0 to 0.25. No open gapless region with c > 1 was seen.
- Coverage of the local runs: 320 uniform points, 11 c-maximization restarts (40 to 90 evaluations each), 3 gap-first restarts with ball tests, 5 XZ-direction points to L = 128, and 6 Kitaev-line fields at rings to 20. This is a pilot, far smaller than the 103,000-point survey. The adversarial coverage is in the Typhon job.

## Pending on Typhon

Job 17594766, /scratch/lustre/jkudlerflam/c155/. Stage E: 104 XZ-neighbourhood points (26 directions times 4 distances, including 12 inside C144's tuned hyperplane) and 23 Kitaev-line points, DMRG to L = 128 and to 256 where an L = 128 fit exceeds 1.05 (xz/). Stage A: 12,000 uniform points (results_crit/uniform_*). Stage B: 16 gap-first tasks for 4 h (results_crit/gap_*). Stage C: 32 c-maximization tasks for 4 h (results/task_*), a quarter of them seeded at XZ, inside the tuned hyperplane, or on the Kitaev line. Stage D: refinement of every cluster above 1.1, plus the 24 best clusters and the two local flags, by rings 16/20 and DMRG to L = 128 (refine/). The digest goes to summary.txt (summarize.py). Collect it with `scp -r typhon:/scratch/lustre/jkudlerflam/c155/{summary.txt,refine,xz} .`

## Next steps

(i) Read summary.txt and refine/ once the job finishes, and take any surviving candidate to L = 256 with chi of 256 or more. (ii) The Kitaev [110] line needs either an effective model on the near-degenerate manifold (C157's construction) or DMRG at L of order 10^3. Neither ring ED nor L of 256 or less can decide it. (iii) The screen should also be run at rings 12/18, which share the period-3 commensuration that biases the 12/16 pair near XZ.
