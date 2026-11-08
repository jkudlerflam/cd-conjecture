This work was done by a Claude Opus agent (2026-11-08), not by Codex. Astra was not used (out of credits).

# C164 result: lattice search for accidental two-valley soliton onsets

Derivations and numbers are in PROOF.md, the log in PROGRESS.md. Scripts in this folder: clat.py (lattice classical spin chain), calib_classical.py, lsw_check.py, solscreen.py (classical screen), qband.py, qscreen.py, qopt.py (adversarial quantum two-valley search), qtune.py (continuation to onset), twobody.py, onset_dmrg.py (past-onset DMRG, version 2), controls.py, masscheck.py, kpik_check.py, test_corr.py, make_endpoints.py, summarize.py. Copied unmodified: est.py, dmrg.py, dmrg_c.py (C155), vumps.py, fes.py (C153), fesx.py, fesfit.py (C158), kfun.py, scatlen.py (C162), soliton_valley.py, soliton_pair.py (C161). dmrg_d.py is dmrg.py with the local dimension inferred (for the d = 4 control). Local pilot data are in local/.

## Verdict

Null so far. No explicit c = 2 bond. The local search covers 188 classical bonds, 610 quantum ring screens and 16 optimizer restarts (about 2400 objective evaluations) with 4 continuations, and it finds accidental two-valley quantum bands in N but no two-valley onset with separated, non-binding valleys. The decisive large search and the past-onset DMRG test are in Typhon job 17601039, which was pending at write-up. Three structural results narrow (A_sol): the mass-imbalance correction has a definite immiscible sign, LSM reduces T to accidental cases, and a new protected stratum (k -> pi - k) turns out to be C152's two-Majorana surface.

## Calibrations

- [numerical] Classical lattice code: C161's tilt eigenvalue is reproduced, -0.0671 at eps = 0.05 (-0.0645 at 0.1, -0.0677 at 0.025) against -0.067. C162's K binding is reproduced, U_same(T = 10)/eps = -0.0084 at eps = 0.05 against -0.0088, with U_alt = +0.0158 against +0.0169 and E_sol + U_same < 0. LSW matches the exact one-magnon band to 1e-14.
- [numerical] DMRG c at the c = 1 point XX: 0.977 at L = 32 locally. L = 64, 128 and 256, the c = 2 control (two decoupled XX chains as d = 4 super-sites, a two-site cell with range-3 couplings, outside B_2) at L = 32 and 64, and FES for both are in the job (ctrl/).
- [numerical] Correlators and the helicity-resolved S_+(q) checked against ED to 6e-14. Ring ED tables exact to 1e-14 (est.py).

## Main findings

1. Classical screen (Task 1) [numerical]. Genuine soliton minima (Hessian-tested) occur in 39 of 121 stable backgrounds, 28 of them below the magnon gap. No bond had a second soliton valley. At eps = 0.3 to 0.5 the K solitons are 10 to 28 sites wide, and the 14-site quantum band puts the soliton well above the classical value (0.72 against 0.31 at eps = 0.5). The classical soliton band is a poor guide at eps of order one [heuristic].
2. Quantum bands (Task 2) [numerical]. Smooth, isolated, near-degenerate two-valley bands above uniform k = 0 backgrounds exist in N. The best has k_A = 0.43, k_B = 2.15, E = 0.273 and 0.275, bandwidth 1.04. Naive optimization finds many more, but rechecks at L = 14 to 18 show them to be period-2 or level-crossing artifacts. They are excluded by a smoothness criterion (high-harmonic weight < 2 percent). Continuation toward onset merges or flattens the valleys. The best endpoint reached gap/bw = 0.066 with a flat bottom (k = 0.45 to 1.8), and its 14-site two-particle levels sit above the free floor (no like-valley bound state), too small a ring to extract J_z.
3. Mass imbalance [numerical, derived]. For impenetrable two-species gases, E(alternating) - E(segregated) > 0, proportional to eta^2 n^2. Each A-B interface costs about 0.3 eta^2 hbar^2 n^2/M. This is the leading dilute term (n^2 against n^3 for scattering lengths), so unequal valley masses force immiscibility as n -> 0. A miscible window needs n Delta a > 0.1 eta^2 roughly.
4. Symmetry [proved or literature, numerical checks]. By LSM with time reversal, T has no unique gapped symmetric background, so T two-valley bands are accidental (ring artifacts with +-k0 valleys confirmed non-uniform). D = 0 bonds show no incommensurate minima (195 samples). The Pi-symmetric surface with traceless transverse exchange has an exact symmetry U = I R_{pi/2}(z) prod_odd R_pi(z) giving E(k) = E(pi - k) on one-particle states (7e-14). This is C152's cos phi = 0 surface, the two-Majorana case, c <= 1 there.
5. Kinematics [derived]. Two emergent valley charges need k_A, k_B not in {0, pi} and k_A + k_B != 0 mod 2 pi. The C161 K valleys (pi -+ delta) sum to 2 pi, so at eps of order one only N_+ - N_- is protected there [heuristic].

## What each estimate probes (Task 4)

The quantum particles found are lattice scale (at most about 3 sites, from band smoothness). Past onset at depth 0.1 bandwidth the density is 0.1 to 0.2, the pseudospin scale is 3e-3 to 3e-2 (scattering lengths) against about 6e-4 from mass imbalance at eta = 0.3, and L = 64 to 128 holds 6 to 26 particles. The bisection in the splitting resolves about 1.6e-3 of the bandwidth. A miscible window of width 3e-3 or more would be visible, and one below about 1e-3 would not. The blow-up regime (eps <= 0.1) stays out of reach (C162: pseudospin scale about 1e-7).

## Typhon job 17601039 (PENDING, Priority, estimated start 2026-11-08 12:43)

Directory /scratch/lustre/jkudlerflam/c164/, script typhon/c164.slurm, task lists stage1.txt and controls.txt. onset_dmrg.py was replaced by version 2 on scratch after submission (copy in typhon/onset_dmrg_v2_uploaded.py, md5 01eed6aa...). Stage 3 reads it at run time.
- controls (7 cores): ctrl/dmrg_xx_{64,128,256}.out, ctrl/dmrg_twoxx_{32,64}.out, ctrl/fes_xx.json, ctrl/fes_twoxx.json, ctrl/fesfit.out. Pass: XX within 0.05 of 1, two-XX within 0.1 of 2, FES slopes about 1 and 2.
- stage 1 (3 h): qopt/q_*.jsonl (31 tasks, hF, dm, S^8) and cl/cl_*.jsonl (10 classical screens: N at eps 0.3, 0.5, 0.8, 1.2, T at 0.5, 1.0, dm x2, hF, S^8).
- stage 2 (2 h): qtune/t_*.jsonl, the continuation of every result with O < 0.3.
- stage 3 (4 h): endpoints.jsonl (41 best by gap/bw), onset/e_*.json with coarse scan, bisection and L = 64, 128 at lam*.
- summary.txt (summarize.py). Collect with `scp -r typhon:/scratch/lustre/jkudlerflam/c164/{summary.txt,onset,ctrl,endpoints.jsonl} .`

Pre-fixed c = 2 criteria: entanglement fits > 1.5 not falling from L = 128 to 256, truncation < 1e-7, two incommensurate S_+ peaks at k_A and k_B with comparable weight at lam*, FES slope > 1.5 at chi >= 256. Stage 3 only promotes points (c_prof(64) > 1.2 goes to L = 128). L = 256 and FES on a survivor are the follow-up.

## Residual assumption (A_sol') [assumed]

No accidental two-valley onset in N (or in T, which LSM reduces to the same accidental type) at eps of order one combines separated, non-binding valleys, a density with n Delta a > about 0.1 eta^2 (or near-equal masses), and J_z > -|J_perp|. The dilute limit with unequal masses is now excluded. The protected staggered case is covered by C152.

## Next steps

(i) Read summary.txt. For any endpoint whose L = 128 fit exceeds 1.3 with two S_+ peaks, run L = 256 and FES with fesx.py gen (two-site cell, chi to 256). (ii) For candidates that are clean, extract scattering lengths from DMRG with a few particles (or an MPS excitation ansatz), to get J_z and J_perp directly. (iii) The quantum soliton band at eps of 0.2 to 0.4 needs rings of 30 sites or more, beyond ED.
