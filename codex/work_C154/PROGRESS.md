# C154 progress log (Claude Opus agent, 2026-11-07)

This run was done by a Claude Opus agent (2026-11-07), not by Codex.

## Setup

- Read specs/C154_generic_bonds.md, specs/_qubit_c_context.md, STATUS items 1, 4, 5, contrib/d2_gap_reduction/README.md, C141, C143, C144 RESULT, C144 PROOF section 1.
- Copied C144's validated DMRG (dmrg.py, xzlib.py) and C155's momentum-resolved ring ED (est.py, saved as est_c155.py, unmodified) into this folder.
- Typhon check at start (18:2x local): C144 jobs 17580192 and 17580432 are still RUNNING (5.5 h and 4.7 h). No yx_*, zy_*, t_L128+ or u_* jsonl exists yet in /scratch/lustre/jkudlerflam/c144/. Only coneBx_L128 finished (gapped, c_fit 0.0). The L = 64 tuned-line files exist.

## Key early observation (to be proved in PROOF section 3)

- XZ + eps YX + eps' ZY = sum (alpha XZ + beta YX + gamma ZY) has J = D P with D = diag(alpha, beta, gamma) and P the cyclic permutation of axes (a 2 pi/3 rotation about (1,1,1)). The period-3 site-dependent rotation sigma_j = P^{-j} tau_j maps it to an XYZ chain whose bond couplings cycle through the cyclic shifts of (alpha, beta, gamma). So the YX and ZY directions at XZ are not symmetry-free: they lie in a 2-dimensional stratum with a twisted (period-3) Klein symmetry Z2 x Z2 and Kramers time reversal, and (alpha, beta, gamma) = (1,1,1) is the antiferromagnetic Heisenberg chain. This explains C144's period-3 commensuration and exact triplets.

## Numerics, round 1

- cyc_check.py [proved, checked to 1e-14]: spectra of the uniform ring h = a XZ + b YX + g ZY equal those of the period-3 XYZ ring with bond couplings (a,b,g), (g,a,b), (b,g,a) (L = 6, 9; four parameter points); the three twisted Klein generators prod_j (P^{-j} e_n).sigma_j and Kramers Theta = prod(iY)K commute with H exactly.
- cyc_spec.py (momentum ED, L = 12, 18; cyc_spec.jsonl, cyc_spec_run1.out): along the YX line (1, eps, 0), eps = 0.05 to 1, through (1,1,0.5) to Heisenberg (1,1,1), every point shows the same tower: unique ground state at momentum 0 (L = 12) or pi (L = 18), an exact triplet at relative momenta pi + {0, +-2 pi/3} (x near 1/2), a singlet at relative momentum pi just above it, and an exact current sextet at relative momenta +-2 pi/L + {0, +-2 pi/3} (x = 1). ZY (1, 0, 0.1) is identical to YX (1, 0.1, 0) (mirror). At eps = 0.05 and L = 18 the singlet (0.840) lies below the triplet (0.915): possible marginally relevant flow (dimerization) at small eps, or A-sector crossover. Needs L = 24 and DMRG.
- Typhon job 17594697 (c154_1, partition typhon, 48 cores, packed): ED at L = 24 in sectors m = 0, 1, 4, 12 for eight cyclic points (ed24.jsonl), and open-chain DMRG for (1, eps, 0) with eps = 0.05, 0.1, 0.2 at L = 48, 96, 144, (1,1,0) and Heisenberg (1,1,1) at L = 48, 96 (dm_*.jsonl). Outputs in /scratch/lustre/jkudlerflam/c154/.

## Coordinator messages (received mid-run)

1. Add the phase-diagram topology angle: characterize the c = 1/2 and c = 1 loci of symmetry-free bonds and show they never meet in the pattern a c > 1 fixed point requires. 2. Treat explicitly (a) isolated c > 1 gap closings with no lower-c loci nearby and (b) a stable c = 2 floating phase, and say whether a NN qubit chain can support two independent incommensurate modes. Both are addressed in PROOF sections 5 and 7.

## Prior work found and used (read only)

codex/work_C67, work_C68 (strata T, K, N; floating-region argument; operator counts), contrib/d2_floating_ising, d2_multicritical, d2_two_fluid, d2_classification. Not repeated; cited.

## Round 2

- Theorem 1 (PROOF 1.3): complete classification of site-dependent on-site symmetries for distinct singular values. Only new twisted stratum: the cyclic one (period-3 Klein, B = 0).
- spiral.py: polar decomposition and U(1) parent; point 1282 has harmonic fraction 0.61, parent XXZ with Delta = -0.414, h/J_perp = 0.22 (a Luttinger liquid).
- lsw_check.py: linear-spin-wave band (5.1) agrees with exact L = 12 ring excitations at strong field, error scaling as 1/|B| (0.033 -> 0.008 for |B| 8 -> 32), with C155's momentum sign convention (k -> -k).
- L = 24 ED (Typhon 17594697, ed24.jsonl): on the YX line the singlet/triplet ratio r_s = x_s/x_t at L = 12, 18, 24 is 1.185, 0.918, 0.787 (eps = 0.05), 1.264, 1.090, 1.020 (eps = 0.1), 1.386, 1.288, 1.251 (eps = 0.2), 1.471, 1.396, 1.362 (eps = 0.3), 1.714, 1.626, 1.574 (Heisenberg). So at small eps the isotropic marginal coupling is negative (singlet below triplet and falling): spontaneous dimerization with an exponentially small gap; for eps >= 0.2 it is positive (SU(2)_1, c = 1). Casimir c (v from the current level) at (18, 24): 1.79 (0.05, A-sector contamination), 1.22 (0.1), 1.08 (0.2), 1.045 (0.3), 1.021 (0.5), 1.010 (1,1,0), 1.004 (Heisenberg).
- Local machine load average is above 100 (other agents), so the grid scan and the 1282 path were moved to Typhon: job 17594764 (c154_2: L = 24 ED at eps = 0.07, 0.085, 0.125, 0.175, 0.25 and three off-line points; DMRG eps = 0.15, 0.3, and eps = 0.05 at L = 192) and job 17594765 (c154_3: the fundamental-triangle grid, cyc_grid.jsonl; the spiral path of 1282 at L = 32, 64, 128, path1282.jsonl). Outputs in /scratch/lustre/jkudlerflam/c154/.

## Round 3 (close-out)

- L = 12, 18 crossing data at eps = 0.07 to 0.25 (cyc_spec_run2.out): r_s = 1 crossing at eps = 0.07 (L = 18), about 0.095 (L = 24); eps_c estimated between 0.1 and 0.15 (heuristic extrapolation).
- Local 20-point grid of the cyclic triangle (cyc_grid_local*.out): SU(2)_1 region for gamma/beta >= -0.5 (c_S 0.98 to 1.07), ferromagnet-like ordered region near the XZ-FM edge; no third phase.
- Typhon DMRG L = 48 finished for eps = 0.05, 0.1, 0.2, (1,1,0), Heisenberg (typhon_out/dm_L48.jsonl, about 33 min each). L = 96, 144 still running.
- Local 1282 path pilot: L = 32 complete (Sbar smooth in lam), L = 64 only lam = 0 (c_mid 0.99) and 0.25 (c_mid 0.51, unexplained) before the machine load made local runs impractical.
- PROOF.md sections 1 to 8 and RESULT.md written.
- Pending at close: 17594697 (running, DMRG L = 96, 144), 17594764 and 17594765 (queued). Outputs in /scratch/lustre/jkudlerflam/c154/. C144's 17580192 and 17580432 still running, no YX/ZY outputs yet.
- 2026-11-08: L = 96 DMRG collected (PROOF 9). Job 17594764 started; 17594765 still queued; L = 144 still running.
