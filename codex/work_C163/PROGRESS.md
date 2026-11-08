# C163 progress log

This run is done by a Claude Opus agent (2026-11-08), not by Codex. Astra was not used (out of credits).

## 2026-11-08

- Read the spec, the shared context, C159 RESULT/PROOF, C161 REFEREE (all sections), C160 RESULT, C154 RESULT, C152 RESULT.
- C159's Typhon job 17597033 has COMPLETED (1 h 51 min, ended 2026-11-08 08:37). summary.txt copied to c159_summary.txt, refine/*.json to c159_refine/. C155's 17594766 is still running; its result folders are empty so far.
- First reading of the C159 refinement: two T points have open-chain entanglement slopes that rise with L (dist_T_50, hi_T_40), both with parent Delta about +1.3 (easy-axis antiferromagnetic parent, harmonic fraction 0.5 to 0.7), not in the attractive window. Everything else in stage B reads c <= 1 at L = 128 or falls with L.
- Copied the C159 scan records (c159_scan/, 48 files). Wrote xpi_c159.py: the scaled level at relative momentum pi on the L = 16 rings is not a usable K estimator (calibration at harm < 0.15 is off by -0.5 in the median), because the spiral twist L theta puts e^{i vartheta} states near pi. Dropped.
- Wrote kest.py, K_S(L) = (L/4) S_mm(2 pi/L) on rings, m the spiral axis. Calibration on XXZ: exact 1.000 at XX, 1.47/1.50 at Delta = -0.5, 1.87/1.98 at -0.7, 2.82/3.48 at -0.9 (L = 16 / exact; approaches K from below).
- kest_c159.out: the C159 stage-B attractive points (TFM and T, parent Delta -0.64 to -1.18, harm 0.1 to 0.58) have K_S(16) = 1.9 to 7.5, inside the window K > 15/8, and their open-chain DMRG reads c_prof(128) = 0.98 to 1.03 with flat two-size slopes (hi_T_37, K_S about 7, slopes 1.42 then 1.15, still falling). The two rising points dist_T_50 and hi_T_40 have K_S = 0.4 to 0.65: repulsive liquids near the KT value 1/2, outside the window (dist_T_50 is also within 0.2 of the XZ point, s = (1, 0.20, 0.04), a.b = 0.04).
- Wrote c163lib.py (bond constructions, J_to_v with 8-site ring spectra agreeing to 7e-14; -R_z(1) equals h_F on open 6-site chains to 2e-14), path_scan.py, refine_big.py, twisted.py, summarize.py, typhon/tasks.txt (68 tasks), typhon/c163.slurm. A 20 s local test of path_scan.py at L = 16 ran end to end. One L = 48, chi = 48 DMRG takes 13 s locally.
- Submitted Typhon job 17600999 (c163, partition typhon, one node, 48 cores, 8 h). Directory /scratch/lustre/jkudlerflam/c163/. Outputs big/*.json, fes/T40.json, fes/T50.json, twisted/*.json, paths/*.json, summary.txt (python3 summarize.py at the end), logs/c163_17600999.out, logs/tasks.log.
- Read the C159 stage-B outliers: ising_TFM_74 (c_two 1.76) is a Theta-broken ferromagnetic cat (S(1) = ln 2, S_mid flat at 32 and 64, jump at 128 from a domain wall); ising_K_14, 15 have c_prof up to 3.6 from oscillation-dominated profile fits but S_mid slopes 0.66 to 0.78.
- Wrote PROOF.md (mode inventory, Lemma L, Proposition P, Lemma U, Proposition G for the twisted curve, numerics, A_W) and RESULT.md.
- Job 17600999 was PENDING (Priority) at the last successful check; later ssh attempts to Typhon were refused ("Connection closed"). Collect with: scp -r 'typhon:/scratch/lustre/jkudlerflam/c163/{summary.txt,big,fes,twisted,paths}' typhon_out/ and read summary.txt against the criteria in PROOF 4.2.
