This log was written by Claude Opus 5.5 (agent C168, 2026-11-08), not by Codex. Astra was not used (out of credits).

# C168 progress log

## Reading (done)

specs/C168_island_H.md and specs/_qubit_c_context.md; C154 RESULT and PROOF (1, 2, 4 to 7); C159 RESULT and PROOF; C165 RESULT and PROOF 1, 4; C160 RESULT; C162, C164 RESULT; C155 RESULT and typhon_out/summary.txt; REPORT 9.110 to 9.126. Mid-run: the coordinator's note relaying C174 (codex/work_C174/REDTEAM.md, R1, R5, R6).

## Log

1. Stage-A data. Copied the 12,000 C155 stage-A records from /scratch/lustre/jkudlerflam/c155/results_crit/uniform_*.jsonl to c155_stageA/ (read only). stageA_analyze.py: 636 points have crit < 0.15, 28 of them have c_obj > 1.1 (stageA_flags.json). Found that C155 stage D never refined these 28 (refine.py collects only search and gap-closing outputs), so "every rising reading traced to the XZ Néel cat" does not cover stage A.
2. refine_stageA.py (local, 6 processes): rings 16/20 and open-chain DMRG L = 32, 64 for the 28 flags. Output refineA/pt_*.json, table in refineA_table.out. The machine was shared with other agents (load average 90 to 180), so this ran slowly.
3. Proxy calibration of the ring screen (calib_c2.py, local, 324 samples before I stopped it to reduce load). Range-3 qubit chains with known c (Jordan-Wigner fermions with one or two Fermi seas plus a random symmetry-breaking NN term of norm 0.1). Result in calib_c2_local.out. The C155 rule (crit < 0.15 and c_obj > 1.1) flags only 5 percent of c = 2 proxies, because the crit test fails on small incommensurate rings. The rule c_S(12,16) > 1.5 flags 61 percent of them and 15 percent of the c = 1 proxies.
4. Lemma D checks (imp_mix.py, imp_swap.py, imp_swap2.py; outputs imp_mix.out, imp_swap.out, imp_swap2.out). Two-species hard-core gas on a ring: the ordering energy E(ABAB) - E(AABB) is positive and equals 19.6 eta^2 / L^2 (L = 16 to 32), i.e. 0.31 eta^2 hbar^2 n^2 / M per interface, matching C164. The exchange energy scales as L^-3 (dE_x L^3 = -104, -112, -116, -118, -119, -120 for L = 12 to 32). With eta = 0.5 the ground-state weight on alternating orders falls from 0.33 to 0.15 (L = 12 to 32), while at eta = 0 it stays at 0.517.
5. Typhon job 17602943 (c168, partition typhon, one 48-core node, 8 h limit), submitted 2026-11-08 about 22:55, PENDING (Priority). Scratch /scratch/lustre/jkudlerflam/c168/. Stages: P (1600 proxy chains, families free and int, rings 12/16 and 16/20), S1 (416 stage-A points with c_S(12,16) > 1.5 or a C155 flag, rings 16/20), S2 (open-chain DMRG L = 32, 64, 128 on every S1 point with c_S(16,20) > 1.5 or c_obj(16,20) > 1.1, and on all 28 C155 flags). c168_job.py was updated on scratch before the job started (S2 selection now includes the 28 flags; md5 272596acbc38dc2d65a3547480dad926, same as typhon/c168_job.py). Outputs out/proxy.jsonl, out/s1.jsonl, out/s2_selection.json, out/s2.jsonl, summary.txt, logs/. Collect with
   scp -r 'typhon:/scratch/lustre/jkudlerflam/c168/{summary.txt,out,logs}' typhon_out/
   then python3 typhon/summarize.py typhon_out/out.
6. Theory, PROOF.md. Lemma 1 (no kinematic route at one qubit per site), Lemma O (onsets in N by codimension, any reference state), Lemma D (dilute double onsets are immiscible or alternating when the masses differ), Lemma P (Z_m tongues melt through c <= 1 at weak locking, for every m), Lemma H (full neighbourhood of the AF Heisenberg point in N), the frontier classification (Proposition F) and the residual (R_N). Lead 3 was dropped in its connectivity form after C174 R1; what remains of it is Proposition F, which only sorts frontier points into dynamical events.

7. Coordinator notes mid-run. (a) C174 R1, R5, R6: lead 3 dropped in its connectivity form, P4 added (Lemma P), K_II declared out of scope. (b) C171 R1 to R4: Potts × boson (codimension 2 in N because of the chiral operator), C row and swapped pairs (only from (O3) pairs), composites outside Lemma L (in (R_N)(i)). PROOF 3, 5, 5b, 6, 7 updated.
8. Job 17602943 COMPLETED 2026-11-09 01:17 after 1 h 06 min. Outputs copied to typhon_out/. Proxy efficiencies of F1 and F2 for c = 2 are 0.341 (near-free) and 0.178 (interacting). The C155 rule detects 0.018 and 0.006 of them. S1 found 9 of 416 with c_S(16,20) > 1.5. S2 ran DMRG to L = 128 on 45 points. Two points meet the first stage of the rule, 10193 (L = 128 energy anomalous, likely unconverged) and 450 (c_prof falling from 1.62 to 1.21). The local refinement was stopped (Typhon covers all 28 flags to L = 128, and the energies agree to 1e-6).
9. Job 17603367 (c168b) submitted 2026-11-09 for the second stage: FES (two- and four-site cells, chi to 256) and DMRG reseeds for 10193, 450 and control 10540. PENDING, estimated start 2026-11-12 14:26. Local small-chi VUMPS pilots for 10193 and 450 started (typhon2/loc/). At chi = 8, 10193 has a transfer eigenvalue near -1 (xi of order 1e8, err 0.14), which hints at a period-doubled or cat structure that a two-site cell cannot represent cleanly. The four-site-cell FES in 17603367 tests this.

10. Local VUMPS pilots stopped at chi = 24, unconverged (err 0.08 to 0.2). 10193 sits at e = -0.96827 against the DMRG slope -0.977, with a transfer eigenvalue at -1. 450 sits at -0.87335 against -0.87869, with a pair slope of about 0.8. These are inconclusive, and the second stage stays with job 17603367.

## Not done

- No targeted numerics on the P1 corner or on guest onsets inside liquids (these are the residual).
- K_II is outside scope (assigned to (A6)/(H*) as C174 R6 asks).
- The second stage of the pre-fixed rule for 10193 and 450 (job 17603367) had not run at write-up.
