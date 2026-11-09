# C181 progress log

This run was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

- Read specs/C181_C185_counterexample_push.md, _qubit_c_context.md, STATUS.md, REPORT 9.137 to 9.140, specs/C175_C180_fibonacci_swarm.md, RESULT of C173, C175, C179, C180, C175 PROOF, C164 RESULT (mass imbalance), C176 RESULT and code, and the collaborator's contrib/cusp_binding_certificate_20261109, cusp_scattering_20261109, cusp_validation_20261109 (read only).
- Copied unmodified: fibdmrg.py, run.py, analyze.py (C176), scat.py, nbody_box.py (C175).
- Few-body predictions F1 to F3 registered in RESULT.md, then fewbody.py run: fewbody_w2_R100.out (|w| = 2, nine φ), fewbody_w1_R400.out (|w| = 1, background, slow at local load average 180 to 350).
- Design pilot local/pilot_w2_L48.out (E(N) at |w| = 2, φ = π/2, L = 48), used only to choose densities. Smoke test local/smoke (complex run.py at L = 16).
- Many-body predictions M1 to M5 and the first-stage criterion registered in RESULT.md before submission.
- Typhon job 17604449 submitted 2026-11-09 (one 48-core node, 36 h limit, 134 tasks, 24 at a time, 2 threads each). Scratch /scratch/lustre/jkudlerflam/c181/, outputs out/<tag>.jsonl (one line per sweep, measurement on the last sweep), logs/<tag>.log. Task list typhon/tasks.txt from make_cfg.py.
- Amendment (job pending, no output): |w| = 2 main density changed from 1/4 to 9/32 (n = 1/4 is commensurate with period 4). Configs on scratch replaced in place, job not resubmitted. Registration amended in RESULT.md. Design pilot killed (local load).
- M6 (FFLO pair momentum) registered in RESULT.md, job still pending, no output.
- Ring magnetization-curve predictions R1 to R4 registered; ring_mag.py runs launched locally (two background chains).
- fewbody |w| = 1 stopped after φ = 1.56 (load). Local DMRG previews launched at |w| = 2, n = 9/32: φ = 1.50 at L = 32, 64, 96 and φ = π/2 at L = 64 (χ to 256), local/dmrg/. These read against M3 and M4 as previews only.
- Job 17604449 cancelled while pending (no output) and resubmitted as 17604471 with --mem=96G instead of 300G (same tasks, same scratch). Estimated start still 2026-11-13 (Priority).
- Local previews at |w| = 1, φ = 1.555, n = 3/16 launched: L = 64 (N = 11, 12, 13) and L = 96 (N = 18).
- Local previews added: |w| = 1, n = 3/16, L = 64, φ = 1.53 and 1.50 (to see the imbalance set in on the open chain).
- Local: |w| = 1, φ = 1.53, L = 96, N = 18, and the c = 1 and c = 2 controls of the Typhon list at L = 64 and 96, n = 3/16 (calibration of the 64 -> 96 estimator).
- Local previews done: |w| = 2 (φ = π/2, 1.50; L = 32, 64, 96), |w| = 1 (φ = 1.555 at L = 64 with N ± 1 and L = 96; φ = 1.53 at L = 64, 96; φ = 1.50 at L = 64), controls at L = 64, 96. summarize.py and relwin.py (L-proportional window, chosen after the controls were seen, flagged as unregistered). Results in PROOF 5.2.
- Ring runs: ringmag_w2_L28/L32, ringmag_w2V0_L32, ringmag_w1_L30/L36 complete; ringmag_w2_L36 partial at write-up (still running locally).
- RESULT.md and PROOF.md written. Typhon job 17604471 still pending.
- Coordinator request: backfill job 17604615 (16 cores, 11 h 55 min, 16 runs: |w| = 1, φ = 1.53 and 1.555, L = 64 and 128 with N ± 1, ctl1 and ctl2 at L = 64 and 128, χ to 384, one thread each), typhon_bf/. Job 17604471 trimmed in place on scratch to 124 tasks (removed w1n316_p1.5550 and ctl1/ctl2). Both recorded in RESULT.md.
