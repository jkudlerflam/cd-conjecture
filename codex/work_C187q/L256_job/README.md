Coordinator job (Claude Opus 5.5, 2026-11-09), Typhon 17604910 (c187b): 32 cores, 12 h, scratch /scratch/lustre/jkudlerflam/c187b/.

The L = 256 stage of C181's partially polarized point (REPORT §9.150). It uses C173's effective model at bare couplings, with |w| = 1, n = 3/16, and φ = 1.53 and 1.555. Each φ is run at N = 47, 48 and 49, alongside the c = 1 and c = 2 controls at N = 48. The code is C176's fibdmrg.py and run.py, copied unmodified from codex/work_C181. The bond dimension χ rises to 768, and measurement is done on the last sweep only.

Registered before any output. Two-size c is taken from L = 128 to 256, pairing C181's backfill L = 128 runs (codex/work_C181/typhon_out_bf/out_bf/bf_*) with these runs. It uses windows of W = 1, 6 and 12 cuts, as in C181's summarize.py.
- Calibration. The c = 1 control must read within 0.15 of 1, and the c = 2 control within 0.25 of 2. If either fails, the run is uncalibrated and nothing below applies.
- Pass, meaning the reading is consistent with c = 2. At both φ, c_two(128 → 256) is at least 1.6 in at least two of the three windows. The truncation at L = 256 must be below 1e-7. L × the parity gap at L = 256 must not exceed 1.5 times its L = 128 value (1.53 at φ = 1.53, 1.50 at φ = 1.555).
- Fail, meaning the reading is a crossover toward c = 1. At either φ, c_two(128 → 256) is at most 1.3 in at least two windows, or L × the parity gap more than doubles from 128 to 256.
- Anything else is inconclusive.

Collect from codex/work_C187q/L256_job/:

    mkdir -p typhon_out && scp -o LogLevel=ERROR 'typhon:/scratch/lustre/jkudlerflam/c187b/out/*.jsonl' typhon_out/ && scp -r -o LogLevel=ERROR typhon:/scratch/lustre/jkudlerflam/c187b/logs typhon_out/

Do not copy out/ckpt.

Amendment (2026-11-09). This was registered before any c estimate was computed. Some log lines had been seen in a progress check: candidate runs at χ ≤ 64, and the ctl1 run's S_mid. cfg/ctl1_L256_N48.json has the same error as the desktop L = 192 config: "t2": 0.0, where C181's L = 128 ctl1 keeps the default t2 = −1. The c = 1 calibration uses ctl1fix_L256_N48 (the default t2), run in a separate small job, and ctl1_L256_N48 is reported as a config error.

Result (2026-11-09, REPORT §9.166). All runs completed the schedule at χ = 768. Calibrated: ctl1fix reads 1.040, 1.002, 1.007 and ctl2 reads 2.085, 1.920, 2.001. φ = 1.53 reads c_two(128 → 256) = 1.743, 1.720, 1.751, and φ = 1.555 reads 1.773, 1.842, 1.851. Truncation is 3.5e-8 to 4.0e-8. L × the parity gap at L = 256 is −0.18 and 0.49. Verdict: pass. The entropies are still rising with χ at 768 (by 0.011 from 512), so the reading is not converged to better than about 0.1 in c. Numbers: `python3 -E -s pair256.py`, output in pair256.out.
