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
