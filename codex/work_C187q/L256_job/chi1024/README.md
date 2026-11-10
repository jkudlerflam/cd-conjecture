# C187d: continuation of the L = 256 stage to χ = 1024 (registered 2026-11-09, before any output)

Coordinator: Claude Fable 5. Jonah approved this continuation on 2026-11-09. Context: REPORT §9.166 and `../README.md`.

## Why

The L = 256 stage of the partially polarized points passed its registered rule at c_two(128 → 256) = 1.72 to 1.85, but the entropies were still rising with χ (by 0.011 from χ = 512 to 768 at L = 256), and the L = 128 partners stopped at χ = 384. More χ at L = 256 raises the reading and more χ at L = 128 lowers it. At matched χ = 384 the reading was 1.47 to 1.61.

## Runs

All runs resume from the existing checkpoints with the unmodified code (`run.py`, `fibdmrg.py`, `analyze.py` of the C187b and C181 backfill jobs). The original outputs are left untouched, and the continuation works on copies in a separate scratch folder.

- The six candidate runs at L = 256 (φ = 1.53 and 1.555, N = 47, 48, 49) get three more sweeps at χ = 1024.
- The c = 2 control ctl2 at L = 256, N = 48 gets the same three sweeps.
- The L = 128, N = 24 partners of both candidates and of ctl2 (C181 backfill tags bf_*) get sweeps at χ = 512, 512, 768, 768, 1024, 1024.

Configs are in `cfg/`. The c = 1 control is not continued. It was converged at χ = 768.

## Readings

S̄_W(L, χ) is the entropy averaged over the cuts L/2 − W to L/2 + W at the last sweep of the stage with bond dimension χ, and c_two = 6 [S̄_W(256) − S̄_W(128)] / ln 2, as in `../pair256.py`, for W = 1, 6, 12.

- R1, matched. Both sizes at χ = 1024.
- R2, extrapolated. For each size S̄_W is extrapolated linearly in the discarded weight of the last sweep of each stage, through the stages χ = 512, 768 and 1024. The uncertainty is the difference from the two-point extrapolation through 768 and 1024.
- Convergence. A size counts as converged in χ if |S̄₆(1024) − S̄₆(768)| is at most 0.005 at L = 256 and at most 0.002 at L = 128. If it is not, R2 with its uncertainty is the number to quote and the reading is called not converged.
- Calibration. R2 counts only if the control ctl2 gives R2 within 0.25 of 2 under the same procedure.
- Parity gap. L·Δ_P at L = 256 from N = 47, 48, 49 at χ = 1024, reported against the L = 128 values (1.53 and 1.50), with the earlier criterion (not above 1.5 times the L = 128 value).

## Verdict rule

The thresholds are those of the C187b registration.

- E1, the reading stays at c ≈ 2. R1 is at least 1.6 in at least two of the three windows at both φ, and R2 is at least 1.75 in at least two of the three windows at both φ.
- E2, fail. R1 is at most 1.3 in two windows at either φ.
- E3, anything else. In particular R1 ≥ 1.6 with R2 below 1.75 is reported as a pass of the original rule with a converged value below the control, with the numbers.

## Where it runs

Scratch folder /scratch/lustre/jkudlerflam/c187d/ on Typhon. The payload runs as two extra steps inside running cd allocations with idle cores (`body.sh` with `tasks_A.txt` and `tasks_B.txt`), and `c187d.slurm` is the queued fallback. `runone.sh` holds a lock with a heartbeat per tag, so a tag never runs twice at once, and a finished tag is skipped because every sweep is already recorded.

Collection, from this folder (never copy out/ckpt):

    mkdir -p typhon_out/logs
    scp -o LogLevel=ERROR 'typhon:/scratch/lustre/jkudlerflam/c187d/out/*.jsonl' typhon_out/
    scp -o LogLevel=ERROR 'typhon:/scratch/lustre/jkudlerflam/c187d/logs/*.log' typhon_out/logs/

A run is complete when its jsonl has 28 lines (L = 256) or 27 lines (L = 128).

## Launch record (2026-11-09 22:52 EST)

- The code in c187b and c181 is identical (same md5 for run.py, fibdmrg.py, analyze.py), so one copy serves both sets of checkpoints.
- Step A, the six candidate runs with 4 threads each, is step 17603971.1 (inside the c177 allocation, node38, 2 d 10 h left at launch).
- Step B, ctl2 at L = 256 (4 threads) and the three L = 128 partners (3 threads each), is step 17604038.1 (inside the c175 allocation, node12, 10 h 48 min left at launch).
- Fallback job 17608223 (c187d, pending). It runs `tasks_all.txt`, and skips every tag that is finished or locked.
- Resume check at 22:54 EST: bf_ctl2_L128_N24 wrote sweep 21 (χ = 512, 112 s) after its 21 existing lines, so the runs continue and do not restart.
- Memory. The earlier job peaked at 15 GB for eight runs at χ = 768, so about 27 GB is expected here, inside allocations of 200 GB that use 32 GB and 8 GB.
- Expected duration: 3 to 4.5 h for step A (three χ = 1024 sweeps at L = 256, plus the measurement), and less for step B.

Check, in one ssh call:

    squeue -h -s -u jkudlerflam -j 17603971,17604038 -o '%i %j %M' | grep '\.1 '; wc -l /scratch/lustre/jkudlerflam/c187d/out/*.jsonl; cat /scratch/lustre/jkudlerflam/c187d/logs/srun_?.out

If a step dies with its donor before its tags are complete, the fallback job picks the tags up, or relaunch the same body in another running c1xx allocation with idle cores (`srun --jobid=<donor> --overlap -n1 -c48 bash body.sh tasks_A.txt 6`). Never two copies of a body on purpose, although the locks protect against it.
