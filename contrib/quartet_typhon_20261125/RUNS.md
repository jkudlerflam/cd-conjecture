# Run ledger

All submissions use `python3 tools/typhon/typhon.py submit RUN RECIPE`.
Source bytes, SHA-256 hashes, seed, resources and Git revision are retained
per run. A successful allocation is not a convergence certificate.

| Run | Slurm ID | Recipe | Resources | State / numerical outcome |
|---|---|---|---|---|
| quartet-pilot-20261125a | 17764244 | pilot.job.json | 4 CPU, 16 GB, 1 h | COMPLETED 0:0 in 50:58, MaxRSS 509872 KB; chi256 checkpoint saved at deadline, unconverged; chi384/512 not reached |
| quartet-threads-20261125a | 17764268 | threads.job.json | 8 CPU, 4 GB, 10 min | COMPLETED 0:0 in 77 s, MaxRSS 211648 KB; actual MKL count 1, not a valid 4/8-thread speedup test |
| quartet-c1control-20261125a | 17764277 | control.job.json | 1 CPU, 4 GB, 1 h | COMPLETED 0:0 in 14:18, MaxRSS 568092 KB; all stages sweep-converged; quartet and matrix at chi128/256 |
| quartet-centered-20261125a | 17764302 | centered.job.json | 1 CPU, 2 GB, 10 min | COMPLETED 0:0 in 1:50, MaxRSS 204512 KB; both saved chi128 states measured, probes aligned at i+1/2 |
| quartet-final-measure-20261125a | 17764463 | final_measure.job.json | 1 CPU, 4 GB, 20 min | COMPLETED 0:0 in 6:32, MaxRSS 238904 KB; full quartet and aligned matrix from both chi256 states; same-state control repeat agrees exactly |

Remote states: `/home/pace/cd-conjecture/runs/quartet-pilot-20261125a/states/`.
Fresh initial state, seed 1. Submission uses working-tree source at initial
upstream 389eadff; exact source hashes are in the uploaded snapshot and
`results/provenance.json`. The new files were uncommitted at submission.
No large MPS files are fetched or committed.

The pilot retained chi64, chi128 and chi256 states (5.6, 17.0 and
61.4 MB respectively, decimal units). The last state has E=-133.56448160561783,
last sweep |Delta E|=1.473e-6 and |Delta S|=.006754, versus requested
1e-9 and 1e-5. It is explicitly `shelved=true, converged=false`.
The ladder is incomplete. Its measured TotalCPU is 3:21:35 for 50:58
elapsed, so the pilot did use approximately four cores on average;
this does not repair the separate timing experiment described below.
The process's own maximum RSS is 542.10 MiB; Slurm's sampled MaxRSS
is 509872 KB. Full per-sweep timings and checkpoint roundtrips are retained.

The thread experiment starts each sweep from the same immutable pilot
chi64 state. It is a timing test on an unconverged low-chi state, not
additional convergence evidence. Per-setting times are 11.39–11.53 s.
The six repeated sweeps consumed 75.3 CPU seconds in 77 elapsed seconds;
the eight allocated CPUs were not used concurrently by MKL. Record actual
thread counts, not just requested values. Complete small artifacts are in
`results/quartet-threads-20261125a/`.

All five allocations are terminal with exit 0:0; the user queue was empty
after the final fetch. Final target/control state SHA-256 hashes, parity,
norm tests, resource records and stage changes are consolidated in
`analysis/final_summary.json`. The target's scientific nonconvergence is
retained explicitly despite its successful allocation exit.
