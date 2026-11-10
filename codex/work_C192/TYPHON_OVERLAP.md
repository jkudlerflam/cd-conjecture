# C190 and C189 payloads running inside other allocations (coordinator, 2026-11-09)

Why. The Typhon queue is memory-bound, and the queued C190 and C189 jobs have a scheduler start estimate that moves between 2026-11-10 and 2026-11-14. Several of our running cd jobs hold 48 cores and 180 to 320 GB but use a few cores and a few GB. So the queued payloads were started at 2026-11-09 20:23 EST as extra steps inside those allocations (`srun --jobid=<donor> --overlap`). Nothing was cancelled, paused, held or released, and no other user's job is affected.

What runs where.

| payload | queued job | donor job (node) | step | donor ends at the latest |
|---|---|---|---|---|
| C190 node A (decisive, φ = 1.53 bond at L = 64, 128, 256), 8 slots × 6 threads | 17606939 | c185 17604455 (node18) | .2 | about 41 h, earlier if its last blockk.py task ends |
| C190 node B (controls, symmetric-line and rational bonds), 6 slots × 6 threads | 17606940 | c166 17602270 (node26) | .1 | about 22 h |
| C189 inf_s316 (10 threads) | 17606915 | c186 17604527 (node11) | .0 | about 18 h |
| C189 inf_ff2v (10 threads) | 17606918 | c186 17604527 (node11) | .1 | about 18 h |
| C189 inf_s18 (8 threads) | 17606916 | c178 17603981 (node39) | .3 | about 9.5 h |
| C189 inf_s16 (6 threads) | 17606917 | c178 17603981 (node39) | .2 | about 9.5 h |
| C189 inf_c1 (4 threads) | 17606919 | c178 17603981 (node39) | .4 | about 9.5 h |

How it hands over. Each step runs the same body as its queued job (launchers `codex/work_C190/typhon/overlap_run.sh` and `codex/work_C189/typhon/overlap_inf.sh`, both also in scratch), and polls the queue. When its queued job is RUNNING the step cancels itself, and the queued job resumes from the checkpoints. Both codes write checkpoints atomically (run190.py already did, run_inf.py was patched for this on 2026-11-09, two lines). When a donor ends, its steps die with it, and the payload then waits for its queued job unless it is relaunched.

Check (one ssh call).

    squeue -s -u jkudlerflam -h -o "%i %M %N" | grep -v "batch\|extern"
    squeue -h -u jkudlerflam -j 17606939,17606940,17606915,17606916,17606917,17606918,17606919 -o "%i %j %T %S"
    tail -n 3 /scratch/lustre/jkudlerflam/c18x/c190/logs/pipeline_A_ov.log /scratch/lustre/jkudlerflam/c18x/c190/logs/overlap_*.log

Relaunch rule. Relaunch a payload only if all three hold: its step is gone from `squeue -s`, its queued job is still PENDING, and its outputs are not complete (C190: `logs/overlap_<A|B>.log` has no "payload finished" line. C189: the tag's last χ is missing from `out/<tag>.jsonl`). Pick as donor one of our own running c1xx jobs with idle cores (node CPULoad well below 48) and the most time left. Never use another user's job, never cancel or pause a running job, never release a held job. Scale the slots so that slots × 6 threads fits in the idle cores. From the login node (the srun client does no computing):

    cd /scratch/lustre/jkudlerflam/c18x/c190
    nohup setsid srun --jobid=<donor> --overlap -n1 -c48 bash overlap_run.sh <A|B> <slots> <queued job id> > logs/srun_<A|B>.out 2>&1 < /dev/null &

    cd /scratch/lustre/jkudlerflam/c18x/c189
    nohup setsid srun --jobid=<donor> --overlap -n1 -c48 bash typhon/overlap_inf.sh <tag> <threads> <queued job id> > logs/srun_<tag>.out 2>&1 < /dev/null &

Put the `cd` and the `nohup` on separate lines (a `cd ... && nohup ... &` list keeps the ssh session open). Never start a second copy of a payload that is still running, and never start one while its queued job is RUNNING.

Collection is unchanged (see each RESULT.md). Analyses must use only runs with status done.
