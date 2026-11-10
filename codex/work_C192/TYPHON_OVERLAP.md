# C190 and C189 payloads running inside other allocations (coordinator, 2026-11-09)

Why. The Typhon queue is memory-bound, and the queued C190 and C189 jobs have a scheduler start estimate that moves between 2026-11-10 and 2026-11-14. Several of our running cd jobs hold 48 cores and 180 to 320 GB but use a few cores and a few GB. So the queued payloads were started at 2026-11-09 20:23 EST as extra steps inside those allocations (`srun --jobid=<donor> --overlap`). Nothing was cancelled, paused, held or released, and no other user's job is affected.

What runs where.

| payload | queued job | donor job (node) | step | donor ends at the latest |
|---|---|---|---|---|
| C190 node A (decisive, φ = 1.53 bond at L = 64, 128, 256), 8 slots × 6 threads | 17606939 | c185 17604455 (node18) | .2 | about 41 h, earlier if its last blockk.py task ends |
| C190 node B (controls, symmetric-line and rational bonds), 6 slots × 6 threads | 17606940 | c166 17602270 (node26) | .1 | about 22 h |
| C189 inf_s316 (10 threads) | 17606915 | c186 17604527 (node11) | .0 | about 18 h |
| C189 inf_ff2v (10 threads) | 17606918 | c186 17604527 (node11) | .1 | about 18 h |
| C189 inf_s18 (8 threads) | 17606916 | c173 17603887 (node19), relaunched 2026-11-10 00:55 EST | .1 | about 23.5 h |
| C189 inf_s16 (6 threads) | 17606917 | c175 17604038 (node12), relaunched 2026-11-10 00:55 EST | .3 | about 8.5 h |
| C189 inf_c1 (4 threads) | 17606919 | c175 17604038 (node12), relaunched 2026-11-10 00:55 EST | .2 | about 8.5 h |

How it hands over. Each step runs the same body as its queued job (launchers `codex/work_C190/typhon/overlap_run.sh` and `codex/work_C189/typhon/overlap_inf.sh`, both also in scratch), and polls the queue. When its queued job is RUNNING the step cancels itself, and the queued job resumes from the checkpoints. Both codes write checkpoints atomically (run190.py already did, run_inf.py was patched for this on 2026-11-09, two lines). When a donor ends, its steps die with it, and the payload then waits for its queued job unless it is relaunched.

Check (one ssh call).

    squeue -s -u jkudlerflam -h -o "%i %M %N" | grep -v "batch\|extern"
    squeue -h -u jkudlerflam -j 17606939,17606940,17606915,17606916,17606917,17606918,17606919 -o "%i %j %T %S"
    tail -n 3 /scratch/lustre/jkudlerflam/c18x/c190/logs/pipeline_A_ov.log /scratch/lustre/jkudlerflam/c18x/c190/logs/overlap_*.log

Relaunch rule. Relaunch a payload only if all three hold: its step is gone from `squeue -s`, its queued job is still PENDING, and its outputs are not complete (C190: `logs/overlap_<A|B>.log` has no "payload finished" line. C189: the tag's last χ is missing from `out/<tag>.jsonl`). Pick as donor one of our own running c1xx jobs with idle cores (node CPULoad well below 48) and the most time left. Never use another user's job, never cancel a running job, never pause one except under the pause rule below, never release a held job. Scale the slots so that slots × 6 threads fits in the idle cores. From the login node (the srun client does no computing):

    cd /scratch/lustre/jkudlerflam/c18x/c190
    nohup setsid srun --jobid=<donor> --overlap -n1 -c48 bash overlap_run.sh <A|B> <slots> <queued job id> > logs/srun_<A|B>.out 2>&1 < /dev/null &

    cd /scratch/lustre/jkudlerflam/c18x/c189
    nohup setsid srun --jobid=<donor> --overlap -n1 -c48 bash typhon/overlap_inf.sh <tag> <threads> <queued job id> > logs/srun_<tag>.out 2>&1 < /dev/null &

Put the `cd` and the `nohup` on separate lines (a `cd ... && nohup ... &` list keeps the ssh session open). Never start a second copy of a payload that is still running, and never start one while its queued job is RUNNING.

Collection is unchanged (see each RESULT.md). Analyses must use only runs with status done.

Pause rule (approved by Jonah, 2026-11-09, about 20:48 EST). It covers one case only, C190 node A, and one job only, c182 (17604443, the phase-diagram screen, node20).

- When. All of these hold: node A's step is gone, its queued job 17606939 is still PENDING, its outputs are not complete, and no running c1xx job has enough idle cores for it (about 24).
- How. Freeze c182 first, then start node A inside it. A step that exists when the signal is sent would be frozen too.

      scancel --signal=STOP --full 17604443
      cd /scratch/lustre/jkudlerflam/c18x/c190
      nohup setsid bash -c 'srun --jobid=17604443 --overlap -n1 -c48 bash overlap_run.sh A 8 17606939; scancel --signal=CONT --full 17604443' > logs/srun_A.out 2>&1 < /dev/null &

  The second half of that command resumes c182 as soon as node A's step ends, whether the payload finished or the queued job took over.
- Check. After the STOP, confirm that c182's processes are stopped (state T) and that node A's log advances. This use of STOP has not been tried on this cluster. If it does not freeze the job, send CONT and fall back to waiting for idle cores.
- Backstop. At every check-in, if c182 was paused and node A's step is no longer in `squeue -s`, send `scancel --signal=CONT --full 17604443`. CONT on a job that is not stopped is harmless.
- Limits. c182 keeps its allocation and memory and loses nothing, but its wall clock keeps running while it is frozen. Never cancel it. Pause nothing else without asking Jonah. Record the pause and resume times here.

Status 2026-11-09 20:54 EST. All seven steps are running and all seven queued jobs are PENDING. Nothing has been paused.

Status 2026-11-09 21:38 EST. All seven steps are still running in the same donors, and all seven queued jobs are PENDING with start estimates on 2026-11-10 after 17:00. Nothing has been paused.

## GPU production of C190 (added 2026-11-09 22:20 EST)

The GPU backend (REPORT §9.168) runs the registered C190 configurations as job 17607561 on apollo02, with the resume job 17607562 waiting on it. Neither is an overlap step and neither is to be touched. Output is in /scratch/lustre/jkudlerflam/c18x/c190g/outG/. Check and collect with the commands in `codex/work_C190/gpu/RESULT.md` (never copy outG/ckpt). A finished run has "status": "done" in its json. When the p153c runs at L = 64, 128 and 256 are done, evaluate C1 to C6 and V1 to V6 of `codex/work_C190/RESULT.md` with the amendment in `gpu/RESULT.md` (lowest energy among CPU and GPU runs of a tag).

## C187d, the χ = 1024 continuation of the L = 256 stage (added 2026-11-09 22:55 EST)

Jonah approved it on 2026-11-09. Registration, launch record, check and collection commands are in `codex/work_C187q/L256_job/chi1024/README.md`. Step A is 17603971.1 (six candidate runs, inside c177), step B is 17604038.1 (ctl2 and the L = 128 partners, inside c175), and the queued fallback is 17608223 (c187d). A run is complete at 28 jsonl lines (L = 256) or 27 (L = 128). When all ten are complete, collect, apply the registered readings R1 and R2 and the verdict rule E1 to E3, write the REPORT section, and cancel the fallback job only if it is still pending and every tag is complete (it is ours and would do nothing).

Log, 2026-11-10 00:55 EST. The donor c178 (17603981) completed at 22:31 EST on 2026-11-09 and took the three C189 steps it hosted with it (inf_s18 had reached χ = 384, inf_s16 χ = 512, inf_c1 χ = 384). All three conditions of the relaunch rule held, so they were relaunched from their checkpoints, inf_s18 inside c173 (node19) and inf_s16 and inf_c1 inside c175 (node12). The check above lists steps but not which payloads should exist, so compare its step list against the table at every check-in.
