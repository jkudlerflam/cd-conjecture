This log was written by a Claude Opus 5.5 agent (C195, 2026-11-10), not by Codex.

# C195 progress log

Times are EST, from Typhon's clock where a job is involved.

1. 07:43. Read the spec, `../work_C190/gpu/RESULT.md`, `../work_C190/RESULT.md` and the C190g code. Typhon state: gpu-medium held by ddL_g (to 10:29) and c193, gpu-short by gpu_lean1 and gpu_lean2, gpu-long by ddH_g (to 13:30). gpu-veryshort is not in this account's QOS list. c193 ended a few minutes later, which freed one gpu-medium slot.
2. 07:46. Folder created, code copied (z2dmrg.py, models.py, pipeline190.py, gpu/z2dmrg_xp.py; ref/ holds unmodified copies of gpu/z2dmrg_xp.py, gpu/run190g.py, z2dmrg.py and models.py, md5 identical to the C190g files on Typhon). Registration T1 to T4 copied into RESULT.md.
3. z2dmrg_xp.py: OFFLOAD added (PROOF.md section 1). Only dmrg() changed. The in-memory path is kept as it was and is used for host sweeps (χ < GPU_MIN_CHI) and when OFFLOAD is off.
4. Local numpy check, validate_off.py at L = 20 (logs/validate_off_L20.log, 4 threads, and logs/validate_off_L20_1thread.log). Three tests against the unmodified code: a run whose first sweeps are in memory and later ones offloaded (GPU_MIN_CHI = 16), a restart whose first sweep is offloaded (the resume path), and the odd sector from a product state. The mixer sweeps are bit identical. The mixer-off sweeps differ by at most 3.7e-14 in E and 4.1e-12 in the entropy at any cut (the same with one BLAS thread), against the gate's 1e-10 and 1e-8. Matvec counts identical everywhere. Passed.
5. run195.py (driver), with stop-ruled intermediate stages ('conv') and the resume logic per stage. Local test at L = 16 with OFFLOAD_TEST: a run interrupted after 5 sweeps and resumed reproduces the uninterrupted one to 5e-14 in every sweep energy, the stages end as their rules say.
6. 07:51. Gate and benchmark job 17611626 (c195_gate, gpu-medium, 1 A100 on apollo02, 160 GB) submitted and started. Smoke test of bench195.py at L = 64 passed.
7. analyze195.py written and tested on the C190 outputs as L = 128 → 256 (logs/analyze_test_256_128.txt). It reproduces the R6 and R7 readings (ground-branch c_two(128 → 256) = 1.902, 1.936, 1.971, c_prof(256) = 2.028, C4 steps 0.0871 and 0.0269).
8. 08:20. Gate passed (typhon_logs/gate_17611626.txt). The offloaded run is bit identical to the unmodified C190g code in all 9 sweeps and in the final measurements, the resumed run differs by at most 6.9e-11 in S after the resume.
9. 08:22. Benchmark at L = 512 done (typhon_logs/bench_L512_17611626.json): 1.42, 2.33, 5.83 s per bond at χ = 768, 1024, 1536 with 100 matvecs, host memory 32.6, 57.5, 128.0 GB, device at most 20 GB. Schedule fixed with final χ = 1536 and written into RESULT.md at 08:25, with the time estimate.
10. 08:25. Production submitted: 17611738 (gpu-medium, 3 days) started at 08:25:44 on apollo02 GPUs 0 and 1, and 17611739 (gpu-long, 7 days, afterany) as the chained resume.
11. Coordinator's addition on the final checkpoints answered in RESULT.md: format unchanged (MPS.save of C190g), kept after the run, paths listed.
12. 08:50. First five sweeps of both runs checked (logs/first_sweeps_0850.txt): residual targets met, bulk energy per site within 1e-7 to 1.1e-6 of the L = 256 run at the same stage, n_mid within 0.0012. The overflow warning at the start (right canonicalization of the random start at L = 512) only touches M₀, which the exact solve of bond 0 replaces, so nothing was changed. Collection commands and analyze195.py tested on the running output. Stopped here as the spec says.
