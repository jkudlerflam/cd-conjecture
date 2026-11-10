This work was done by a Claude Opus 5.5 agent (C195, 2026-11-10), not by Codex. Astra was not used (out of credits).

# C195: p153c at L = 512 on the GPU

Spec: specs/C195_C196_third_size_and_infinite_chain.md (section C195 and the common rules). Log: PROGRESS.md. Derivations and memory estimates: PROOF.md.

## Registration (copied from the spec on 2026-11-10 07:46 EST, before any code was changed or any job submitted)

Key p153c (the collaborator's field h_z = −19.97589125), both parity sectors, the GPU code `codex/work_C190/gpu/z2dmrg_xp.py` with the conventions of `pipeline190.SCHED`, started from scratch with the registered initialization. At most 2 GPUs.

Known obstacle: at χ = 768 and L = 256 the device holds 27 GB of environments on a 40 GB card. L = 512 at χ ≥ 1024 does not fit. Keep only the active environments on the device and the rest on the host (up to 160 GB of host memory per run is allowed for this task) or on scratch. Gate before production: the modified code reproduces an L = 128 run of p153c sweep by sweep to 1e-10 in energy and 1e-8 in entropy against the unmodified GPU code.

Schedule: χ stages rising to a final χ of at least 1024, with a 1536 stage added if the benchmark says it completes inside the 7-day limit. The job must resume from its checkpoint, and a chained resubmission must be in place.

- (T1) Convergence at the final χ: residual ≤ 1e-8, last two sweeps |ΔE| ≤ 2e-8 (the L = 256 bound scaled with L), |ΔS̄_W| ≤ 1e-3, |Δn_mid| ≤ 5e-4, discarded weight ≤ 1e-7. The step of S̄₆ between the last two χ stages is reported against 0.01 (C4), and a failure of C4 alone is reported as in §9.178 and does not by itself void T2 or T3.
- (T2) Pass. On the ground branch c_two(256 → 512) ≥ 1.6 at W = 6 and W = 12, with S̄(256) from the highest χ available at L = 256 for this key, and T1 holds. Also reported, not deciding: c_two at matched χ, c_two from S̄ extrapolated in discarded weight at both sizes, and c_prof(512).
- (T3) Fail. c_two(256 → 512) ≤ 1.3 at W = 6 and W = 12 with T1 holding.
- (T4) Density. n_mid at L = 512 within 0.006 of the L = 256 value of the same key.
- Anything between T2 and T3 is undecided and reported with the numbers.

Definitions used for T1 to T4 are those of `../work_C190/RESULT.md` (S̄_W, c_two, c_prof, n_mid, ground branch), unchanged.

## Code (2026-11-10)

- `z2dmrg_xp.py`. The C190g solver with one change, OFFLOAD (PROOF.md section 1). A GPU sweep keeps on the device only the arrays of the bond being updated (its two environments and two MPS tensors, the MPO, the local solve), and every other environment and tensor on the host. Per bond it uploads one environment and one tensor and downloads one of each. Every arithmetic call receives the same arrays in the same order as in the C190g code. Only dmrg() changed. The in-memory path is unchanged and still runs the host sweeps (χ < 128). No tolerance or algorithmic setting changed, complex128 throughout.
- `run195.py`. run190g.py with the offloaded solver, with intermediate stages that run to a stop rule (key 'conv' of the configuration, same rule as the final stage of run190g), and with the final measurements on the returned host MPS without a further copy. A configuration without 'conv' runs as one dmrg() call, as in run190g. Output JSON, checkpoint format (MPS.save of C190g, unchanged) and resume behaviour are those of run190g. Every history record carries 'stage'.
- `pipeline195.py`. The two production runs, one per GPU, resume and skip as pipeline_x1024.py.
- `bench195.py`, `gate195.py`, `validate_off.py`, `analyze195.py`. Benchmark, gate comparison, local check, analysis.
- `ref/`. Unmodified copies of the C190g files used as the reference in the gate (md5 equal to the files in /scratch/lustre/jkudlerflam/c18x/c190g/).

## Gates

Local numpy check (validate_off.py, L = 20, p153c, logs/validate_off_L20.log). Against the unmodified code, with in-memory sweeps followed by offloaded ones, with a resume that starts offloaded, and in the odd sector from a product state. The mixer sweeps are bit identical, the mixer-off sweeps differ by at most 3.7e-14 in E and 4.1e-12 in the entropy at any cut, also with one BLAS thread. Matvec counts identical. Passed.

GPU gate (job 17611626, one A100-PCIE-40GB on apollo02, typhon_logs/gate_17611626.txt, run logs typhon_logs/ref.log, off.log, res1.log, res2.log). p153c, L = 128, even sector, rand:1, a short version of pipeline190.SCHED[128] (χ = 32, 32, 64, 64 with mixer on the host, 128 and 192 with mixer and 256 without on the GPU, then two sweeps at χ = 384 with residual 1e-9; cfg/gate_p153c_L128_P0.json), production settings (GPU from χ = 128, device SVD). Three runs on the same GPU: the unmodified C190g code (ref/run190g.py), the C195 code, and the C195 code stopped after 8 sweeps and resumed from its checkpoint.

| comparison with the unmodified code | worst |ΔE| over 9 sweeps | worst |ΔS| at any cut | worst Δ⟨Z⟩ | final exact entropies | final density |
|---|---|---|---|---|---|
| C195, uninterrupted | 0 | 0 | 0 | 0 | 0 |
| C195, stopped after 8 sweeps and resumed | 0 | 6.9e-11 | 2.5e-11 | 5.4e-11 | 1.2e-11 |

The offloaded run is bit identical to the unmodified one in every sweep, including the energy variance (6.077053350850292e-05 in both). The resumed run differs only after the resume, at the level of rounding (the resume re-canonicalizes the checkpoint state). The gate (1e-10 in E, 1e-8 in S) passes. Final values E = −1290.1134150841, S̄₆ = 1.854346. Device memory with offloading at L = 128, χ = 384 was 15 to 19 GB in use on the card, shared by the three processes. Sweep times with three processes on one card were 233 to 300 s (offloaded) against 258 to 265 s (unmodified) at χ = 384.

## Checkpoints of the final states (coordinator's addition, 2026-11-10)

Each run writes its checkpoint after every sweep with MPS.save of the C190g code (unchanged, numpy arrays M0 … M511, n0, P, meta, uncompressed npz), so the offloading does not change the format, and the final-state checkpoint is the one written after the last sweep, which is the state that is measured. Nothing deletes it when the run ends, and a resubmission of a finished run only prints "already done". The run json has the keys meas193.py reads (status, cfg with model, hz, gauge, L, P, meas.gwin, E_off), plus an extra 'conv' entry in cfg. The final states will be at

    /scratch/lustre/jkudlerflam/c19x/c195/out512/p153c_L512_P0_rand1.json  /scratch/lustre/jkudlerflam/c19x/c195/out512/ckpt/p153c_L512_P0_rand1.npz
    /scratch/lustre/jkudlerflam/c19x/c195/out512/p153c_L512_P1_rand1.json  /scratch/lustre/jkudlerflam/c19x/c195/out512/ckpt/p153c_L512_P1_rand1.npz

(the two lines are in the "<run.json> <ckpt.npz>" form of meas193.py batch). At χ = 1536 each npz is about 37 GB. Before a run is done, the npz holds the last finished sweep, and the run json has status running.

## Benchmark at L = 512 (job 17611626, typhon_logs/bench_L512_17611626.json and .log)

bench195.py: p153c at L = 512, a random Z₂-block MPS with every bond at min(χ, 2^dist), canonicalized and with all environments built on the GPU and stored on the host as in an offloaded sweep, then three consecutive mid-chain bond updates (bonds 255 to 257) exactly as in the offloaded sweep (uploads, TwoSite, kron4, Davidson with 100 matvecs, device SVD split, env_left, downloads) plus the entropy SVD of a left step. One A100-PCIE-40GB, 4 BLAS threads. Seconds per bond update, mean of the three.

| χ | per bond | matvec (100) | preconditioner | Davidson bookkeeping | split | kron4 eigh | uploads + downloads | entropy SVD | device memory in use (peak) | host memory (VmRSS) | of which environments, MPS |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 768 | 1.42 | 0.38 | 0.27 | 0.45 | 0.20 | 0.05 | 0.022 | 0.05 | 11.7 GB | 32.6 GB | 23.3, 9.3 GB |
| 1024 | 2.33 | 0.76 | 0.31 | 0.75 | 0.33 | 0.07 | 0.034 | 0.06 | 20.3 GB | 57.5 GB | 41.3, 16.6 GB |
| 1536 | 5.83 | 2.27 | 0.86 | 1.73 | 0.66 | 0.11 | 0.075 | 0.11 | 14.1 GB | 128.0 GB | 92.7, 37.1 GB |

- The host memory agrees with the estimate of PROOF.md section 2 (33, 58, 130 GB). At χ = 1536 a run needs 128 GB on the host, inside the 160 GB allowed. χ = 2048 would need about 230 GB and is not possible within the allowance.
- The device memory is well inside the 40 GB card at every χ. (The 1024 row reads higher than the 1536 row because the cupy pool kept blocks of the earlier χ, the pool was freed between rows but not defragmented. Both are peaks of the pool, not of live arrays.)
- The transfers cost 1.3 % to 1.6 % of a bond update.
- Check against production. At χ = 768 the per-bond cost is 0.32 s plus 11 ms per matvec. The L = 256 production runs of C190g took 594 to 893 s per χ = 768 sweep (510 bond updates, 84 to 150 matvecs on average), that is 1.16 to 1.75 s per bond, and the formula gives 1.24 to 1.97 s. So the benchmark is a fair guide.

## Schedule (fixed 2026-11-10 08:25 EST, before production)

Final χ = 1536. Both sectors, rand:1, field −19.97589125, model p153 without gauge, Davidson settings and tolerances of C190 (maxit 300, mmax 32, residual targets as below), GPU from χ = 128 with the device SVD, gwin = 48, no variance (as pipeline190 for L > 128). pipeline195.SCHED512:

- Fixed stages as pipeline190.SCHED[256], two sweeps each: χ = 32 (mixer 1e-3, residual 1e-7), 64 (1e-4, 1e-8), 128 (1e-4, 1e-8), 192 (1e-5, 1e-8), 256 (1e-5, 1e-8), 384 (no mixer, 1e-8), 512 (no mixer, 1e-9).
- χ = 768 and χ = 1024, each run to the stop rule (residual 1e-9, at least 3 and at most 6 sweeps, stop when |ΔE| < 4e-9 and |ΔS̄₆| < 5e-4 between consecutive sweeps).
- Final χ = 1536, the same rule with at least 3 and at most 8 sweeps.

Why stages run to a stop rule. R8 of C190g showed that two sweeps per stage leave the stage unconverged in sweeps, so that the C4 step from the previous stage mixed sweep convergence with χ dependence (at L = 256 the first χ = 768 sweep alone moved S̄₆ by 0.06). With χ = 768 and 1024 converged in sweeps, the C4 step 1024 → 1536 measures the χ dependence alone, the discarded-weight extrapolation of T2 has two converged stages at L = 512, and the matched-χ comparison with L = 256 is available at χ = 768 (outG) and χ = 1024 (outX1024, if it is done). The stop tolerance is the L = 256 one (2e-9) doubled with L, as T1 doubles C2.

Time estimate (from the benchmark and the C190g sweep times, for one run per A100). A sweep is 1022 bond updates. With 100 to 150 matvecs per solve on average (84 to 150 in the L = 256 production at χ = 768):
- χ = 768: 1.4 to 2.0 s per bond, 24 to 34 min per sweep, 3 to 6 sweeps, 1.2 to 3.4 h.
- χ = 1024: 2.3 to 3.2 s per bond, 40 to 55 min per sweep, 3 to 6 sweeps, 2.0 to 5.5 h.
- χ = 1536: 5.8 to 8.2 s per bond, 1.6 to 2.3 h per sweep, 3 to 8 sweeps, 4.9 to 18.5 h.
- Fixed stages χ ≤ 512: twice the L = 256 sweep times, about 2.4 h.
- Checkpoints: 37 GB per sweep at χ = 1536, about a minute each at 1 GB/s.
- Final host measurements (entropies, one-body matrices, C_P, 4 threads): 775 to 808 s at L = 256, χ = 768, scaled by L χ³ (a factor 16), about 3.5 h.

Total about 14 to 33 h per run. The two runs go in parallel on two GPUs. Even with a factor 2 for longer Davidson solves at L = 512, this is inside the 7-day limit, so the 1536 stage is included. The estimate is uncertain mainly through the matvec counts at L = 512 and the number of sweeps the stop rules take.

Jobs. The production job uses the gpu-medium QOS (3 days, the slot this account can get first), and a chained resubmission of the same script in the gpu-long QOS (7 days) with --dependency=afterany resumes whatever is unfinished.

## Production (status 2026-11-10 08:50 EST)

Jobs. 17611738 (c195_prod, partition apollo, QOS gpu-medium, 3-day limit, node apollo02, 2 A100-PCIE-40GB, 8 cores, 320 GB), started 08:25:44, runs p153c_L512_P0_rand1 on GPU 0 and p153c_L512_P1_rand1 on GPU 1 (pipeline195.py, 4 threads each). 17611739 is the same script in QOS gpu-long (7-day limit) with --dependency=afterany:17611738. It resumes every unfinished run from its checkpoint (at most one sweep is lost) and skips finished ones. Script typhon/c195_prod.slurm, settings C190G_GPU_MIN_CHI=128, C190G_SVD_DEV=1, C195_OFFLOAD=1. The configurations are in out512/cfg/ on Typhon.

Output layout on Typhon, /scratch/lustre/jkudlerflam/c19x/c195/: out512/<tag>.json (run190g format plus 'stage', 'host_rss_GB', 'dev_used_GB' per sweep and 'conv' in cfg), out512/logs/<tag>.log (one JSON line per sweep, '# ckpt' lines with the checkpoint time), out512/ckpt/<tag>.npz (up to about 37 GB each, do not copy), logs/pipeline_<jobid>.log, logs/c195_prod_<jobid>.out, logs/prod_smi_<jobid>.log (GPU memory and load and host memory every 5 minutes). Gate and benchmark files are in gate/ and logs/.

First sweeps (logs/first_sweeps_0850.txt, made with `python3 -E -s check_first.py typhon_out/out512/logs`). Every local solve met its residual target in all five sweeps of both runs (1e-7 at χ = 32, 1e-8 from χ = 64 on). The comparison is with the C190g L = 256 and L = 128 runs of the same key, initialization and schedule at the same sweep index. The bulk energy per site is (E(512) − E(256))/256 against (E(256) − E(128))/128.

| sector | sweep (χ) | E(512) | bulk (512 − 256) | bulk (256 − 128) | difference | n_mid(512) | n_mid(256) | sweep time |
|---|---|---|---|---|---|---|---|---|
| even | 2 (64) | −5130.54321948 | −10.001119767 | −10.001119876 | 1.1e-7 | 0.2140 | 0.2140 | 241 s |
| even | 3 (64) | −5130.54312963 | −10.001119455 | −10.001119978 | 5.2e-7 | 0.2143 | 0.2143 | 317 s |
| even | 4 (128, first GPU sweep) | −5130.54409020 | −10.001121496 | −10.001122036 | 5.4e-7 | 0.2152 | 0.2153 | 640 s |
| odd | 2 (64) | −5130.54315837 | −10.001119411 | −10.001120110 | 7.0e-7 | 0.2149 | 0.2140 | 244 s |
| odd | 3 (64) | −5130.54304017 | −10.001119114 | −10.001119986 | 8.7e-7 | 0.2154 | 0.2142 | 318 s |
| odd | 4 (128, first GPU sweep) | −5130.54396004 | −10.001120986 | −10.001122054 | 1.1e-6 | 0.2163 | 0.2153 | 652 s |

- The bulk energy per site agrees with the L = 256 run at the same stage to 1e-7 to 1e-6, the size of the unconverged-stage differences, and n_mid agrees with L = 256 to 0.001 (T4 asks for 0.006 at the end). At χ = 32 the states are still far from converged at both sizes (n_mid 0.18 against 0.19) as at L = 256.
- Sweep times at χ = 128 (640 to 652 s) are twice the L = 256 times (300 to 327 s), as expected for twice the bonds, so the offloading costs nothing visible. Device memory at χ = 128 is 1 GB per run and host memory 1.6 GB.
- One warning at the start of both runs, which needs no fix. At L = 512 the right canonicalization of the random start (rand:1, z2dmrg.random_mps with D0 = 8, the registered initialization) overflows. The leftward-pushed R factor ends with a norm above 1e154, its squared norm overflows, and the final normalization of M₀ divides by inf ("overflow encountered in dot", "invalid value encountered in divide"). This happens in the unchanged C190g and CPU code alike (L = 256 stays below the overflow), and it touches M₀ only. M₁ to M₅₁₁ are orthonormal QR factors and finite. The first local problem of the first sweep (bond 0, dimension below 128) is solved by exact diagonalization, which does not read M₀, so M₀ is replaced before it is ever used. The finite energies, residuals and densities of all sweeps confirm that no NaN entered.

Estimated finish. From the benchmark and the measured sweep times (the χ = 128 sweeps match twice the L = 256 times), each run should reach the end of χ = 512 at about 11:00 EST, the end of χ = 768 between about 12:00 and 14:30, of χ = 1024 between about 14:00 and 20:00, finish the χ = 1536 stage between about 19:00 on 2026-11-10 and 14:30 on 2026-11-11, and finish the host measurements about 3.5 h later. So both runs should be done between about 22:30 EST on 2026-11-10 and 18:00 EST on 2026-11-11. This is an extrapolation from benchmark times with 100 to 150 matvecs per local solve and from the stop rules taking 3 to 6 (8 at the final χ) sweeps. Longer solves at L = 512 could stretch it by up to a factor 2 (to about 2026-11-12 evening), still inside the 3-day limit of the first job, and the chained job covers anything beyond.

## Collection and analysis (tested 2026-11-10 08:50 on the running output, 132 kB)

From inside codex/work_C195, small JSON and log files only (no checkpoints). The remote globs are quoted for zsh.

    mkdir -p typhon_out/out512/logs typhon_out/logs
    scp 'typhon:/scratch/lustre/jkudlerflam/c19x/c195/out512/*.json' typhon_out/out512/
    scp 'typhon:/scratch/lustre/jkudlerflam/c19x/c195/out512/logs/*.log' typhon_out/out512/logs/
    scp 'typhon:/scratch/lustre/jkudlerflam/c19x/c195/logs/pipeline_*.log' 'typhon:/scratch/lustre/jkudlerflam/c19x/c195/logs/c195_prod_*.out' 'typhon:/scratch/lustre/jkudlerflam/c19x/c195/logs/prod_smi_*.log' typhon_out/logs/
    python3 -E -s analyze195.py > typhon_out/analyze195.txt

Add --partial to read running runs at their last mixer-off sweep (the χ = 384 stage is the first one with entropies, about 2.5 h after the start). Without it only runs with status done are used. The L = 256 side is read from ../work_C190/gpu/typhon_out/outG and ../work_C190/gpu/typhon_out/outX1024, so the χ = 1024 continuation of C190g R8 enters once its outputs are collected there (`scp 'typhon:/scratch/lustre/jkudlerflam/c18x/c190g/outX1024/*.json' ../work_C190/gpu/typhon_out/outX1024/`, read only on Typhon). `python3 -E -s check_first.py typhon_out/out512/logs` prints the sweep-by-sweep comparison with L = 256.

What analyze195.py prints: per run the stage history (last sweep of each χ), T1 with every number and its bound, the C4 step, the ground branch at both sizes, c_two(256 → 512) for the ground branch, each sector and all four parity pairings, the T2/T3 verdict as registered (T1 required for the ground-branch run at L = 512, S̄(256) per sector from the highest-χ finished run), c_prof(512), c_two at matched χ, c_two from the R8 extrapolation in the discarded weight (L = 512 stages 1024 and 1536, L = 256 stages 768 and 1024, the latter only when outX1024 has the run), and T4. Tests: on the C190 outputs as L = 128 → 256 it reproduces R6 and R7 (ground-branch c_two 1.902, 1.936, 1.971, c_prof 2.028; logs/analyze_test_256_128.txt), on the gate output of run195.py as L = 64 → 128 it reads the new 'stage' records (logs/analyze_test_gate128_64.txt), and on the running L = 512 output it reports that no mixer-off sweep exists yet.

## What is left

- Production runs to completion (both sectors). The coordinator collects and runs analyze195.py, then measures the final checkpoints with meas193.py (paths above).
- The L = 256 side at χ = 1024 (C190g R8, job 17609630) was still pending at 08:50. Without it S̄(256) comes from χ = 768 and the extrapolated c_two is not available.
