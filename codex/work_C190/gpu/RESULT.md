This work was done by a Claude Opus 5.5 agent (C190g, 2026-11-09), not by Codex. Astra was not used (out of credits).

# C190g: a GPU backend for the C190 DMRG

Spec: specs/C190g_gpu_port.md. Log: PROGRESS.md in this folder. Registration of C190: ../RESULT.md (unchanged, amended only by the paragraph below).

## Code

- z2dmrg_xp.py. Array-module port of z2dmrg.py. Every sweep routine runs in the module of its inputs (numpy or cupy), with the same operation order as z2dmrg.py. On the GPU the whole local solve stays on the device (two-site tensor, both environments, MPO tensors, the eigenbases of the (s1, s2)-block Kronecker preconditioner, the Davidson basis). Inside the Davidson loop only the m × m Ritz matrix (m ≤ 32), its eigenvectors and the residual norm cross to the host. Hamiltonian classes and final measurements (entropies, measure, variance, parity) are the unchanged host code of z2dmrg.py. complex128 throughout, float64 for the real-gauge bonds, no tolerance changed.
- run190g.py. run190.py with the backend switch. Same output JSON (same keys, checked on the G3 output) and the same checkpoint format (numpy arrays), so a run can resume on either backend and analyze190.py reads the output folder unchanged.
- pipeline190g.py. One run per GPU, configurations built by pipeline190.cand_task (same schedules, seeds, initializations).
- Production settings. Sweeps with χ < 128 run on the host (C190G_GPU_MIN_CHI=128), the truncation SVD, mixer eigh and entropy SVD run on the device (C190G_SVD_DEV=1). Both choices follow from G4 below.
- Typhon environment. cupy 13.6.0 (cupy-cuda12x, ~/.local, python 3.12) works on apollo02 only after CUDA_PATH points at cudahdr/, a symlink to the CUDA 12.6 runtime headers that already exist in ~/.local/lib/python3.13 (nvidia-cuda-runtime-cu12 wheel). Without it NVRTC fails on cuda_fp16.h. Nothing was installed. On apollo01 the python process hung at the cupy import in two attempts, so all jobs are pinned to apollo02. GPU: NVIDIA A100-PCIE-40GB, driver 580.105.08.

## Gates

G1, numpy backend against z2dmrg.py (local, logs/G1_numpy.log). Passed.
- L = 10, validate.py cases (p153, rat, sym40 with and without the staggered gauge, xx1, xx2), both sectors, random and product starts. Every difference between z2dmrg_xp (numpy) and z2dmrg.py is exactly 0.0 (energy, entropy at every cut, ⟨Z⟩, one-body matrices, per-sweep energies). Errors against ED are those of validate.py (worst candidate bond 3.1e-12 in E, 6.3e-10 in S). The xx2 odd-sector row is the known degeneracy (../RESULT.md).
- L = 32, p153c, 6 sweeps (χ = 32 and 64 with mixer, 96 without): worst difference 1.8e-11 in E and 1.1e-11 in the entropy at any cut, identical matvec counts. Two runs of z2dmrg.py itself on the same input differ by 4.3e-11 after the first mixer sweep (multithreaded BLAS), so this is the reproducibility floor of the CPU code.

G2, cupy backend against ED at L = 10 on an A100 (typhon_logs/G2_*.log). Passed.
- Host SVD: worst candidate-bond errors against ED 6.5e-12 in E and 6.0e-10 in S (p153 6.5e-12, 3.8e-10; rat 4.6e-12, 3.9e-10; sym40/stag 3.9e-12, 6.0e-10; sym40 4.1e-12, 5.6e-10), one-body matrices 3.5e-10, C_P 1.2e-8, variance below 6e-11. Against z2dmrg.py on the same input, 2.3e-14 in E and 3.9e-12 in S.
- Device SVD: against ED 6.4e-12 in E and 5.9e-10 in S, against z2dmrg.py 2.7e-14 in E and 2.6e-12 in S.
- xx1 agrees with ED to 7e-15. xx2 even sector to 1e-15, odd sector the known degeneracy, as on the CPU.
- L = 32 sweep by sweep (GPU at every stage, host SVD). After the first mixer sweep at χ = 32 (residual target 1e-7) the energies differ by 1.4e-7, then 6.9e-10, 4.2e-11, 8.4e-11, and at the two χ = 96 sweeps without mixer 1.7e-12 and 1.0e-11 in E and 2.1e-8 in the entropy at any cut, with identical matvec counts except one sweep (48.06 against 48.08). The early difference is a different but equally converged state at a loose residual and truncation, which the later sweeps remove. This is larger than the CPU run-to-run floor, but it is the case only in the growth stages, and G3 tests what matters (the converged state).
- Production path at L = 32 (host for χ < 64, then GPU with device SVD): bit-identical on the host sweeps, then 1.0e-13 in E and 6.6e-11 in the entropy at any cut.

G3, the L = 64 pilot at the collaborator's field with the pilot's schedule (local/cfg), GPU at every stage, host SVD (typhon_logs/G3_*.log and p153c_L64_P*.json). Passed. Both runs stop after the same 13 sweeps as the pilot.

| sector | E last sweep, GPU | reference | ΔE | S̄₁ GPU (ref) | S̄₆ GPU (ref) | ΔS̄ max |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| even | −650.0414852453 | −650.0414852455 | 2.5e-10 | 1.614034 (1.614034) | 1.608546 (1.608546) | 5.9e-7 |
| odd | −650.0414362272 | −650.0414362271 | 1.8e-10 | 1.593724 (1.593724) | 1.588129 (1.588129) | 7.1e-7 |

Tolerances 2e-8 in E and 5e-4 in S̄. Every local solve met its residual target, and n_mid agrees to 1e-6 (0.218296, 0.201770). Wall time 697 s and 746 s on one GPU, against 2104 s and 2392 s for the pilot on the local machine.

G4, benchmark on L = 256 p153c (complex), bench.py (typhon_logs/bench_L256*.json). Passed, the speed-up at χ = 768 is 48 with the host SVD and 60 with the device SVD.

Method. Random Z₂-block MPS with every bond at min(χ, 2^distance) in mixed canonical form, exact environments, three consecutive mid-chain bond updates (TwoSite, preconditioner eigh, Davidson, split, environment update, entropy SVD) with 80 matvecs each. The CPU reference is the original z2dmrg.py on the same environments and tensors (one bond), 6 BLAS threads (as in production) on the GPU node's AMD EPYC 7742. The matvec of the two codes agrees to 1e-15 relative on the same vector. Seconds per bond update and per matvec.

| χ | GPU per bond | CPU per bond | speed-up per bond | GPU per matvec | CPU per matvec | speed-up per matvec |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| 32 | 0.62 | 0.12 | 0.19 | 0.0028 | 0.00056 | 0.2 |
| 64 | 0.55 | 0.33 | 0.61 | 0.0029 | 0.0017 | 0.6 |
| 128 | 0.57 | 1.48 | 2.6 | 0.0029 | 0.0076 | 2.6 |
| 192 | 0.63 | 3.29 | 5.2 | 0.0029 | 0.015 | 5.3 |
| 256 | 0.96 | 6.31 | 6.6 | 0.0031 | 0.029 | 9.3 |
| 384 | 0.73 | 15.5 | 21 | 0.0030 | 0.076 | 25 |
| 512 | 0.94 | 29.2 | 31 | 0.0029 | 0.149 | 51 |
| 768 | 1.57 (1.25 device SVD) | 74.9 | 48 (60) | 0.0037 | 0.415 | 111 |

- The GPU matvec costs about 3 ms at every χ up to 768. It is bound by kernel launches (a few dozen small block GEMMs per matvec), not by arithmetic, so the GPU gains grow with χ and the host is faster below χ = 128. Hence GPU_MIN_CHI = 128.
- At χ = 768 the GPU bond time is matvec 0.30 s, preconditioner 0.20 s, Davidson bookkeeping 0.37 s, preconditioner eigh 0.05 s, split 0.54 s with the host SVD (0.21 s on the device), entropy SVD 0.10 s (0.05 s on the device). SVD on a 768 × 768 complex block: cusolver 0.10 s, LAPACK 0.28 s, both with reconstruction error below 7e-15 and singular values equal to 6.6e-15 relative. Hence the device SVD in production.
- The CPU numbers match the CPU production logs: bench gives 3.29 s per bond at χ = 192 with 80 matvecs, about 1700 s per sweep, and the CPU runs p153c_L256_P0/P1 took 1764 s and 2686 s for their first χ = 192 sweep (95 matvecs on average).
- Device memory at χ = 768, L = 256 (full set of environments, MPS and Davidson basis with m ≤ 32): 27.2 GB in use of 42.4 GB (cupy pool 26.7 GB). One run per GPU.

Profile of the CPU code (step 2, local, 6 threads, L = 64, 60 matvecs per bond; PROGRESS 5). Per bond at χ = 256: matvec 24 %, preconditioner 5 %, Davidson bookkeeping 68 %, everything else 3 %. The bookkeeping is mostly the full recomputation of V^dag A V each iteration (a cost m²n). On the CPU node at χ = 768 the shares are matvec 44 %, preconditioner 15 %, bookkeeping 38 %.

## Registration amendment (written before any GPU production output)

The GPU runs are additional runs of the registered configurations (same model, field, schedule, seed and initialization as pipeline190, started from scratch, never from a CPU checkpoint). For each key, size and parity sector, the run with the lowest final energy among the runs (CPU or GPU) that pass C1 to C4 is used. Where a CPU and a GPU run of the same configuration agree in energy within 1e-7, their S̄₆ must agree within 0.01 (the existing rule C5), and otherwise the higher one is reported as metastable. Rules V1 to V6, the convergence targets C1 to C6 and the estimators are unchanged.

## Production (status 2026-11-09 22:14 EST)

Jobs. 17607561 (c190g_prod, partition apollo, QOS gpu-medium, node apollo02, 6 A100 GPUs, 28 cores, 320 GB, 2-day limit), started 21:48. 17607562 is the same script with --dependency=afterany:17607561. It resumes every unfinished run from its checkpoint (one sweep is lost at most) and skips finished ones. Script typhon/c190g_prod.slurm, settings C190G_GPU_MIN_CHI=128, C190G_SVD_DEV=1, 4 threads per run. Only one gpu-medium job of this user can run beside ddL_g, so the resume job waits for the first to end.

Order (pipeline190g.py, one run per GPU, the next task starts when a GPU frees). (a) p153c L = 256 P0 and P1 rand:1. (b) p153c L = 128 P0 and P1, rand:1 and prod:28/29, then L = 64 P0 and P1, rand:1 and prod:14/15. (c) sym40_n316 at its tuned field (−19.976047250613814 from outB/tune) at L = 256, 128, 64, rand:1, plus prod at L ≥ 128 as in the CPU pipeline. (d) xx1 and xx2 controls, skipped when the CPU run of the same tag is done at the moment a GPU frees. (e) rat_n316 (−19.975988986408243) and sym40_n18 (−19.976333693659278) were queued at start from the CPU tune files. p153_n316 and sym80_n316 are added when their tune files appear, as long as the pipeline is still running.

Output layout on Typhon, /scratch/lustre/jkudlerflam/c18x/c190g/: outG/<tag>.json (same format as outA), outG/logs/<tag>.log (JSON line per sweep, preceded by a "# c190g" line with backend, node and GPU), outG/cfg/<tag>.json, outG/ckpt/<tag>.npz (up to about 5 GB at L = 256, do not copy), logs/pipeline_G_<jobid>.log, logs/prod_smi_<jobid>.log (GPU memory and load every 5 minutes), logs/c190g_prod_<jobid>.out.

Collection, from inside codex/work_C190/gpu (tested on the partial output at 22:00, 156 kB):

    mkdir -p typhon_out/outG/logs
    scp 'typhon:/scratch/lustre/jkudlerflam/c18x/c190g/outG/*.json' typhon_out/outG/
    scp 'typhon:/scratch/lustre/jkudlerflam/c18x/c190g/outG/logs/*.log' typhon_out/outG/logs/
    scp 'typhon:/scratch/lustre/jkudlerflam/c18x/c190g/logs/pipeline_G_*.log' typhon_out/
    python3 ../analyze190.py typhon_out/outG --partial

analyze190.py reads outG unchanged (it ran on the partial output and printed the header with no rows yet, since no run had a mixer-off sweep). It analyzes one folder at a time. The selection of the amendment (lowest energy among CPU and GPU runs of a tag) needs the two folders side by side, and no merge script was written.

First-sweep check against the CPU runs (outA/logs, read only). Every local solve of every sweep so far meets its residual target (res_max ≤ res_tol). Energies E_phys at the same sweep:

| run | sweep (χ) | GPU | CPU | difference |
| :-- | :-- | :-- | :-- | :-- |
| p153c_L256_P0_rand1 | 3 (64) | −2570.2565490656 | −2570.2565489877 | 7.8e-8 |
| p153c_L256_P0_rand1 | 4 (128) | −2570.2569873415 | −2570.2569873075 | 3.4e-8 |
| p153c_L256_P0_rand1 | 6 (192) | −2570.2570582039 | −2570.2570582777 | 7.4e-8 |
| p153c_L256_P1_rand1 | 5 (128) | −2570.2569788372 | −2570.2569785149 | 3.2e-7 |
| p153c_L256_P1_rand1 | 6 (192) | −2570.2570630230 | −2570.2570629764 | 4.7e-8 |
| p153c_L128_P0_rand1 | 7 (192) | −1290.1134051929 | −1290.1134051345 | 5.8e-8 |
| p153c_L128_P0_rand1 | 8 (256) | −1290.1134126066 | −1290.1134126251 | 1.9e-8 |
| p153c_L128_P1_rand1 | 8 (256) | −1290.1134079109 | −1290.1134079099 | 1.0e-9 |
| p153c_L128_P1_rand1 | 9 (256) | −1290.1134075302 | −1290.1134075210 | 9.2e-9 |
| p153c_L128_P0_prod28 | 8 (256) | −1290.1134126249 | −1290.1134126228 | 2.1e-9 |
| p153c_L128_P0_prod28 | 9 (256) | −1290.1134121865 | −1290.1134121773 | 9.2e-9 |
| p153c_L128_P1_prod29 | 5 (128) | −1290.1133655558 | −1290.1133655547 | 1.1e-9 |

At the first mixer-off sweep (χ = 256, L = 128) S̄₆ agrees to 1.3e-5 (P1 rand1, 1.832528 against 1.832541) and 2.0e-5 (P0 prod28, 1.839789 against 1.839769).

The very first sweeps (χ = 32, residual 1e-7, from the random start) differ by up to 5e-4 (p153c_L256_P1 sweep 0). This is the CPU code's own sensitivity to the BLAS build, not the port. On apollo02 the original run190.py on the first four sweeps of p153c_L128_P0_rand1 with 4 and with 6 BLAS threads gives energies equal to the GPU pipeline's host stage in all printed digits (−1290.0992102557, −1290.1120509129, −1290.1132150438, −1290.1131919250), while the CPU production run on its own node gives −1290.0992279136 at sweep 0 (1.8e-5 away) and the same code with 1 thread −1290.0992104585. From χ = 128 on, the GPU and CPU runs agree to 3e-8 to 1.2e-6 in E within the mixer stages (the largest, 1.2e-6, at p153c_L256_P1 sweep 4) and to 1e-9 to 2e-8 once the mixer is off.

Measured sweep times on the GPU (CPU production in brackets): L = 128, χ = 128: 151 to 167 s (308 to 421 s); χ = 192: 147 to 166 s (724 to 942 s); χ = 256: 129 to 143 s (1207 to 1354 s); χ = 384: 166 to 179 s (not reached on the CPU). L = 256, χ = 64 on the host: 126 to 223 s (124 to 179 s); χ = 128: 300 to 327 s (570 to 863 s); χ = 192: 332 to 337 s (1764 to 2686 s). Device memory per run at L = 128, χ = 384: 6 to 10 GB.

Estimate for the L = 256 runs. They finished their first χ = 192 sweep at 22:13. The remaining pre-final sweeps (one at 192, two each at 256, 384 and 512) should take about 45 to 50 minutes, extrapolated from the measured L = 128 sweep times at the same χ (per bond 0.55 s at 256 and 0.69 s at 384, times 510 bond updates) and the bench ratio 1.3 between χ = 512 and 384. That puts the start of the final χ = 768 stage at about 23:00 to 23:15 on 2026-11-09. The bench ratio 1.7 between 768 and 384 gives about 10 minutes per final sweep, so the 2 to 6 final sweeps end between about 23:20 and 00:15. The final host measurements (entropies, one-body matrices, 4 threads) at χ = 768 are not measured and may add tens of minutes. This is an extrapolation, uncertain by perhaps a factor 1.5, from per-bond times that do not yet include a sweep at χ ≥ 512 in production.

## R5. p153c at L = 64 and 128, registered final χ (collected 2026-11-09 22:57 EST)

All eight GPU runs finished (two initializations per parity sector and size, χ = 384 at L = 64 and 512 at L = 128). Every run passes C1 to C4. The two initializations agree to 1e-9 in energy and to four decimals in S̄ (C5). The one CPU run finished so far at a final χ, p153c L = 64 even, has E = −650.0414853399 and S̄₁ = 1.61456, against −650.0414853398 and 1.6146 here.

| L | sector | E | S̄₁ | S̄₆ | S̄₁₂ | S̄₆ step from previous χ | n_mid | W_min/W_max |
|---|---|---|---|---|---|---|---|---|
| 64 | even | −650.0414853399 | 1.6146 | 1.6090 | 1.5933 | 0.0005 | 0.2183 | 0.783 |
| 64 | odd | −650.0414362973 | 1.5941 | 1.5885 | 1.5727 | 0.0004 | 0.2018 | 0.751 |
| 128 | even | −1290.1134157022 | 1.8682 | 1.8591 | 1.8515 | 0.0046 | 0.2211 | 0.779 |
| 128 | odd | −1290.1134107394 | 1.8598 | 1.8510 | 1.8436 | 0.0044 | 0.2130 | 0.764 |

The even sector is the ground branch at both sizes. c_two(64 → 128) on the ground branch is 2.196, 2.165, 2.234 for W = 1, 6, 12, and 2.12 to 2.37 over all parity pairings at W = 1. By V3 the collaborator's downturn does not reproduce (the threshold is 1.6). L·Δ_GC/|t₂| is +0.251 at L = 64 and +0.051 at L = 128. The valley weight ratio at L = 128 is 0.78 and 0.76, a little under the 0.8 of V6. The analysis output is `typhon_out/analyze_G_partial.txt`. The L = 256 runs were at χ = 384 (even) and 512 (odd) of 768 at collection, and no reading is taken from them.

## R6. L = 256 at χ = 768: p153c odd sector and sym40_n316 (collected 2026-11-10 00:45 EST)

Runs with status done. The analysis is in typhon_out/analyze_G_partial.txt.

| key | L | sector | init | E | S̄₁ | S̄₆ | S̄₁₂ | c_prof | n_mid | S̄₆ step from χ = 512 | C1 to C4 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| p153c | 256 | odd | rand1 | −2570.2571184393 | 2.0879 | 2.0827 | 2.0791 | 2.028 | 0.2153 | 0.0269 | C4 fails |
| sym40_n316 | 256 | even | rand1 | −2570.2770003174 | 2.0633 | 2.0579 | 2.0547 | 2.032 | 0.1956 | 0.0249 | C4 fails |
| sym40_n316 | 256 | even | prod50 | −2570.2770003188 | 2.0633 | 2.0579 | 2.0547 | 2.032 | 0.1956 | 0.0246 | C4 fails |
| sym40_n316 | 256 | odd | rand1 | −2570.2770066290 | 2.0606 | 2.0558 | 2.0527 | 2.045 | 0.1917 | 0.0237 | C4 fails |
| sym40_n316 | 256 | odd | prod51 | −2570.2770066288 | 2.0606 | 2.0558 | 2.0527 | 2.045 | 0.1917 | 0.0237 | C4 fails |
| sym40_n316 | 128 | even | rand1 | −1290.1231026155 | 1.8265 | 1.8183 | 1.8112 | 2.251 | 0.1889 | 0.0034 | pass |

Every run passes C1 (residual 1e-9), C2 (energy steps 5e-10 to 1.4e-9, entropy steps 1e-6 to 1e-5) and C3 (discarded weight 3e-10 to 3.5e-10). C5 passes for sym40_n316 at L = 256 in both sectors (energies agree to 1.5e-9, S̄₆ to 1e-4). C6 passes for sym40_n316 (n_mid 0.1917 on the ground branch against the target 0.1957).

S̄₆ by stage, last sweep with the mixer off, with the discarded weight.

| run | χ = 384 | χ = 512 | χ = 768 |
|---|---|---|---|
| p153c L = 256 odd | 2.0203 (2.6e-9) | 2.0559 (1.2e-9) | 2.0828 (3.5e-10) |
| sym40_n316 L = 256 odd | 2.0013 (2.2e-9) | 2.0322 (1.0e-9) | 2.0559 (2.9e-10) |
| sym40_n316 L = 256 even | 1.9961 (2.3e-9) | 2.0331 (1.1e-9) | 2.0579 (3.0e-10) |

Readings.
- p153c, c_two(128 → 256) on the ground branch (even at L = 128, odd at L = 256) is 1.902, 1.936, 1.971 for W = 1, 6, 12. The even-sector run at L = 256 is still running at χ = 768, with energy −2570.2571001943 (1.8e-5 above the odd sector) and S̄₆ = 2.0830. With it the parity pairings at W = 1 span 1.90 to 1.99.
- p153c, parity splitting L·Δ_GC/|t₂| = +0.251, +0.051, −0.374 at L = 64, 128, 256 (the last with the unfinished even run).
- p153c, effective-model offset ΔS̄₆ = −0.2315, −0.2251, −0.2129 at L = 64, 128, 256. V4's difference is +0.019, inside ±0.05.
- p153c, valley ratio W_min/W_max = 0.772 at L = 256, under V6's 0.8. On the symmetric line it is 0.999.
- sym40_n316, c_two(128 → 256) from the two finished even-sector runs is 2.049, 2.074, 2.108 for W = 1, 6, 12. The odd sector is the ground branch at L = 256 (L·Δ_GC/|t₂| = −0.129), and its L = 128 partner is still at χ = 256, so the ground-branch value is not available yet. Effective-model offset at L = 256 is −0.2048.

Against the registration.
- V2 (falls toward 1) needs c_two(128 → 256) ≤ 1.3 at W = 6 and 12, or c_two(64 → 128) ≤ 1.3. The readings are 1.94 and 1.97, and 2.17 and 2.23. V2 does not hold.
- V1's thresholds hold: c_two ≥ 1.6 at W = 6 and 12 for both size pairs, and c_prof(256) = 2.03. V1 also requires C1 to C4 for every run and the two-leg control at the same size pairs. C4 fails at L = 256 for every run (S̄₆ rises by 0.024 to 0.027 from χ = 512 to 768, against 0.01), and the L = 256 two-leg control has not finished. V1 is not declared.
- As registered, the outcome so far is "unresolved", for the reason that the L = 256 entropy is not converged in χ.

Direction of the bias. The entropy rises with χ at fixed size, and the L = 128 partner passes C4 (step 0.0046). An unconverged L = 256 state therefore lowers c_two. Outside the registration, a linear extrapolation in the discarded weight through the last two stages gives S̄₆ = 2.094 at L = 256 and 1.861 at L = 128, and c_two = 2.01.

## R7. p153c even sector at L = 256 finished (collected 2026-11-10 00:50 EST)

p153c_L256_P0_rand1, χ = 768, 20 sweeps. E = −2570.2571001943, S̄₁ = 2.0892, S̄₆ = 2.0830, S̄₁₂ = 2.0793, n_mid = 0.2192, c_prof = 2.011, W_a/W_b = 0.779. C1 to C3 pass (residual 1e-9, last two sweeps |ΔE| = 1e-9, |ΔS̄₆| = 1.3e-4, |Δn_mid| = 2e-6, discarded weight 3.7e-10). C4 fails, step 0.0871 from χ = 512.

| χ | S̄₆ even | S̄₆ odd | discarded weight (even, odd) |
|---|---|---|---|
| 384 | 1.9122 | 2.0203 | 2.1e-9, 2.6e-9 |
| 512 | 1.9959 | 2.0559 | 1.1e-9, 1.2e-9 |
| 768 | 2.0830 | 2.0828 | 3.7e-10, 3.5e-10 |

The sectors agree at χ = 768 to 2e-4 in S̄₆ after different histories. Ground sector odd, L·Δ_GC/|t₂| = −0.374. Ground-branch c_two(128 → 256) = 1.902, 1.936, 1.971 (unchanged from R6), all parity pairings at W = 1 from 1.90 to 1.99. The verdict of R6 stands: V2 does not hold, and V1 waits on C4 at L = 256 and on the two-leg control. Updated table in typhon_out/analyze_G_partial.txt.

## R8. Registration of a χ = 1024 stage at L = 256 (written 2026-11-10 03:10 EST, before any χ = 1024 output)

Coordinator, Claude Fable 5. Jonah approved this stage on 2026-11-10. It is written after the χ = 768 results of R6 and R7 were seen, and before any χ = 1024 output exists. No criterion of `../RESULT.md` changes. The L = 256 runs are extended by one stage, because C4 is the one criterion that fails at χ = 768.

Why C4 fails at χ = 768. The stages at χ = 384 and 512 have two sweeps each and are not converged in sweeps. In the even sector the first sweep at χ = 768 alone moves S̄₆ from 1.9959 to 2.0582, and the six sweeps at χ = 768 then read 2.05815, 2.07751, 2.08143, 2.08254, 2.08290, 2.08303. So the registered step from χ = 512 mixes the sweep convergence of the earlier stage with the dependence on χ. A stage at χ = 1024 that starts from the converged χ = 768 state measures the second one alone.

Runs. Each continues a finished run of a registered configuration from its final χ = 768 state (the checkpoint of that run, read only), with the same model, field, gauge, parity sector, Davidson settings and tolerances, and sweeps at χ = 1024 with residual target 1e-9 and no subspace expansion. The stop rule is that of the registered final stage (|ΔE| < 2e-9 and |ΔS̄₆| < 5e-4 between the last two sweeps, at most 6 sweeps), with at least 3 sweeps instead of 2, since the first sweep grows the bonds.

| run | continued from | χ = 768 values (E, S̄₁, S̄₆, S̄₁₂, discarded weight) |
|---|---|---|
| p153c_L256_P0_rand1 | 20 sweeps | −2570.2571001943, 2.0892, 2.0830, 2.0793, 3.7e-10 |
| p153c_L256_P1_rand1 | 18 sweeps | −2570.2571184393, 2.0880, 2.0828, 2.0792, 3.5e-10 |
| sym40_n316_L256_P0_rand1 | 18 sweeps | −2570.2770003174, 2.0633, 2.0579, 2.0547, 3.0e-10 |
| sym40_n316_L256_P1_rand1 | 18 sweeps | −2570.2770066290, 2.0607, 2.0559, 2.0527, 2.9e-10 |
| xx2_L256_P0_rand1 (two-leg control) | 20 sweeps | −12.2485203313, 2.0307, 2.0104, 2.0158, 2.0e-9 (from the run's own record, output not collected yet) |

For sym40_n316 the two initializations of each sector agree at χ = 768 (energies within 2e-9, S̄₆ within 1e-4, so C5 passes), and one per sector is continued. The control is continued so that candidate and control are read at the same χ.

How the criteria are read.
- C1 to C3 apply to the χ = 1024 sweeps as written.
- C4 at L = 256 is |S̄₆(1024) − S̄₆(768)| ≤ 0.01, with each value taken at the last sweep of its stage. This is the registered rule with χ = 1024 as the final stage.
- The run used for a key and sector is the one of lowest final energy among those passing C1 to C4, as in the amendment above. A continuation has lower energy than the run it starts from, so it is the run used if it passes.
- c_two(128 → 256) pairs the χ = 1024 value at L = 256 with the L = 128 run at its registered final χ = 512. c_two(64 → 128) is unchanged.
- V1 and V2 are evaluated as written. The control condition of V1 at 128 → 256 is evaluated with the control at χ = 1024, and the χ = 768 value is quoted next to it.
- If C4 still fails at χ = 1024 for a run that V1 needs, V1 is not declared from these runs. The readings at χ = 768 and 1024 and the extrapolation below are then reported as they are, without a verdict.

Extrapolation in the discarded weight (reported next to the registered reading, not used for V1 or V2). With ε the maximum discarded weight of the last sweep of a stage,

S̄_W(ε → 0) = S̄_W(1024) + [S̄_W(1024) − S̄_W(768)] ε₁₀₂₄ / (ε₇₆₈ − ε₁₀₂₄),

quoted with an uncertainty equal to the whole increment |S̄_W(ε → 0) − S̄_W(1024)|, and not quoted if ε₁₀₂₄ ≥ ε₇₆₈. Only the two stages that are converged in sweeps enter. The L = 128 partner is not extrapolated, and its own C4 step is quoted as the size of what that neglects. Since the entropy grows with χ at both sizes, the extrapolated c_two is an upper estimate. A shift δ of S̄_W at L = 256 alone moves c_two(128 → 256) by 6δ/ln 2 = 8.66 δ.

Code and job. `pipeline_x1024.py` (one run per GPU, output folder `outX1024/` on Typhon scratch with the layout of `outG/`, nothing in `outG/` is written), `typhon/c190g_x1024.slurm` (partition h200, 3 GPUs, 12 cores, 8 h, cupy smoke test on the h200 node before production, since cupy was tested on apollo02 only). The complex runs need about 48 GB of device memory at χ = 1024, which the 40 GB cards of apollo do not have. A resubmission resumes unfinished runs from `outX1024/ckpt` and skips finished ones.

Submitted 2026-11-10 03:20 EST as job 17609630 (c190g_x1024, QOS gpu-h200). The dry listing of `pipeline_x1024.py` on Typhon shows the five runs with the χ = 768 values above. The QOS allows one running job per user on h200, and another of Jonah's jobs holds that slot until 10:16 EST at the latest, so the job is pending with reason QOSMaxJobsPerUserLimit. Estimated run time once started, from the χ = 768 sweep times on the A100 cards scaled by χ³, is at most about 5 h for the five runs.

Collection, once runs are done (status done only, no checkpoints): `scp 'typhon:/scratch/lustre/jkudlerflam/c18x/c190g/outX1024/*.json' typhon_out/outX1024/` and `scp 'typhon:/scratch/lustre/jkudlerflam/c18x/c190g/outX1024/logs/*.log' typhon_out/outX1024/logs/`.

## R9. Controls at χ = 768 collected: both calibrate both size pairs (2026-11-10 03:45 EST)

GPU runs xx2 at L = 64, 128, 256 and xx1 at L = 64, 128, all with status done, in typhon_out/outG/. The XX chain at L = 256 is the CPU run of node B (`../typhon_out/outB/xx1_L256_P0_rand1.json`).

| run | final χ | E | S̄₁, S̄₆, S̄₁₂ | C4 step | discarded weight | last ΔE |
|---|---|---|---|---|---|---|
| xx2_L64_P0 | 384 | −2.8537224973 | 1.4960, 1.5699, 1.5273 | 0.0000 | 2.1e-11 | 3.7e-10 |
| xx2_L128_P0 | 512 | −6.9480312045 | 1.8304, 1.7913, 1.7998 | 0.0003 | 9.9e-10 | 6.1e-9 |
| xx2_L256_P0 | 768 | −12.2485203313 | 2.0307, 2.0104, 2.0158 | 0.0026 | 2.0e-9 | 3.2e-8 |
| xx1_L64_P1 | 384 | −3.4740158336 | 0.9103, 0.8997, 0.8889 | 0.0000 | 5.8e-14 | 1e-14 |
| xx1_L128_P1 | 512 | −6.1242613702 | 1.0132, 1.0082, 1.0050 | 0.0000 | 1.0e-13 | 2e-14 |
| xx1_L256_P0 (CPU) | 768 | | 1.1186, 1.1206, 1.1207 | 0.0000 | 1.5e-13 | 2e-13 |

| c_two at W = 1, 6, 12 | 64 → 128 | 128 → 256 |
|---|---|---|
| XX chain (c = 1) | 0.891, 0.939, 1.005 | 0.912, 0.974, 1.002 |
| two-leg XX (c = 2) | 2.895, 1.917, 2.359 | 1.734, 1.896, 1.869 |
| p153c ground branch, for comparison (R5 to R7) | 2.196, 2.165, 2.234 | 1.902, 1.936, 1.971 |

- The registered calibration is at W = 6 (XX within 1 ± 0.15, two-leg within 2 ± 0.25). It holds for both size pairs: 0.939 and 0.974, 1.917 and 1.896.
- The two-leg control depends strongly on the window at 64 → 128 (2.90 at W = 1, 2.36 at W = 12), which is the oscillation of a free two-band entropy with the cut position. The candidate does not show a spread of that size.
- The two-leg run at L = 256 agrees with the exact free-fermion values quoted in `../RESULT.md` (2.0311, 2.0108, 2.0162) to 4e-4. Its residual, entropy steps, density step, discarded weight and C4 pass. Its last energy step, 3.2e-8 after the maximum of six sweeps, is above the 1e-8 of C2. The χ = 1024 continuation of R8 adds sweeps to this run.
- With this, V1 at L = 256 waits on C4 for the candidate runs, and on the energy step of the control.

## R10. Two more tuned keys at χ = 768: rat_n316 and sym40_n18 (collected 2026-11-10 03:40 EST)

GPU runs of the collaborator's rational bond tuned to n = 3/16 (rat_n316, h_z = −19.975988986) and of the symmetric-line bond tuned to n = 1/8 (sym40_n18, h_z = −19.976333694), at L = 64, 128 and 256, both parity sectors, initialization rand1, all with status done. Files in typhon_out/outG/. The analysis is typhon_out/analyze_G_rat_sym18.txt, made with `python3 -E -s analyze_keys.py outG` from typhon_out/ (the script is `../analyze_keys.py` there). It applies C1 to C4 and C6 and computes c_two for the ground branch and for each sector.

| key | L | sector | final χ | E | S̄₁, S̄₆, S̄₁₂ | n_mid | discarded weight | C4 step | C1 to C4 |
|---|---|---|---|---|---|---|---|---|---|
| rat_n316 | 64 | even | 384 | −650.0447774541 | 1.5732, 1.5658, 1.5494 | 0.1852 | 1.5e-11 | 0.0003 | pass |
| rat_n316 | 64 | odd | 384 | −650.0448097690 | 1.5957, 1.5897, 1.5736 | 0.2016 | 2.4e-11 | 0.0004 | pass |
| rat_n316 | 128 | even | 512 | −1290.1204465833 | 1.8505, 1.8415, 1.8343 | 0.2048 | 2.2e-10 | 0.0041 | pass |
| rat_n316 | 128 | odd | 512 | −1290.1204594204 | 1.8404, 1.8325, 1.8255 | 0.1967 | 2.0e-10 | 0.0037 | pass |
| rat_n316 | 256 | even | 768 | −2570.2715836167 | 2.0675, 2.0629, 2.0597 | 0.1955 | 3.0e-10 | 0.0247 | C4 fails |
| rat_n316 | 256 | odd | 768 | −2570.2715871821 | 2.0720, 2.0673, 2.0641 | 0.1995 | 3.2e-10 | 0.0252 | C4 fails |
| sym40_n18 | 64 | even | 384 | −650.0571061504 | 1.3738, 1.3874, 1.3815 | 0.1230 | 1.4e-12 | 0.0000 | pass |
| sym40_n18 | 64 | odd | 384 | −650.0571528391 | 1.4747, 1.4705, 1.4531 | 0.1348 | 2.2e-12 | 0.0001 | pass |
| sym40_n18 | 128 | even | 512 | −1290.1464438772 | 1.7465, 1.7373, 1.7305 | 0.1393 | 6.1e-11 | 0.0016 | pass |
| sym40_n18 | 128 | odd | 512 | −1290.1464423236 | 1.7221, 1.7158, 1.7102 | 0.1315 | 4.8e-11 | 0.0013 | pass |
| sym40_n18 | 256 | even | 768 | −2570.3249046022 | 1.9591, 1.9531, 1.9502 | 0.1324 | 1.1e-10 | 0.0149 | C4 fails |
| sym40_n18 | 256 | odd | 768 | −2570.3249096989 | 1.9694, 1.9658, 1.9636 | 0.1359 | 1.2e-10 | 0.0140 | C4 fails |

Every run meets C1 to C3 (residual 1e-9, last energy step at most 1.8e-9, entropy step at most 3.7e-5). The L = 256 runs have 18 sweeps, 4 at χ = 768.

| c_two at W = 1, 6, 12 | 64 → 128 | 128 → 256 |
|---|---|---|
| rat_n316, ground branch (odd at all three sizes) | 2.119, 2.102, 2.180 | 2.004, 2.033, 2.065 |
| rat_n316, even sector | 2.401, 2.386, 2.466 | 1.878, 1.917, 1.952 |
| sym40_n18, ground branch (odd, even, odd) | 2.353, 2.309, 2.401 | 1.930, 1.978, 2.018 |
| sym40_n18, even sector | 3.226, 3.028, 3.021 | 1.840, 1.868, 1.902 |
| sym40_n18, odd sector | 2.141, 2.124, 2.225 | 2.141, 2.164, 2.193 |
| p153c ground branch (R5 to R7) | 2.196, 2.165, 2.234 | 1.902, 1.936, 1.971 |
| sym40_n316 ground branch (R6) | | 2.049, 2.074, 2.108 |

- Sector splittings E(even) − E(odd): rat_n316 +3.23e-5, +1.28e-5, +3.57e-6 at L = 64, 128, 256, and sym40_n18 +4.67e-5, −1.55e-6, +5.10e-6. For sym40_n18 the ground sector alternates, so its ground branch mixes sectors and the per-sector rows are the cleaner comparison.
- C6 passes for both keys. The ground-sector n_mid at L = 256 is 0.1995 against the target 0.1957, and 0.1359 against 0.1326.
- Against the registration. V2 is excluded for both keys, since c_two(128 → 256) is 1.87 to 2.19 at W = 6 and 12 in every sector, and c_two(64 → 128) is above 2. The values meet V1's thresholds, but V1 is not declared, because C4 fails at L = 256 (steps 0.025 for rat_n316 and 0.014 to 0.015 for sym40_n18 against the allowed 0.01). This is the same pattern as p153c and sym40_n316.
- Four keys now read about 2 from 128 to 256 at χ = 768, at two densities (n = 3/16 and 1/8) and at three bonds (the collaborator's field, their rational bond, the symmetric line). The n = 1/8 point has the smaller C4 step.
- These two keys are not in the χ = 1024 job of R8. A χ = 1024 stage for them would need its own registration and submission.

## R11. Production complete: p153_n316 and sym80_n316 at χ = 768, and the CPU–GPU cross-check (collected 2026-11-10 06:40 EST)

The production job has finished, and all 51 runs in typhon_out/outG/ have status done. The last two keys are the collaborator's φ = 1.53 bond tuned to n = 3/16 (p153_n316, h_z = −19.975988793) and the symmetric-line bond at V/Δ = 80 tuned to n = 3/16 (sym80_n316, h_z = −39.988023246). The analysis is typhon_out/analyze_G_p153n316_sym80.txt, made with `python3 -E -s ../analyze_keys.py outG p153_n316:0.1957 sym80_n316:0.1957` from typhon_out/. The script now takes the keys as arguments (with the R10 keys as the default), lists every finished run of a configuration, and uses the lowest-energy one among those passing C1 to C4, as the amendment says.

| key | L | sector | final χ | E | S̄₁, S̄₆, S̄₁₂ | n_mid | discarded weight | C4 step | C1 to C4 |
|---|---|---|---|---|---|---|---|---|---|
| p153_n316 | 64 | even | 384 | −650.0447248647 | 1.5735, 1.5659, 1.5495 | 0.1851 | 1.5e-11 | 0.0003 | pass |
| p153_n316 | 64 | odd | 384 | −650.0447490177 | 1.5954, 1.5894, 1.5734 | 0.2016 | 2.4e-11 | 0.0004 | pass |
| p153_n316 | 128 | even | 512 | −1290.1203340670 | 1.8493, 1.8399, 1.8327 | 0.2047 | 2.2e-10 | 0.0042 | pass |
| p153_n316 | 128 | odd | 512 | −1290.1203554204 | 1.8403, 1.8324, 1.8254 | 0.1967 | 2.0e-10 | 0.0037 | pass |
| p153_n316 | 256 | even | 768 | −2570.2713859128 | 2.0675, 2.0630, 2.0598 | 0.1955 | 3.0e-10 | 0.0247 | C4 fails |
| p153_n316 | 256 | odd | 768 | −2570.2713807247 | 2.0718, 2.0670, 2.0638 | 0.1995 | 3.1e-10 | 0.0253 | C4 fails |
| sym80_n316 | 64 | even | 384 | −1300.0228725231 | 1.5741, 1.5701, 1.5527 | 0.1850 | 9.9e-12 | 0.0006 | pass |
| sym80_n316 | 64 | odd | 384 | −1300.0228656912 | 1.5961, 1.5901, 1.5741 | 0.2016 | 1.5e-11 | 0.0007 | pass |
| sym80_n316 | 128 | even | 512 | −2580.0612141617 | 1.8261, 1.8190, 1.8120 | 0.1884 | 9.7e-11 | 0.0050 | pass |
| sym80_n316 | 128 | odd | 512 | −2580.0612270693 | 1.8382, 1.8305, 1.8235 | 0.1965 | 1.1e-10 | 0.0057 | pass |
| sym80_n316 | 256 | even | 768 | −5140.1378350102 | 2.0594, 2.0548, 2.0516 | 0.1952 | 1.5e-10 | 0.0346 | C4 fails |
| sym80_n316 | 256 | odd | 768 | −5140.1378342986 | 2.0548, 2.0505, 2.0474 | 0.1912 | 1.5e-10 | 0.0339 | C4 fails |

Every run meets C1 to C3 (residual 1e-9, last energy step at most 1.9e-9, entropy step at most 4.7e-5, discarded weight at most 3.1e-10).

| c_two at W = 1, 6, 12 | 64 → 128 | 128 → 256 |
|---|---|---|
| p153_n316, ground branch (odd, odd, even) | 2.120, 2.104, 2.182 | 1.966, 1.996, 2.029 |
| p153_n316, even sector | 2.388, 2.371, 2.451 | 1.889, 1.931, 1.966 |
| p153_n316, odd sector | 2.120, 2.104, 2.182 | 2.003, 2.030, 2.063 |
| sym80_n316, ground branch (even, odd, even) | 2.286, 2.253, 2.344 | 1.915, 1.942, 1.975 |
| sym80_n316, even sector | 2.181, 2.154, 2.245 | 2.020, 2.042, 2.074 |
| sym80_n316, odd sector | 2.096, 2.081, 2.159 | 1.876, 1.904, 1.938 |

- Sector splittings E(even) − E(odd): p153_n316 +2.42e-5, +2.14e-5, −5.19e-6 at L = 64, 128, 256, and sym80_n316 −6.83e-6, +1.29e-5, −7.1e-7. The ground sector changes with L for both keys, so the per-sector rows are the cleaner comparison. At L = 256 the sym80_n316 sectors are 7e-7 apart, which is resolved against the 2e-9 sweep steps but small.
- Two initializations. For p153_n316 at L = 128 the product start and the random start end at the same state in both sectors. The energies agree to 1e-10 and S̄₆ to 1e-4, so C5 holds there.
- C6 passes for both keys, with ground-sector n_mid at L = 256 of 0.1955 and 0.1952 against the target 0.1957.
- Against the registration. V2 is excluded for both keys, since c_two(128 → 256) is 1.90 to 2.07 at W = 6 and 12 in every sector, and c_two(64 → 128) is above 2. The values meet V1's thresholds, but V1 is not declared, because C4 fails at L = 256 (steps 0.025 for p153_n316 and 0.034 to 0.035 for sym80_n316 against the allowed 0.01).
- All six keys now read between 1.94 and 2.03 at W = 6 on the ground branch from 128 to 256 (p153c 1.936, sym40_n316 1.946, rat_n316 2.033, sym40_n18 1.978, p153_n316 1.996, sym80_n316 1.942). The reading does not move between V/Δ = 40 and 80, where the corrections to the effective model differ by a factor of four. The C4 step is largest at V/Δ = 80.

CPU against GPU (rule C5). Three CPU runs of p153c have finished with the same schedule as their GPU counterparts. The comparison is typhon_out/cross_cpu_gpu.txt, made with `python3 -E -s cross_cpu_gpu.py ../typhon_out/outA typhon_out/outG` from codex/work_C190/gpu/.

| run | χ | E (CPU) | E (GPU) | E difference | S̄₁, S̄₆, S̄₁₂ (both) | S̄₆ difference |
|---|---|---|---|---|---|---|
| p153c, L = 64, even, rand1 | 384 | −650.0414853399 | −650.0414853398 | 4e-12 | 1.6146, 1.6090, 1.5933 | 1.4e-7 |
| p153c, L = 128, even, prod28 | 512 | −1290.1134157016 | −1290.1134157013 | 3.3e-10 | 1.8682, 1.8591, 1.8515 | 1.4e-7 |
| p153c, L = 128, odd, rand1 | 512 | −1290.1134107386 | −1290.1134107388 | 1.8e-10 | 1.8598, 1.8510, 1.8436 | 8.7e-7 |

The two codes land on the same states. The energies agree far inside C5's 1e-7 and S̄₆ agrees to 1e-6 against the allowed 0.01. The CPU runs at L = 256 are still going, so the L = 256 readings rest on the GPU runs alone for now.

What is left for the registration: the χ = 1024 stage of R8 for C4 at L = 256 (queued, p153c and sym40_n316 only), and the CPU L = 256 runs as a further C5 check.

## Ground-branch summary correction (OpenAI Codex audit, 2026-11-10)

The completed sym40_n316 files select odd parity at both 128 and 256 sites. Its ground-branch W=6 slope is 1.946, not 2.074; the latter is the even-sector comparison from the earlier partial collection. R11’s current six-key summary is corrected to 1.94–2.03. Historical partial-stage tables are retained. C4 still fails and no registered pass is declared. See [reproducible audit](../../../contrib/d2_status_audit_20261110/README.md).
