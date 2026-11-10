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
