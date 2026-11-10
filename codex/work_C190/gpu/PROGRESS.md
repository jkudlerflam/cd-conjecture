This work is by a Claude Opus 5.5 agent (C190g, 2026-11-09), not by Codex. Astra was not used (out of credits).

# C190g progress log (GPU backend for the C190 DMRG)

Spec: specs/C190g_gpu_port.md. Local work in codex/work_C190/gpu/, Typhon work in /scratch/lustre/jkudlerflam/c18x/c190g/. Everything else is read-only.

1. Read the spec, z2dmrg.py, run190.py, pipeline190.py, models.py, validate.py, analyze190.py (reads outdir/*.json and outdir/tune/*.json, so a GPU output folder works unchanged), RESULT.md and PROOF section 5.

2. Wrote z2dmrg_xp.py (array-module port: every sweep routine runs in the module of its inputs, numpy or cupy; Hamiltonian classes and final measurements are imported unchanged from z2dmrg.py), run190g.py (run190.py with the backend switch, same JSON and checkpoint format), validate_xp.py (gates G1, G2) and bench.py (profile and gate G4).

3. G1 locally (logs/G1_numpy.log, 4 threads). At L = 10, all six cases, both sectors, both starts, the numpy path of z2dmrg_xp.py is bit-identical to z2dmrg.py (every difference 0.0) and its errors against ED are those of validate.py (worst candidate-bond errors 3.1e-12 in E, 6.3e-10 in S; the xx2 odd row is the known degeneracy). The L = 32 run (p153c, chi 32, 64, 96 with mixer then without, 6 sweeps) agrees sweep by sweep to 1.8e-11 in E and 1.1e-11 in the entropy at every cut, with identical matvec counts. That residual difference is the run-to-run noise of the original code itself: running z2dmrg.py twice on the same input gives energies differing by 4.3e-11 after the first mixer sweep (multithreaded BLAS is not bitwise reproducible). G1 passes.

4. Typhon folder c190g/ created. z2dmrg.py, models.py, validate.py, pipeline190.py, analyze190.py, run190.py and c176_ref.json copied from c190/ (cp -p, md5 identical to the local copies).

5. Profile of the CPU path (step 2), local, 6 threads, bench.py numpy at L = 64 (p153c), two mid-chain bonds per chi, 60 matvecs per solve, random mixed-canonical state with exact environments (logs/profile_L64_numpy.log). Seconds per bond, mean of the two bonds.

   | chi | total | matvec | preconditioner | Davidson bookkeeping | kron4 eigh | split SVD | env update | entropy SVD |
   | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
   | 192 | 1.73 | 0.38 (22 %) | 0.085 (5 %) | 1.23 (71 %) | 0.014 | 0.018 | 0.006 | 0.005 |
   | 256 | 3.56 | 0.87 (24 %) | 0.17 (5 %) | 2.43 (68 %) | 0.028 | 0.041 | 0.011 | 0.010 |

   The Davidson bookkeeping, not the matvec, dominates on the CPU. It is mostly the full recomputation of the m x m Rayleigh-Ritz matrix V^dag A V every iteration (m up to 32, a cost m^2 n, plus a conjugated copy of V each time) and the two classical Gram-Schmidt passes. That cost grows like chi^2 times m^2, which is the part of the measured chi^2-to-chi^3 growth outside the contractions. On the GPU these are large GEMMs and GEMVs and move with everything else. The algorithm is kept unchanged so that the numbers stay those of the CPU code (an incremental Ritz matrix would change them at roundoff level).

6. First test job (17607414, one 4-GPU job, because the gpu-short QOS allows two jobs per user and Jonah's gpu_lean1 and gpu_lean2 hold both, so four separate 1-GPU jobs sat in QOSMaxJobsPerUserLimit and were cancelled, they were mine). The GPUs are NVIDIA A100-PCIE-40GB (apollo02), driver 580.105, CUDA 13.0 driver, cupy 13.6.0 (cupy-cuda12x in ~/.local for python 3.12, runtime 12.9). Every cupy path that compiles a kernel failed: NVRTC cannot find cuda_fp16.h, because /usr/local/cuda-11.x holds only man pages and cupy found no matching nvidia-cuda-runtime-cu12 headers for python 3.12. The test job ended within a minute with no results.

7. Fix for the headers: cudahdr/include on Typhon is a symlink to the CUDA 12.6 runtime headers of the nvidia-cuda-runtime-cu12 wheel already present in ~/.local/lib/python3.13 (installed earlier for another project, nothing was installed here), and the jobs export CUDA_PATH to cudahdr. cupy then reports NVRTC 12.9 and compiles. Smoke test (cupy_smoke.py, job 17607463 on apollo02) passes for concatenate, QR, eigh and SVD in complex128. On apollo01 the python process hung at the cupy import for more than 5 minutes in two attempts (17607423, 17607434, both cancelled, mine), so every job is pinned to apollo02, which has 7 free GPUs.

8. QOS: GPU jobs of this user fall in gpu-short (2 jobs per user, full with gpu_lean1 and gpu_lean2), gpu-medium (2 per user, ddL_g holds one) or gpu-long (1, held by ddH_g). So this task can hold one gpu-medium job at a time, with up to 8 GPUs in it. The tests run as one 4-GPU job (17607465, apollo02): G2 on one GPU, G3 P0 and P1 on two, the benchmark on the fourth. The benchmark first crashed (random 256-site test state overflowed, the pushed QR factor grew at each step; fixed by rescaling, bench.py only) and was relaunched as an overlap step inside 17607465 (typhon/bench_step.sh, GPU chosen by UUID).

9. pipeline190g.py written for production (one run per GPU, priority order of the spec, configurations from pipeline190.cand_task). Field of sym40_n316 read from the CPU tune file at full precision (−19.976047250613814, RESULT.md prints the rounded −19.976047251).

10. Test job 17607465 (4 GPUs on apollo02, 13 min) plus two overlap steps in it (benchmark, and an L = 32 check of the production path, host below chi 64 then GPU with device SVD). G2, G3 and G4 pass, numbers in RESULT.md. The L = 32 GPU run differs from the CPU run by 1.4e-7 in E after the first mixer sweep at chi 32 (a different equally converged state at residual 1e-7) and by 1e-11 after the chi 96 sweeps. Production settings chosen from G4: GPU_MIN_CHI = 128 (the host is faster below), SVD on the device (2.6 times faster than LAPACK at 768, same accuracy). Measured speed-up per bond at chi 768: 48 with the host SVD, 60 with the device SVD.

11. RESULT.md written with the gates and the registration amendment, before submitting production.

12. Production submitted 21:48: job 17607561 (c190g_prod, apollo02, 6 GPUs, gpu-medium, 2 days, 28 cores, 320 GB) and the resume job 17607562 (same script, --dependency=afterany:17607561). Started at once on GPUs 0 to 5 with p153c L = 256 P0 and P1, and the four p153c L = 128 runs. The pipeline also picked up the CPU tune files that exist (rat_n316 at h_z = −19.975988986408243, sym40_n18 at −19.976333693659278), queued after (a) to (d).

13. The first sweeps (chi 32, host stage, numpy path) differ from the CPU production runs by 1e-5 to 5e-4 in E (for example p153c_L256_P1 sweep 0: −2570.2286139 against −2570.2291280). To see whether this is the ordinary sensitivity of a random-start growth stage to rounding, the original run190.py was run on the first four sweeps of p153c_L128_P0_rand1 with 1, 4 and 6 BLAS threads (chaos/, overlap step in 17607561).

14. The thread test settles item 13: on apollo02 the original CPU code with 4 or 6 threads reproduces the GPU pipeline's host stage digit for digit, and differs from the CPU production run by 1.8e-5 at sweep 0, so the early-stage spread belongs to the CPU code itself. From chi 128 on, GPU and CPU runs agree to 3e-8 to 1.2e-6 in E in the mixer stages and 1e-9 to 2e-8 without mixer (table in RESULT.md). All residual targets met.

15. Production section of RESULT.md written (job ids, order, layout, tested collection commands, first-sweep table, sweep times, ETA). Stopped here as the spec says, with production running.
