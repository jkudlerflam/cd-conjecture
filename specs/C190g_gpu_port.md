# C190g: a GPU backend for the C190 full-qubit DMRG (2026-11-09)

Coordinator: Claude Fable 5. Jonah approved starting this port on 2026-11-09. Context: REPORT §9.164, `codex/work_C190/` (RESULT.md holds the registration), `specs/C189_C191_decisive_tests.md`.

## Why

C190 is the converged test of the nearest-neighbour qubit chain at the c = 2 candidate. Its decisive reading is c_two from L = 128 to 256, which needs χ = 768 at L = 256. On the CPU path a complex L = 256 run is estimated at 10 to 40 h in total (3 to 6 h per χ = 768 sweep, uncertain by about a factor of 2), and the CPU queue is days deep, so the runs ride inside other allocations.

Measured cost per matvec in z2dmrg.py with 6 threads: 8 ms at χ = 64, 32 to 45 ms at 128, 98 ms at 192 (Typhon), 125 ms at 256 (local). There are 50 to 90 matvecs per bond. The growth sits between χ² and χ³.

Typhon's default python (anaconda3/2024.06, python 3.12, numpy 2.4.6) has cupy 13.6.0. Partition apollo has two nodes (apollo01, apollo02) with 8 GPUs, 256 CPUs and 1 TB each, a 7-day limit, and nothing pending at 21:00 EST on 2026-11-09, with 8 GPUs free. CUDA 11.4 and 11.7 are under /usr/local. Partition h200 has a 1-day limit and a queue.

## Goal

A second, independent and faster path to the registered C190 readings. Same model, same schedules, same convergence targets and the same verdict rules as `codex/work_C190/RESULT.md`. Nothing about the registration changes.

## Hard rules

- Work only in `codex/work_C190/gpu/` (local) and `/scratch/lustre/jkudlerflam/c18x/c190g/` (Typhon). Everything else is read-only. Never edit, move or delete anything in `codex/work_C190/` outside `gpu/`, or anything under `/scratch/lustre/jkudlerflam/c18x/c190/`. The CPU runs are live and import those files, and the C190 agent is still working in `codex/work_C190/`.
- Never cancel, hold, release or signal a job that this task did not submit. Jonah has other GPU jobs on these partitions (gpu_lean1, gpu_lean2, ddL_g, ddH_g, ddS_g). Leave them alone.
- SLURM only. GPU tests go through sbatch or srun on apollo. Never compute on the login node (an import or version check is fine).
- Job names start with c190g. One GPU per run, at most 6 GPUs at once, about 4 CPU cores and at most 64 GB per GPU. Use apollo. Use h200 only if apollo has no free GPU.
- Submitting to Typhon is always OK and needs no approval.
- No pip or conda install. If cupy 13.6.0 does not work on the GPU nodes, try the other anaconda modules (2025.06, 2020.07) and the cuda modules. If nothing works, stop and report.
- complex128 throughout (float64 where the CPU code already uses a real gauge). No single precision and no change to any tolerance.
- Local machine: no GPU, numpy backend only, at most about 8 cores.
- Do not commit or push. No probabilities in any file. No absolute local paths in files you write. Prose with commas, not em-dashes, few colons, no italics.
- No Fable subagents. Astra (gpt_solve.sh) is out of credits, do not call it.
- ssh: `ssh -o LogLevel=ERROR -o ConnectTimeout=20 typhon`, retried with waits of about 30 s, because the bastion rate-limits. A local hook blocks `sleep N; command` chains in Bash, so put waits inside the ssh command string. Another hook blocks `grep -r`, `find .` and `ls -R`, so use explicit paths. Quote remote globs in scp (the local shell is zsh).

## Steps

1. Read `z2dmrg.py`, `run190.py`, `pipeline190.py`, `models.py`, `validate.py`, `analyze190.py`, RESULT.md (registration and R0 to R3) and PROOF.md section 5 (the preconditioner).

2. Profile first, locally with numpy. On the L = 64 p153c problem at χ = 192 to 256, split the time of a mid-chain bond into matvec, preconditioner, Davidson bookkeeping (orthogonalization, Rayleigh–Ritz), SVD and environment updates. Record the table. The χ² to χ³ growth above says part of the time is outside the χ³ contractions, and the port should move that part too.

3. Backend. Copy `z2dmrg.py` to `gpu/z2dmrg_xp.py` with an array-module switch (numpy or cupy).
   - Everything in the local solve lives on the device: the two-site tensor, both environments, the MPO tensors, the preconditioner data and the Davidson basis. No host–device transfers inside the Davidson loop. The small dense eigenproblem of the Davidson subspace may stay on the host.
   - SVD: compare cupy's with numpy's on a production-size block for accuracy and time. Fall back to the host SVD if either is worse.
   - Small blocks (early stages) may be faster on the host. Measure the crossover and switch backend by stage if that helps.
   - Keep checkpoints backend-neutral (numpy arrays), so a run can stop on one backend and resume on the other.
   - Keep the output JSON identical in format to the CPU code, so `analyze190.py` reads it unchanged.

4. Gates. Record each with numbers in `gpu/RESULT.md` before production starts.
   - G1. The numpy backend of `z2dmrg_xp.py` reproduces `z2dmrg.py` at L = 10 (validate.py, all bonds, both sectors) to 1e-10 in energy and 1e-8 in entropy, and one L = 32 run sweep by sweep to 1e-10.
   - G2. The cupy backend passes the same L = 10 validation against ED on an apollo GPU, to the same tolerances.
   - G3. The cupy backend reproduces the L = 64 pilot at the collaborator's field with the pilot's schedule. Reference values at χ = 256, last sweep, from `codex/work_C190/local/out/p153c_L64_P0.json` and `_P1.json`: even sector E = −650.0414852455, S̄₁ = 1.6140, S̄₆ = 1.6085, odd sector E = −650.0414362271, S̄₁ = 1.5937, S̄₆ = 1.5881. Tolerances 2e-8 in E and 5e-4 in S̄.
   - G4. Benchmark on an L = 256 complex problem: time per matvec and per bond at χ = 256, 384, 512 and 768 on the GPU, against the CPU code on the same problem where that is affordable and against the numbers above otherwise. A few mid-chain bonds suffice, and a partially converged state is fine for timing. Report the measured speed-up. If it is below 3 at χ = 768, stop and report, and do not launch production.
   - Note the GPU model (nvidia-smi) and the device memory in use at χ = 768, L = 256.

5. Registration amendment, written in `gpu/RESULT.md` before any GPU production output. The GPU runs are additional runs of the registered configurations. For each key, size and parity sector, the run with the lowest final energy among those passing C1 to C4 is used. Where a CPU and a GPU run agree in energy within 1e-7, their S̄₆ must agree within 0.01 (the existing rule C5), and otherwise the higher one is reported as metastable. Rules V1 to V6 are unchanged.

6. Production, in this order of priority, with the schedules of `pipeline190.SCHED` and the same seeds and initializations, starting from scratch (do not read the CPU checkpoints).
   - (a) p153c (the collaborator's field h_z = −19.97589125) at L = 256, both parity sectors, rand:1.
   - (b) p153c at L = 128, both sectors (rand:1 and the product start), and at L = 64, both sectors.
   - (c) sym40_n316 (h_z = −19.976047251) at L = 64, 128 and 256, both sectors.
   - (d) The XX and two-leg XX controls at the same sizes, only where the CPU controls (node B) have not finished them. They are needed for calibration.
   - (e) The other tuned keys (p153_n316, rat_n316, sym40_n18, sym80_n316) once their fields exist in the CPU tune files under `/scratch/lustre/jkudlerflam/c18x/c190/outA/tune` and `outB/tune` (read only).
   - Ask for the time the benchmark implies plus a margin, at most 2 days, and make sure a job that ends early resumes from its checkpoint.

7. Check the first sweeps of every production run. Residual targets must be met, and energies must be consistent with the CPU run of the same configuration at the same stage (`/scratch/lustre/jkudlerflam/c18x/c190/outA/logs/<run>.log`, read only). Record job ids, the output layout and tested collection commands in `gpu/RESULT.md`. Give an estimate of when the L = 256 runs reach their final stage, with its basis.

8. Stop and report once production is running and checked. Do not wait for the runs to finish. The coordinator's check-in collects and analyzes with `analyze190.py`. If that script needs anything to read the GPU output folder, say exactly what.

## Report

What passed and what failed at each gate, with numbers. The measured speed-up table. Job ids and where the output lands. The time estimate. What is left.
