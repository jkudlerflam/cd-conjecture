#!/bin/bash
#SBATCH --job-name=c165
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=14:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c165/logs/c165_%j.out
# C165 packed job (Claude Opus agent, 2026-11-08).
#  long tasks (22): DMRG L = 128, 256 on the Q2 line (a=0.9, g=0.3, hx=0.2, hz=0.05, Dx=0.03, Dz=1.15..1.3), the
#                   N-stratum floating control and XX; TFI L = 128; Ising x XX ladder L = 64, 128; VUMPS FES to
#                   chi = 256 for Ising x XX, XX, TFI, the floating control and Q2 at Dz = 1.2, 1.25, 1.3.
#  stage A (252):   two-site VUMPS screen chi <= 64 (families A1, A2, A3 of make_tasks.py)
#  stage B:         stageB.py selects up to 24 screen points (criteria in its docstring) -> chi 96..256, DMRG for 4.
cd /scratch/lustre/jkudlerflam/c165
mkdir -p logs fes fes_cal npz dm
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONNOUSERSITE=1
run() { xargs -d '\n' -P "$1" -n 1 bash -c 'eval "$0" >> logs/tasks.log 2>&1'; }
( cat tasks_long.txt | run 22 ) &
cat tasks_A.txt | run 26
python3 stageB.py > logs/stageB.out 2>&1
cat tasks_B.txt | run 26
wait
python3 summarize.py > summary.txt 2>&1
