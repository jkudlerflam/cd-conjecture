#!/bin/bash
#SBATCH --job-name=d2dmrgP
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=40
#SBATCH --mem=40G
#SBATCH --time=10:00:00
#SBATCH --output=logs/dmrgP_%A_%a.out
# Packed DMRG: array task i runs lines i, i+P, i+2P, ... of tasks_dmrg.txt (offset 4),
# ten at a time with 4 BLAS threads each.
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
P=${SLURM_ARRAY_TASK_COUNT}
awk -v i=$SLURM_ARRAY_TASK_ID -v P=$P 'NR>4 && ((NR-5)%P)==i' tasks_dmrg.txt | \
  xargs -P 10 -L 1 sh -c 'python3 dmrg.py "$0" --L "$1" --seed "$2" --chi 32 64 128 192 --sweeps 2 --out "results/dmrg_$0_L$1_s$2.json" > "logs/dmrg_$0_L$1_s$2.log" 2>&1'
