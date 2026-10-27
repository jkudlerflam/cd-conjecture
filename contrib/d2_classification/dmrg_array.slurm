#!/bin/bash
#SBATCH --job-name=d2dmrg
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=6:00:00
#SBATCH --output=logs/dmrg_%A_%a.out
# Open-chain two-site DMRG, one (chain, L, seed) per task.  chi = 32, 64, 128, 192, two sweeps each.
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
read NAME L SEED <<< "$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" tasks_dmrg.txt)"
python3 dmrg.py "$NAME" --L "$L" --seed "$SEED" --chi 32 64 128 192 --sweeps 2 --out "results/dmrg_${NAME}_L${L}_s${SEED}.json"
