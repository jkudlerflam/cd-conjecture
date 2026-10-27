#!/bin/bash
#SBATCH --job-name=cdidmrg
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/%x_%A_%a.out
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
cd /scratch/lustre/jkudlerflam/cdsearch
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$1")
echo "task $SLURM_ARRAY_TASK_ID on $(hostname): $LINE"
python3 -u idmrg.py $LINE
