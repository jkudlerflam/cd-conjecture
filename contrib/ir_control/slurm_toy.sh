#!/bin/bash
#SBATCH --job-name=irtoy
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
# usage: sbatch --array=0-14 slurm_toy.sh tasks.txt   (run from the folder holding toy_ir.py and idmrg.py)
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
cd /scratch/lustre/jkudlerflam/ir_control
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$1")
echo "task $SLURM_ARRAY_TASK_ID on $(hostname): $LINE"
python3 -u toy_ir.py $LINE --out results/task_${SLURM_ARRAY_TASK_ID}.jsonl
