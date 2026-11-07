#!/bin/bash
#SBATCH --job-name=C135-deeper
#SBATCH --array=0-61
#SBATCH --cpus-per-task=4
#SBATCH --mem=192G
#SBATCH --time=3-00:00:00
#SBATCH --output=slurm-deeper-%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OPENBLAS_NUM_THREADS=4
export OMP_NUM_THREADS=4
# Optional expensive extension for the n=7,8 two-level differences.
# Generic complex n=10 is not benchmarked and may remain computationally infeasible.
python3 campaign.py --large --index "${SLURM_ARRAY_TASK_ID}" --nmax 10
