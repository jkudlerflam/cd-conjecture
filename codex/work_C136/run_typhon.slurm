#!/bin/bash
#SBATCH --job-name=C136-all-fillings
#SBATCH --array=0-5
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=C136-%A_%a.out

# Submit from this directory with Python, NumPy and SciPy available.
# Supply a site-appropriate CPU partition/account if required. No GPU is used.
# Every L=2,...,256 and every n=1,...,floor(L/2), six anisotropies.
# Paired Bethe rings reach 512. ED cross-checks every available ring <=12.
# Failed solves remain in the output and produce a nonzero exit status.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python3 check_c136.py --mode full --index "${SLURM_ARRAY_TASK_ID}" --output "full_${SLURM_ARRAY_TASK_ID}.jsonl"
