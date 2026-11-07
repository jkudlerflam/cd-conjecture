#!/bin/bash
#SBATCH --job-name=C134-half
#SBATCH --array=0-7
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=24:00:00
#SBATCH --output=C134-half-%A_%a.out

# Submit from work_C134 in a Python environment with NumPy and SciPy.
# Supply a site-appropriate CPU partition/account to sbatch if required.
# No GPU is used. Failed solves are recorded and make the job fail.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
python3 half_scan.py --mode large --index "${SLURM_ARRAY_TASK_ID}" --output "half_large_${SLURM_ARRAY_TASK_ID}.jsonl"
