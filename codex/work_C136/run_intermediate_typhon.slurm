#!/bin/bash
#SBATCH --job-name=C136-dressing
#SBATCH --array=0-5
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=C136-dressing-%A_%a.out

# Refines the narrow-kernel Fredholm quadrature, including failed local cases.
# The high-order dense systems belong on Typhon, not the local workstation.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python3 check_c136.py --mode intermediate --index "${SLURM_ARRAY_TASK_ID}" --max-order 4374 --output "dressing_${SLURM_ARRAY_TASK_ID}.json"
