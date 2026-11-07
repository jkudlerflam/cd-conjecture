#!/bin/bash
#SBATCH --job-name=C129-bethe
#SBATCH --array=0-10
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=24:00:00
#SBATCH --output=C129-bethe-%A_%a.out

# Submit from the directory containing these scripts with a NumPy/SciPy environment.
# This is a CPU job. Set an available CPU partition/account at submission if required.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
python3 bethe_scan.py --index "${SLURM_ARRAY_TASK_ID}" --output "bethe_large_${SLURM_ARRAY_TASK_ID}.jsonl"
