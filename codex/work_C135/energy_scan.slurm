#!/bin/bash
#SBATCH --job-name=C135-energy
#SBATCH --array=0-61
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm-energy-%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OPENBLAS_NUM_THREADS=4
export OMP_NUM_THREADS=4
# Uses only already installed numpy, scipy, cvxpy and solvers.
# Each task resumes existing results. Generic n=9,10 may exceed this allocation;
# first complete n<=8, then request the optional deeper run separately.
python3 campaign.py --large --index "${SLURM_ARRAY_TASK_ID}" --nmax 8
