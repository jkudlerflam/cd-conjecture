#!/bin/bash
#SBATCH --job-name=C114-ring
#SBATCH --partition=apollo
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=C114-%j.out
#SBATCH --error=C114-%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python3 checks.py --full
