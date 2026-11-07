#!/bin/bash
#SBATCH --job-name=C128-histories
#SBATCH --partition=apollo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=C128-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OPENBLAS_NUM_THREADS=4
export OMP_NUM_THREADS=4
python3 checks.py --cluster
