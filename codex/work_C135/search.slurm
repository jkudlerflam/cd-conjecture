#!/bin/bash
#SBATCH --job-name=C135-search
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm-search-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OPENBLAS_NUM_THREADS=4
export OMP_NUM_THREADS=4
# Full eight-angle searches. Generic n=6,8 SDPs are the expensive part.
python3 optimize.py --large --starts 8 --budget 300
python3 report.py
