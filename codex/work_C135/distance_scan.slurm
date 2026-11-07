#!/bin/bash
#SBATCH --job-name=C135-distance
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm-distance-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OPENBLAS_NUM_THREADS=4
export OMP_NUM_THREADS=4
# Run after the ring dictionary jobs have finished. Includes n=9 for L=8.
python3 distances.py --large
python3 report.py
