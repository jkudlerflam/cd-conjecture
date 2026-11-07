#!/bin/bash
#SBATCH --job-name=C126_ed
#SBATCH --partition=apollo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=C126_%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
python3 checks.py --cluster
python3 algebra_checks.py
python3 copy_rotation_checks.py
