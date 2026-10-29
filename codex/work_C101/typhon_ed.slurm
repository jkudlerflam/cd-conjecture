#!/bin/bash
#SBATCH --job-name=C101-ed
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=C101-ed-%j.out
set -euo pipefail
: "${SLURM_JOB_ID:?This script must run as a SLURM job}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "${SLURM_SUBMIT_DIR}"
python3 adversarial_ed.py --cluster
