#!/bin/bash
#SBATCH --job-name=cd-D3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%A_%a.out
set -euo pipefail
# Submit from this directory after copying it to Typhon.
# Choose an allowed partition/account at submission, according to local policy.
# This CPU sparse-ED implementation does not use a GPU.
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
python3 cluster_job.py
