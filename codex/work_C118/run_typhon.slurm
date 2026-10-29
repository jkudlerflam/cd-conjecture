#!/bin/bash
#SBATCH --job-name=C118-ratio
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=C118-%j.out
#SBATCH --error=C118-%j.err
set -eu
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python3 checks.py --cluster
