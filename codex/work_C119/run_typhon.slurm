#!/bin/bash
#SBATCH --job-name=C119-ratio
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=C119-%j.out
#SBATCH --error=C119-%j.err
set -eu
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python3 checks.py --cluster
