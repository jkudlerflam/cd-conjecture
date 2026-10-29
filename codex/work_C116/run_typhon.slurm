#!/bin/bash
#SBATCH --job-name=C116-rings
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=C116-%j.out
#SBATCH --error=C116-%j.err
set -eu
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python3 check_rings.py --cluster
