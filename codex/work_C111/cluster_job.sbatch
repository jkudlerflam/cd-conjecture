#!/bin/bash
#SBATCH --job-name=C111-XX
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=C111-XX-%j.out

# Submit from a staged copy of work_C111 on Typhon, choosing a valid partition.
# No large computation is permitted on the login node.
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python3 large_gaussian.py
