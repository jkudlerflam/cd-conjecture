#!/bin/bash
#SBATCH --job-name=C105-XX
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=C105-XX-%j.out

# Submit from the staged work_C105 directory on Typhon.
# The partition must be selected from the site's available CPU partitions.
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python3 checks.py --cluster --output large_gaussian_results.json
