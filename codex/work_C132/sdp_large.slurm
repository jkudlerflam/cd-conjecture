#!/bin/bash
#SBATCH --job-name=C132-sdp
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=C132-sdp-%j.out
export OPENBLAS_NUM_THREADS=4
export OMP_NUM_THREADS=4
python3 scan_sdp.py 1024 --long
