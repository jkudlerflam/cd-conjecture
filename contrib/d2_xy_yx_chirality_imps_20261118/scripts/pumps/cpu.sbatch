#!/bin/bash
#SBATCH --partition=typhon
#SBATCH --qos=short
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=23:00:00
#SBATCH --output=logs/%x_%j.out
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=32 MKL_NUM_THREADS=32 PUMPS_BACKEND=numpy
bash runp.sh "$@"
