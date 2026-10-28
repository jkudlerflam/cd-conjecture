#!/bin/bash
#SBATCH --job-name=tf_ring
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=03:00:00
#SBATCH --output=logs/ring_%A_%a.out
# Ring ED, one (spec, L) per array task from the task file $1.
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
read NAME SPEC L <<< "$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $1)"
python3 ring_fd.py "$SPEC" $L --nev 6 --out data/ring_${NAME}_L${L}.json
