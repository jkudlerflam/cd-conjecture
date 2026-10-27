#!/bin/bash
#SBATCH --job-name=d2fu_ed
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=16
#SBATCH --mem=370G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/%A_%a.out
# Momentum-resolved ring ED of points 220 and 1282, six levels per sector.
# L = 8..26 with ground-state entanglement, then L = 27, 28 without it (memory).
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
NAMES=(220 1282)
NAME=${NAMES[$SLURM_ARRAY_TASK_ID]}
python3 run_point.py "$NAME" --k 6 --procs 16 --out results/${NAME}_k6.json --sizes 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
python3 run_point.py "$NAME" --k 6 --procs 16 --no-ent --out results/${NAME}_k6.json --sizes 27
python3 run_point.py "$NAME" --k 6 --procs 12 --no-ent --out results/${NAME}_k6.json --sizes 28
