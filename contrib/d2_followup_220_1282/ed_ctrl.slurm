#!/bin/bash
#SBATCH --job-name=d2fu_edc
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/edc_%j.out
# Ring ED of two c = 1 incommensurate controls, six levels per sector, L = 8..24.
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
for NAME in xydm xxzm07_dm_h; do
  python3 run_point.py "$NAME" --k 6 --procs 24 --out results/${NAME}_k6.json --sizes 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
done
