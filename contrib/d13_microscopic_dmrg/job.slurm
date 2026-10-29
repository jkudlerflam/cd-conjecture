#!/bin/bash
#SBATCH -J d13dmrg
#SBATCH -p typhon
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 4-00:00:00
#SBATCH --array=1-9
#SBATCH -o logs/%x_%a.out
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
cd "$SLURM_SUBMIT_DIR"
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" tasks.txt)
python3 -B run_mobile.py $line
