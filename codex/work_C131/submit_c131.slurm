#!/bin/bash
#SBATCH --job-name=C131-cover
#SBATCH --partition=apollo
#SBATCH --array=0-223%8
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=24:00:00
#SBATCH --output=c131-%A_%a.out
#SBATCH --error=c131-%A_%a.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
# Requires python3 with numpy and scipy in the submitting environment.
# This is a CPU sparse-ED calculation; no GPU code or GPU allocation is needed.
python3 optimize_q.py --task "${SLURM_ARRAY_TASK_ID:?}" --budget 4000
