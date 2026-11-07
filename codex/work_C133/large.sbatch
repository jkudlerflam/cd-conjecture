#!/bin/bash
#SBATCH --job-name=C133-rings
#SBATCH --partition=apollo
#SBATCH --array=0-4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=C133-%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
families=(xxz ashkin_teller spin_orbital spin1_blbq fz_parafermion)
python3 -B scan.py --large --family "${families[$SLURM_ARRAY_TASK_ID]}"
