#!/bin/bash
#SBATCH -J opus_qubit
#SBATCH -p typhon
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=40G
#SBATCH -t 08:00:00
#SBATCH -o opus_qubit_%j.out
# Full-spectrum ED, mode typhon_qubit (qubit: L=8 vs 16; su3: L=6 vs 12).
module load anaconda3/2024.06 2>/dev/null
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
cd "$SLURM_SUBMIT_DIR"
python3 ed_run.py typhon_qubit
