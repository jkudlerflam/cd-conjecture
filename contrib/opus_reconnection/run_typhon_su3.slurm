#!/bin/bash
#SBATCH -J opus_su3
#SBATCH -p typhon
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=110G
#SBATCH -t 20:00:00
#SBATCH -o opus_su3_%j.out
# Full-spectrum ED, mode typhon_su3 (qubit: L=8 vs 16; su3: L=6 vs 12).
module load anaconda3/2024.06 2>/dev/null
export OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=32 MKL_NUM_THREADS=32
cd "$SLURM_SUBMIT_DIR"
python3 ed_run.py typhon_su3
