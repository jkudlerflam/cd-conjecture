#!/bin/bash
#SBATCH --job-name=encobc
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/obc_%A_%a.out
# Open-chain DMRG, one line of obc_tasks.txt per array task (submit with --array=0-N).
cd /scratch/lustre/jkudlerflam/encoded_chain_numerics/encoded_chain_numerics
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" obc_tasks.txt)
echo "task $SLURM_ARRAY_TASK_ID on $(hostname): $LINE"
python3 -u fdmrg_md.py $LINE
