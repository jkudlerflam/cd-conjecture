#!/bin/bash
#SBATCH --job-name=rA_n2
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=14:00:00
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
# One line of tasks.txt per array index (0-based).
cd /scratch/lustre/jkudlerflam/routeA_n2_chains
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" slurm/$TASKS)
echo "task $SLURM_ARRAY_TASK_ID: $LINE"
python3 -B run_chain.py $LINE
