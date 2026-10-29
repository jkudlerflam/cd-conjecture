#!/bin/bash
# usage: sbatch node.slurm TASKFILE  (one node, every task line in parallel; resumes from state/ checkpoints)
#SBATCH --partition=typhon
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=300G
#SBATCH --time=0-03:00:00
#SBATCH --output=data/logs/node_%j.out
export OMP_NUM_THREADS=9 OPENBLAS_NUM_THREADS=9 MKL_NUM_THREADS=9
cd $SLURM_SUBMIT_DIR
while IFS= read -r line; do name=$(echo "$line" | awk "{print \$4}"); bash -c "$line" > data/logs/$name.log 2>&1 & done < "${1:-tasks_obc.txt}"
wait
