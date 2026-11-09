#!/bin/bash
#SBATCH --job-name=c184
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=14:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c18x/c184/slurm_%j.out
cd /scratch/lustre/jkudlerflam/c18x/c184
mkdir -p out
# longest tasks first (tasks.txt is sorted by L descending)
cat tasks.txt | xargs -P 48 -L 1 bash run_task.sh
python3 summ.py out > out/summary.txt
