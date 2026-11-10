#!/bin/bash
#SBATCH --job-name=c188
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=90G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c18x/c188/slurm_%j.out
cd /scratch/lustre/jkudlerflam/c18x/c188
mkdir -p out
( while read -r line; do ./run_task.sh "$line"; done < tasks_b.txt; echo "stream b done $(date)" ) &
xargs -P 16 -I{} ./run_task.sh {} < tasks_a.txt
echo "stream a done $(date)"
wait
echo done
