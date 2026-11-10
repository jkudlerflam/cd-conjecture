#!/bin/bash
#SBATCH --job-name=c188
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=90G
#SBATCH --time=11:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c18x/c188/slurm_%j.out
cd /scratch/lustre/jkudlerflam/c18x/c188
mkdir -p out2
( for i in 1 2 3 4; do ./run_task2.sh "ring30 x"; done; echo "stream b done $(date)" ) &
xargs -P 16 -I{} ./run_task2.sh {} < tasks_a2.txt
echo "stream a done $(date)"
wait
echo done
