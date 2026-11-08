#!/bin/bash
#SBATCH --job-name=c153
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c153/slurm_%j.out
cd /scratch/lustre/jkudlerflam/c153
mkdir -p logs
xargs -P 48 -I{} ./run_task.sh {} < tasks.txt
echo done
