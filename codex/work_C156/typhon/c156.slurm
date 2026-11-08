#!/bin/bash
#SBATCH --job-name=c156_k110
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c156/slurm_%j.out
cd /scratch/lustre/jkudlerflam/c156
xargs -P 20 -I{} ./run.sh {} < tasks.txt
echo done
