#!/bin/bash
#SBATCH --job-name=c157
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=120G
#SBATCH --time=40:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c157/slurm_%j.out
cd /scratch/lustre/jkudlerflam/c157
xargs -P 24 -I{} ./run.sh {} < tasks.txt
echo done
