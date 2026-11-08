#!/bin/bash
#SBATCH --job-name=c152_j2
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=14:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c152/j2/slurm_%j.out
cd /scratch/lustre/jkudlerflam/c152/j2
xargs -P 48 -L 1 ./run_d2.sh < dtasks2.txt
echo done
