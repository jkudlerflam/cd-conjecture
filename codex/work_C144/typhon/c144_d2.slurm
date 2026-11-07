#!/bin/bash
#SBATCH --job-name=c144_d22
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c144/slurm_%j.out
cd /scratch/lustre/jkudlerflam/c144
xargs -P 24 -I{} ./run_d2.sh {} < dtasks2.txt
echo done
