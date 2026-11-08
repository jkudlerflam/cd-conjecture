#!/bin/bash
#SBATCH --job-name=c154_1
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c154/slurm_%j.out
cd /scratch/lustre/jkudlerflam/c154
xargs -P 21 -I{} ./run1.sh {} < tasks1.txt
echo done
