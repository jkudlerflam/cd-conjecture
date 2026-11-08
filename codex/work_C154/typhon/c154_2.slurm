#!/bin/bash
#SBATCH --job-name=c154_2
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c154/slurm_%j.out
cd /scratch/lustre/jkudlerflam/c154
xargs -P 12 -I{} ./run1.sh {} < tasks2.txt
echo done
