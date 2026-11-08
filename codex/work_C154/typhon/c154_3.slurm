#!/bin/bash
#SBATCH --job-name=c154_3
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c154/slurm_%j.out
cd /scratch/lustre/jkudlerflam/c154
xargs -P 46 -I{} ./run3.sh {} < tasks3.txt
echo done
