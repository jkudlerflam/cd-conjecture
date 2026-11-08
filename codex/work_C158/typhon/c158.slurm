#!/bin/bash
#SBATCH --job-name=c158_fes
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --time=23:55:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c158/slurm_%j.out
cd /scratch/lustre/jkudlerflam/c158
mkdir -p fes state dmrg logs
xargs -P 48 -d '\n' -I{} ./run.sh {} < tasks.txt
echo done
