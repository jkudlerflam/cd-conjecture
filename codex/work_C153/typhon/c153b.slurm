#!/bin/bash
#SBATCH --job-name=c153b
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --time=20:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c153/slurm_%j.out
cd /scratch/lustre/jkudlerflam/c153
mkdir -p logs fes_data strip_data grid_data
xargs -P 48 -I{} ./run_task.sh {} < tasks2.txt
echo done
