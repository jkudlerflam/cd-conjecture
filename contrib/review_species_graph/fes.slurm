#!/bin/bash
#SBATCH -J rsg_fes
#SBATCH -p typhon
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 23:00:00
#SBATCH -a 1-8
#SBATCH -o slurm_%A_%a.out
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
cd /scratch/lustre/jkudlerflam/review_species_graph
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" tasks.txt)
name=$(echo $line | cut -d' ' -f1); args=$(echo $line | cut -d' ' -f2-)
python3 fes_effective.py $name $args --out data/$name.json
