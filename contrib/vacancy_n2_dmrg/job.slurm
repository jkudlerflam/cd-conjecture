#!/bin/bash
#SBATCH --partition=typhon
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --output=data/logs/%x_%a.out
# usage: sbatch --job-name=NAME --array=0-N --cpus-per-task=C --mem=XG job.slurm TASKFILE
# first pass: --partition=apollo --qos=short (CPU only), preempted after 26 min; FES: typhon, long
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
cd "${SLURM_SUBMIT_DIR:-.}"
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$1")
echo "task $SLURM_ARRAY_TASK_ID: $line"; date
eval "$line"
date
