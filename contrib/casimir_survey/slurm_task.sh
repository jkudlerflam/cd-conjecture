#!/bin/bash
#SBATCH --job-name=casimir
#SBATCH --partition=typhon
#SBATCH --output=logs/%x_%A_%a.out
# usage: sbatch --array=0-N --cpus-per-task=C --mem=XXG --time=HH:MM:SS slurm_task.sh tasks.txt
# each line of tasks.txt is a full python3 command line run inside casimir_survey/
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK PYTHONNOUSERSITE=1
cd /scratch/lustre/jkudlerflam/casimir_survey/casimir_survey
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$1")
echo "task $SLURM_ARRAY_TASK_ID on $(hostname): $LINE"
/usr/bin/time -v $LINE 2>&1 | grep -v "^\s*\(Page\|Socket\|Signals\|Swaps\|File\|Average\|Exit\|Minor\|Major\|Volun\|Invol\)"
