#!/bin/bash
#SBATCH -J d2mc
#SBATCH -N 1
#SBATCH -o /scratch/lustre/jkudlerflam/d2_multicritical/logs/%x_%A_%a.out
# usage: sbatch -p PART -c NCPU --mem=M -t T --array=0-(n-1) batch.slurm TASKFILE
# every line of TASKFILE is "OUTNAME ARGS-for-scan.py"; array task i runs chunk i/n of each line
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /scratch/lustre/jkudlerflam/d2_multicritical
C=$SLURM_ARRAY_TASK_ID/$SLURM_ARRAY_TASK_COUNT
while read -r NAME ARGS; do
  [ -z "$NAME" ] && continue
  python3 scan.py $ARGS --procs $SLURM_CPUS_PER_TASK --chunk $C --out data/maps/${NAME}.c${SLURM_ARRAY_TASK_ID}.jsonl
done < $1
