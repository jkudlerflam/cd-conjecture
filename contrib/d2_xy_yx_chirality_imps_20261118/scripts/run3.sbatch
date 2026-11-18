#!/bin/bash
#SBATCH -p typhon
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH --time-min=01:30:00
#SBATCH --signal=B:USR1@420
#SBATCH -o logs3/%x_%j.out
# self-resubmitting, resumable:  sbatch -J name run3.sbatch key=val ...   (stages resume from out3/*.h5)
cd ~/cc_spinchain/v3
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export JOB_END=$(squeue -j $SLURM_JOB_ID -h -o %e); export JOB_END_EPOCH=$(TZ=America/New_York date -d "$JOB_END" +%s)
echo "# $(hostname) $(date) end=$JOB_END epoch_left_h=$(( (JOB_END_EPOCH-$(date +%s))/3600 )) : $@"
resub() { echo "# time limit near: resubmitting $(date)"; kill $PID 2>/dev/null; sleep 5; sbatch -J $SLURM_JOB_NAME run3.sbatch "$@"; exit 0; }
trap 'resub "$@"' USR1
PYTHONPATH= ~/cc_spinchain/venv/bin/python -s -u ims.py "$@" &
PID=$!
wait $PID
rc=$?
if [ $rc -eq 3 ]; then echo "# stage incomplete: resubmitting"; sbatch -J $SLURM_JOB_NAME run3.sbatch "$@"; fi
if [ $rc -ne 0 ] && [ $rc -ne 3 ] && [ $rc -lt 128 ]; then echo "# python failed rc=$rc, not resubmitting"; fi
echo "# finished $(date) rc=$rc"
