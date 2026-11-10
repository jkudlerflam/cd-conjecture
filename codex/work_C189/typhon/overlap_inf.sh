#!/bin/bash
# Coordinator addition (2026-11-09). Runs one C189 infinite-chain job body as an extra step inside another of our
# own running allocations that has idle cores, so that it starts before the queued job does. run_inf.py resumes
# from its checkpoint (written atomically). The step cancels itself as soon as the queued job is running.
# usage, from the login node, in the scratch directory:
#   nohup setsid srun --jobid=<donor job> --overlap -n1 -c48 bash typhon/overlap_inf.sh <tag> <threads> <queued job id> \
#       > logs/srun_<tag>.out 2>&1 < /dev/null &
TAG=$1; TH=$2; REAL=$3
cd /scratch/lustre/jkudlerflam/c18x/c189 || exit 1
mkdir -p logs out
export OMP_NUM_THREADS=$TH OPENBLAS_NUM_THREADS=$TH MKL_NUM_THREADS=$TH
echo "$(date) start $TAG, $TH threads, step $SLURM_JOB_ID.$SLURM_STEP_ID on $(hostname)" >> logs/overlap_$TAG.log
(
  python3 run_inf.py typhon/cfg/$TAG.json out > logs/${TAG}_ov.log 2>&1
  python3 summarize_inf.py out > summary_inf_$TAG.txt 2>&1
  echo "$(date) payload finished" >> logs/overlap_$TAG.log
) &
PAY=$!
while kill -0 $PAY 2>/dev/null; do
  st=$(squeue -h -j $REAL -o %T 2>/dev/null)
  if [ "$st" = "RUNNING" ]; then
    echo "$(date) queued job $REAL is running, cancelling this step" >> logs/overlap_$TAG.log
    scancel $SLURM_JOB_ID.$SLURM_STEP_ID
    sleep 120; exit 0
  fi
  sleep 5
done
