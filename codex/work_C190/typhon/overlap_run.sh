#!/bin/bash
# Coordinator addition (2026-11-09). Runs a C190 node payload as an extra step inside another of our own running
# allocations that has idle cores, so that it starts before the queued job does. The body is the same as
# c190a.slurm / c190b.slurm. The step cancels itself as soon as the queued job is running, so two pipelines
# never write the same output directory for more than a few seconds (checkpoints are written atomically).
# usage, from the login node, in the scratch directory:
#   nohup setsid srun --jobid=<donor job> --overlap -n1 -c48 bash overlap_run.sh <A|B> <slots> <queued job id> \
#       > logs/srun_<A|B>.out 2>&1 < /dev/null &
NODE=$1; SLOTS=$2; REAL=$3
cd /scratch/lustre/jkudlerflam/c18x/c190 || exit 1
mkdir -p logs out$NODE
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "$(date) start node $NODE, $SLOTS slots, step $SLURM_JOB_ID.$SLURM_STEP_ID on $(hostname)" >> logs/overlap_$NODE.log
(
  python3 validate.py 10 > logs/validate_${NODE}_ov.log 2>&1
  python3 pipeline190.py --node $NODE --slots $SLOTS --threads 6 --out out$NODE > logs/pipeline_${NODE}_ov.log 2>&1
  python3 analyze190.py out$NODE > out$NODE/summary.txt 2>&1
  echo "$(date) payload finished" >> logs/overlap_$NODE.log
) &
PAY=$!
while kill -0 $PAY 2>/dev/null; do
  st=$(squeue -h -j $REAL -o %T 2>/dev/null)
  if [ "$st" = "RUNNING" ]; then
    echo "$(date) queued job $REAL is running, cancelling this step" >> logs/overlap_$NODE.log
    scancel $SLURM_JOB_ID.$SLURM_STEP_ID
    sleep 120; exit 0
  fi
  sleep 10
done
