#!/bin/bash
# usage: run_task.sh THREADS TAG L N PHI WABS MODEL NLEV NTHETA TODO
TH=$1; TAG=$2; shift 2
cd /scratch/lustre/jkudlerflam/c178
case "$TAG" in *_L72_*|cusp_p1.54_w1.0_L48_N8|cusp_p1.50_w1.0_L48_N8|cusp_ppi2_w0.5_L48_N8|rods_ppi2_w1.0_L48_N8|rods_ppi2_w1.0_L40_N10) echo "$TAG skipped (priority/size) $(date)" >> logs/done.txt; exit 0;; esac
mkdir -p out/locks
mkdir out/locks/$TAG 2>/dev/null || exit 0   # another worker pool already took this task
export NUMBA_NUM_THREADS=$TH OMP_NUM_THREADS=$TH OPENBLAS_NUM_THREADS=$TH MKL_NUM_THREADS=$TH PYTHONNOUSERSITE=1
/opt/sns/anaconda3/2024.06/bin/python3 -s -B spec2.py "$@" out/$TAG.jsonl > out/$TAG.out 2>&1
echo "$TAG exit $? $(date)" >> logs/done.txt
