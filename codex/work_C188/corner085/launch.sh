#!/bin/bash
# C188 corner follow-up, local run: L = 48 then 96 with 2 threads, L = 192 with 3 threads.
cd "$(dirname "$0")"
( export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2
  python3 run085.py 48 1,0.2,-0.17 > log_L48.txt 2>&1
  python3 run085.py 96 1,0.2,-0.17 > log_L96.txt 2>&1 ) &
echo $! > pid_A.txt
( export OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 MKL_NUM_THREADS=3 VECLIB_MAXIMUM_THREADS=3
  python3 run085.py 192 1,0.2,-0.17 > log_L192.txt 2>&1 ) &
echo $! > pid_B.txt
