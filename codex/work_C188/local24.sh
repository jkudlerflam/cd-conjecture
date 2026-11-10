#!/bin/bash
# C188: local L = 24 rings on the beta = 1 cut, branch sectors S (0), N (4), F (8) and current (1), plus Heisenberg control
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=2
python3 ringk_run.py 24 1,1,1 0,1,4,8,12 200 local24.jsonl > /dev/null
for t in 0.50 0.55 0.53 0.56 0.57 0.58 0.60 0.45; do
  g=$(python3 -c "print(round(-$t,6))")
  python3 ringk_run.py 24 1,1.0,$g 0,1,4,8,12 200 local24.jsonl > /dev/null
done
echo done
