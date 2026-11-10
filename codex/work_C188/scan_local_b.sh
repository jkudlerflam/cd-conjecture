#!/bin/bash
# C188: local L = 12, 18 scan of the ferromagnetic half of the fundamental triangle, (1, b, -t b)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
out=scan_L12_18_b.jsonl
for b in 1.0 0.9 0.8 0.7 0.6; do
 for t in 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0; do
  g=$(python3 -c "print(round(-$t*$b,6))")
  python3 ringk_run.py 12 1,$b,$g 0,1,2,3,4,5,6 150 $out > /dev/null
  python3 ringk_run.py 18 1,$b,$g 0,1,2,3,4,5,6,7,8,9 200 $out > /dev/null
 done
done
echo done
