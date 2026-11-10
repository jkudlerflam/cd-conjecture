#!/bin/bash
# C188: fine t scans across the boundary, L = 12, 18, branch sectors only (S, N, F).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
out=fine_cross.jsonl
run() { b=$1; t0=$2; t1=$3
  for t in $(python3 -c "import numpy as np; print(' '.join(f'{x:.3f}' for x in np.arange($t0,$t1+1e-9,0.005)))"); do
    g=$(python3 -c "print(round(-$t*$b,6))")
    python3 ringk_run.py 12 1,$b,$g 0,2,4 150 $out > /dev/null
    python3 ringk_run.py 18 1,$b,$g 9,3,6 200 $out > /dev/null
  done; }
run 1.0 0.50 0.66
run 0.7 0.60 0.76
run 0.4 0.70 0.86
echo done
