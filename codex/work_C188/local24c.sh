#!/bin/bash
# C188: corner check at L = 24 (beta = 0.1, 0.2), after the local DMRG finishes
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=2
while pgrep -f "dmrg_tau.py 48" > /dev/null; do sleep 30; done
for bt in "0.2 0.3" "0.2 0.6" "0.2 0.8" "0.1 0.3" "0.1 0.6" "0.2 0.85" "0.1 0.8"; do
  b=${bt% *}; t=${bt#* }
  g=$(python3 -c "print(round(-$t*$b,6))")
  python3 ringk_run.py 24 1,$b,$g 0,1,4,8,12 250 local24.jsonl > /dev/null
done
echo done3
