#!/bin/bash
# C188: local L = 24 rings on the beta = 0.7 and 0.4 cuts near the crossing (sectors 0, 1, 4, 8, 12)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=2
while ! grep -q done local24.log; do sleep 20; done
for bt in "0.7 0.64" "0.7 0.655" "0.7 0.665" "0.7 0.675" "0.7 0.69" "0.4 0.76" "0.4 0.775" "0.4 0.785" "0.4 0.795" "0.4 0.81" "0.7 0.60" "0.4 0.70"; do
  b=${bt% *}; t=${bt#* }
  g=$(python3 -c "print(round(-$t*$b,6))")
  python3 ringk_run.py 24 1,$b,$g 0,1,4,8,12 200 local24.jsonl > /dev/null
done
echo done2
