#!/bin/bash
# C188: corner point at L = 192 after local96.sh finishes
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
while ! grep -q done96 dmrg_local96.log; do sleep 15; done
python3 dmrg_tau.py 192 16,32,64,128,160,200,200 dmrg_local_L192.jsonl 1,0.2,-0.12
echo done192b
