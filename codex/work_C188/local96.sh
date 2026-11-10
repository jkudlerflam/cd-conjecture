#!/bin/bash
# C188: local L = 96 DMRG at the points job 2 had not reached, after the local L = 192 run ends
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
while kill -0 69985 2>/dev/null; do sleep 15; done
python3 dmrg_tau.py 96 16,32,64,128,160,160 dmrg_local_L96.jsonl 1,0.2,-0.12 1,0.2,-0.17 1,0.7,-0.462 1,0.4,-0.316
echo done96
