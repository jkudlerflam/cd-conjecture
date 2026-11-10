#!/bin/bash
# C188: corner point beta = 0.2, t = 0.8 at L = 48, 96 after the other local runs
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
while ! grep -q done192b dmrg_local192b.log; do sleep 15; done
python3 dmrg_tau.py 48 16,32,64,128,160,160 dmrg_local_L48.jsonl 1,0.2,-0.16
python3 dmrg_tau.py 96 16,32,64,128,160,160 dmrg_local_L96.jsonl 1,0.2,-0.16
echo doneextra
