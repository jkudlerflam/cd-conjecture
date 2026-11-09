#!/bin/bash
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2
for L in 20 24; do
  python3 onset.py $L 1.0 0.1,0.25,0.4,0.6,1.0 locked
  python3 onset.py $L 2.0 0.5,1.0,2.0 locked
  python3 onset.py $L 0.5 0.03,0.0625,0.1,0.2 locked
  python3 onset.py $L 1.0 0.1,0.2,0.3 zero
done
python3 onset.py 28 1.0 0.25,0.4 locked
python3 onset.py 28 1.0 0.0,0.35 zero
