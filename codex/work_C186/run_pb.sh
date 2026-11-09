#!/bin/bash
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2
python3 pairbind.py 2.0 zero 0.0,0.5,1.0,1.5 40,48 2
python3 pairbind.py 2.0 locked 0.25,0.5,1.0,1.5,2.0 40,48 4,5
python3 pairbind.py 1.5 zero 0.0,0.3,0.5625,0.8 40,48 2
python3 pairbind.py 1.5 locked 0.3,0.5625,0.8,1.2 40,48 4,5
