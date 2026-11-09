#!/bin/bash
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2
U=1.035276180410083

python3 pairbind2.py $U locked 0.05,0.1,0.2,0.27,0.35,0.5 48 3,4
python3 pairbind2.py 2.0 locked 0.25,0.5 48 3,4
