#!/bin/bash
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2
U=1.035276180410083
python3 pairbind2.py $U locked 0.3,0.4 48 4 >> pb2_raw.out 2>> pb2.err
python3 pair_conv.py pb2_raw.out 0.3 6 > pairconv_v03.out 2> pairconv_v03.err
python3 pair_conv.py pb2_raw.out 0.4 6 > pairconv_v04.out 2> pairconv_v04.err
python3 pair_conv.py pb2_raw.out 0.1 6 > pairconv_v01.out 2> pairconv_v01.err
python3 pair_conv.py pb2_raw.out 0.05 6 > pairconv_v005.out 2> pairconv_v005.err
