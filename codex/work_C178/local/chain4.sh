#!/bin/bash
# N +- 1 at n = 1/6 for the parity (odd-even) gap, twist-relaxed; runs after chain2
cd "$(dirname "$0")/.."
export NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
while pgrep -f chain2.sh > /dev/null; do sleep 30; done
for LN in 24:3 24:5 36:5 36:7 20:3 20:5 30:5 30:7; do L=${LN%:*}; N=${LN#*:}
  python3 spec2.py $L $N 1.5707963267948966 1.0 cusp 2 8 half local/v2_cusp_dN.jsonl >> local/v2_cusp_dN.out 2>&1; done
echo chain4 done >> local/v2_cusp_dN.out
