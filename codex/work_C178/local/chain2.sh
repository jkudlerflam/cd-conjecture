#!/bin/bash
# runs after the old L=42 spec.py run finishes: legs control (twist-resolved), and N +- 2 runs for compressibility
cd "$(dirname "$0")/.."
export NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
while pgrep -f "spec.py 42 7" > /dev/null; do sleep 30; done
for LN in 24:4 30:5 36:6 20:4 30:6 16:4 24:6 32:8 32:4; do L=${LN%:*}; N=${LN#*:}
  python3 spec2.py $L $N 1.5707963267948966 1.0 legs 4 8 half local/v2_legs.jsonl >> local/v2_legs.out 2>&1; done
for LN in 24:2 24:6 30:3 30:7 36:4 36:8 20:2 20:6 30:4 30:8 16:2 16:6 24:4 24:8; do L=${LN%:*}; N=${LN#*:}
  python3 spec2.py $L $N 1.5707963267948966 1.0 cusp 2 8 half local/v2_cusp_dN.jsonl >> local/v2_cusp_dN.out 2>&1; done
echo chain2 done >> local/v2_legs.out
