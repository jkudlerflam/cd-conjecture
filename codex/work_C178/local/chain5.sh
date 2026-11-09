#!/bin/bash
cd "$(dirname "$0")/.."
export NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
P=1.5707963267948966
run(){ python3 spec2.py $1 $2 $P $3 cusp $4 8 half local/$5.jsonl >> local/$5.out 2>&1; }
for W in 1.0 2.0; do for LN in 24:2 36:3 48:4 48:3 48:5 60:5; do run ${LN%:*} ${LN#*:} $W 8 dil; done; done
W6=$(python3 -c "import math;print(4*math.sin(math.pi/12))"); W5=$(python3 -c "import math;print(4*math.sin(math.pi/10))"); W4=$(python3 -c "import math;print(4*math.sin(math.pi/8))")
for LN in 20:4 30:6; do run ${LN%:*} ${LN#*:} $W5 8 tf; done
for LN in 24:4 16:4 24:6; do run ${LN%:*} ${LN#*:} $W6 8 tf; done
for LN in 16:4 24:6; do run ${LN%:*} ${LN#*:} $W4 8 tf; done
for LN in 24:3 24:5 36:5 36:7 24:2 24:6; do run ${LN%:*} ${LN#*:} 1.0 2 dN; done
for LN in 24:4 24:3 24:5 36:6 36:5 36:7; do run ${LN%:*} ${LN#*:} 2.0 4 dN; done
run 36 6 $W6 8 tf
for LN in 36:4 36:8; do run ${LN%:*} ${LN#*:} 1.0 2 dN; done
echo chain5 done >> local/dil.out
