#!/bin/bash
# twist-free valley spacing: |w| chosen so that (k_b - k_a) L = 0 mod 2 pi with the bare valleys (s = pi j / L)
cd "$(dirname "$0")/.."
export NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
while pgrep -f "spec2.py 24 6 1.5707963267948966 1.0 cusp" > /dev/null || pgrep -f "local/v2_cusp.jsonl" > /dev/null; do sleep 30; done
W6=$(python3 -c "import math;print(4*math.sin(math.pi/12))")   # n = 1/6, L = 24, 36, 48
W5=$(python3 -c "import math;print(4*math.sin(math.pi/10))")   # n = 1/5, L = 20, 30, 40, 50
W4=$(python3 -c "import math;print(4*math.sin(math.pi/8))")    # n = 1/4, L = 16, 24, 32, 40
for L in 24 36; do python3 spec2.py $L $((L/6)) 1.5707963267948966 $W6 cusp 8 8 half local/v2_tf.jsonl >> local/v2_tf.out 2>&1; done
for L in 20 30; do python3 spec2.py $L $((L/5)) 1.5707963267948966 $W5 cusp 8 8 half local/v2_tf.jsonl >> local/v2_tf.out 2>&1; done
for L in 16 24 32; do python3 spec2.py $L $((L/4)) 1.5707963267948966 $W4 cusp 8 8 half local/v2_tf.jsonl >> local/v2_tf.out 2>&1; done
echo chain3 done >> local/v2_tf.out
