#!/bin/bash
# C189 infinite-chain runs for the desktop (see README.txt). Two processes, 4 threads each (about 8 cores).
# Usage: bash run_desktop.sh      (resumable: rerunning skips finished chi values and restarts from out/ckpt/)
cd "$(dirname "$0")"
mkdir -p out logs
python3 validate_inf.py > logs/validate_inf.out 2>&1
echo "validation written to logs/validate_inf.out (energy differences should be below 1e-9)"
cat tasks.txt | xargs -P 2 -I CMD sh -c 'CMD'
python3 summarize_inf.py out | tee summary.txt
