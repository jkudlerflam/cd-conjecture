#!/bin/bash
# usage: runone.sh <tag> <threads>
# Continue one run from its checkpoint. A lock with a heartbeat keeps two copies of a tag from running at once.
tag=$1; thr=$2
cd /scratch/lustre/jkudlerflam/c187d || exit 1
mkdir -p locks logs
lock=locks/$tag
if ! mkdir "$lock" 2>/dev/null; then
  age=$(( $(date +%s) - $(stat -c %Y "$lock/hb" 2>/dev/null || echo 0) ))
  if [ "$age" -lt 600 ]; then echo "$tag is running elsewhere, skipped"; exit 0; fi
  echo "$tag has a stale lock, taking over"
fi
touch "$lock/hb"
( while true; do sleep 60; touch "$lock/hb"; done ) &
hb=$!
trap 'kill $hb 2>/dev/null; rm -rf "$lock"' EXIT
export OMP_NUM_THREADS=$thr OPENBLAS_NUM_THREADS=$thr MKL_NUM_THREADS=$thr
python3 run.py "cfg/$tag.json" out >> "logs/$tag.log" 2>&1
echo "$tag finished with exit code $? on $(hostname) at $(date)"
