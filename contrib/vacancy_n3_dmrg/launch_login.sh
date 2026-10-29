#!/bin/bash
# Run every line of a task file in the background on the Typhon login node at nice 19, T threads each.
# usage: bash launch_login.sh TASKFILE [T]
T=${2:-3}
export OMP_NUM_THREADS=$T OPENBLAS_NUM_THREADS=$T MKL_NUM_THREADS=$T
mkdir -p data/logs state
while IFS= read -r line; do
  [ -z "$line" ] && continue
  name=$(echo "$line" | awk '{print $4}')
  nohup nice -n 19 bash -c "$line" > data/logs/$name.log 2>&1 &
done < "$1"
