#!/bin/bash
# Prepare points first, then run on the cluster, for example:
# bash submit_cluster.sh --partition=YOUR_AUTHORIZED_PARTITION
# This wrapper is supplied for the user. It was not executed in the local run.
set -euo pipefail
point_file=${D2_POINTS:-results/global_points.json}
total_points=$(python3 -c 'import json,sys; print(len(json.load(open(sys.argv[1]))))' "$point_file")
batch_size=${D2_BATCH_SIZE:-1000}
concurrency=${D2_CONCURRENCY:-48}
for (( point_offset=0; point_offset<total_points; point_offset+=batch_size )); do
    remaining=$((total_points-point_offset))
    count=$batch_size
    if (( remaining<count )); then count=$remaining; fi
    D2_OFFSET="$point_offset" D2_POINTS="$point_file" \
        sbatch "$@" --array="0-$((count-1))%${concurrency}" cluster_array.slurm
done
