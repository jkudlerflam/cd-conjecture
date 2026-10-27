#!/bin/bash
#SBATCH --job-name=d2-qubit
#SBATCH --array=0-1337%48
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=results/slurm_%A_%a.out

# Submit from this directory after transferring d2.py and the generated point file.
# The CPU implementation does not use GPUs. Select an authorized Typhon partition
# with sbatch --partition=...; no partition or Python module is assumed here.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
task_index=$(( ${SLURM_ARRAY_TASK_ID:?} + ${D2_OFFSET:-0} ))
point_file=${D2_POINTS:-results/global_points.json}
point_id=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))[int(sys.argv[2])]["id"])' "$point_file" "$task_index")
python3 d2.py scan --points "$point_file" --index "$task_index" \
    --sizes 14 16 18 --output "results/cluster_scan_${task_index}.json"
python3 d2.py refine --input "results/cluster_scan_${task_index}.json" \
    --ids "$point_id" --sizes 14 16 18 20 22 --levels 16 \
    --output "results/cluster_refined_${task_index}.json"
