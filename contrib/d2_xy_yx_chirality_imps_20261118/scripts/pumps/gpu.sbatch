#!/bin/bash
#SBATCH --partition=apollo
#SBATCH --qos=gpu-long
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=4-00:00:00
#SBATCH --output=logs/%x_%j.out
# packed puMPS runs on one GPU: each line of the list file is "<tag> <Dlist> <pumps args...>"
cd "$SLURM_SUBMIT_DIR"
export CUDA_PATH=${CUDA_PATH:-/path/to/cuda} OMP_NUM_THREADS=8 PUMPS_BACKEND=cupy
export CUPY_CACHE_DIR=$SLURM_SUBMIT_DIR/cupy_cache_$(hostname -s); mkdir -p "$CUPY_CACHE_DIR"
echo "GPU $CUDA_VISIBLE_DEVICES on $(hostname) $(date)"; nvidia-smi --query-gpu=name,memory.total --format=csv
timeout 600 python3 pumps.py model=chir delta=0.1875 h=1.783 L=16 Dsched=8,16,24 ps=-2:2 ns=2 nev=8 save=out/gpusmoke > out/gpusmoke.txt 2>&1 || { echo SMOKE FAILED; tail out/gpusmoke.txt; exit 3; }
grep "^D=24" out/gpusmoke.txt
while read -r line; do [ -z "$line" ] && continue; bash runp.sh $line; done < "$1"
