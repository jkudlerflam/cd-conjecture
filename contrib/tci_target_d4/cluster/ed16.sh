#!/bin/bash
#SBATCH --partition=typhon
#SBATCH --qos=short
#SBATCH --output=logs/%x_%j.out
# usage: sbatch --job-name=X --mem=M --time=T ed16.sh K NLEV [--entropy]
cd /scratch/lustre/jkudlerflam/tci_d4/tci_target_d4
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
/usr/bin/time -v python3 -u ed_lean.py --L 16 --points cluster/l16_points.json --momenta $1 --nlev $2 $3 --out data/l16_k$1.jsonl 2>&1 | grep -v "^\s*\(Page\|Socket\|Signals\|Swaps\|File\|Average\|Exit\|Minor\|Major\|Volun\|Invol\)"
