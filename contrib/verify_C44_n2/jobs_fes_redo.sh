#!/bin/bash
# reruns of VUMPS points that did not converge, from a different length/payload split of chi
#SBATCH -J vC44n2fesR
#SBATCH -p typhon
#SBATCH -c 24
#SBATCH --mem=90G
#SBATCH -t 2-00:00:00
#SBATCH -o logs/%x_%j.out
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
cd $SLURM_SUBMIT_DIR
python3 -B -c "import run_fes; run_fes.payload_state(60)" > logs/cache_60.log 2>&1 &
python3 -B -c "import run_fes; run_fes.payload_state(36)" > logs/cache_36.log 2>&1 &
wait
echo "M0 240 4|M1 144 4|M0 96 2|M2 96 2|M1 288 4|M2 288 4" | tr '|' '\n' |
  xargs -P 6 -n 3 sh -c 'python3 -B run_fes.py $0 $1 $2 > logs/fesR_$0_$1_$2.log 2>&1'
