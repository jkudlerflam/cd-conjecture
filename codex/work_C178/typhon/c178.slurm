#!/bin/bash
#SBATCH --job-name=c178
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=300G
#SBATCH --time=20:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c178/logs/c178_%j.out
# C178 packed job (Claude Opus 5.5 agent, 2026-11-09). Twist-resolved, valley-labelled ring ED (spec2.py, kring.py)
# of C173's model (3.2) in the blockaded space. One heavy task (L = 50, N = 10, 24 threads) and a medium queue
# (6 workers x 4 threads). Outputs out/<tag>.jsonl and out/<tag>.out, completion log logs/done.txt.
cd /scratch/lustre/jkudlerflam/c178
mkdir -p out logs
export PYTHONNOUSERSITE=1 NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
/opt/sns/anaconda3/2024.06/bin/python3 -s -B spec2.py 16 2 1.5707963267948966 1.0 cusp 4 4 half out/warm.jsonl > out/warm.out 2>&1
(sed 's/^/24 /' heavy.txt | xargs -P 1 -L 1 ./run_task.sh) &
(sed 's/^/4 /' medium.txt | xargs -P 6 -L 1 ./run_task.sh) &
wait
echo all done $(date) >> logs/done.txt
