#!/bin/bash
#SBATCH --job-name=c155adv
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=11:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c155/logs/c155_%j.out
# C155 packed job (Claude Opus agent, 2026-11-07). All stages on one 48-core node, every stage resumable.
#   E  directed XZ-neighbourhood + Kitaev [110]-field scan with DMRG to L = 128/256 (xz/pt_*.json)
#   A  uniform volume sample of S^8: 48 x 250 points            (results_crit/uniform_*.jsonl)
#   B  gap-closing minimization + ball test + slide: 16 tasks    (results_crit/gap_*.jsonl)
#   C  c maximization (Nelder-Mead restarts): 32 tasks           (results/task_*.jsonl)
#   D  refinement of candidates (ED 16/20, DMRG 32/64/128)       (refine/cand_*.json)
cd /scratch/lustre/jkudlerflam/c155
mkdir -p results results_crit refine logs xz
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
HOURS=${HOURS:-4}
timeout 2h python3 xz_scan.py 48 128 > logs/E.log 2>&1
seq 0 47 | xargs -P 48 -I{} sh -c 'python3 crit.py --stage A --task {} --npts 250 > logs/A_{}.log 2>&1'
( for t in $(seq 0 15); do echo "B $t"; done; for t in $(seq 0 31); do echo "C $t"; done ) | \
  xargs -P 48 -L 1 sh -c 'if [ "$0" = B ]; then python3 crit.py --stage B --task $1 --hours '"$HOURS"' > logs/B_$1.log 2>&1; else python3 search.py --task $1 --hours '"$HOURS"' --maxfev 300 > logs/C_$1.log 2>&1; fi'
timeout 4h python3 refine.py --thr 1.1 --top 24 --chi 128 --procs 48 > logs/D.log 2>&1
python3 summarize.py > summary.txt 2>&1
