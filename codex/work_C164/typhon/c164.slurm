#!/bin/bash
#SBATCH --job-name=c164
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=11:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c164/logs/c164_%j.out
# C164 packed job (Claude Opus agent, 2026-11-08). One 48-core node.
#   controls (7 cores, whole job): DMRG c at XX (c = 1, L = 64, 128, 256), two decoupled XX chains as d = 4
#            super-sites (c = 2 control outside B_2, L = 32, 64), two-site VUMPS FES for both     -> ctrl/
#   stage 1 (41 cores, 3 h): 31 adversarial two-valley optimizer tasks (qopt.py; hF, dm, S^8)  -> qopt/
#            10 classical lattice soliton screens (solscreen.py; N eps .3 .5 .8 1.2, T eps .5 1, dm x2, hF, S^8) -> cl/
#   stage 2 (41 cores, 2 h): continuation of every qopt result with O < 0.3 to its onset (qtune.py) -> qtune/
#   stage 3 (41 cores, 4 h): the 41 best tuned endpoints (make_endpoints.py), past-onset open-chain DMRG at
#            L = 64 over 7 valley splittings, L = 128 where c_prof(64) > 1.2 (onset_dmrg.py)     -> onset/
#   summary.txt (summarize.py)
cd /scratch/lustre/jkudlerflam/c164
mkdir -p logs qopt cl qtune onset ctrl
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cat controls.txt | xargs -P 7 -I CMD sh -c 'CMD' &
CTRL=$!
cat stage1.txt | xargs -P 41 -I CMD sh -c 'CMD'
cat qopt/q_*.jsonl > qopt_all.jsonl
seq 0 40 | xargs -P 41 -I{} sh -c 'timeout 2h python3 qtune.py qopt_all.jsonl qtune/t_{}.jsonl 0.3 {}/41 > logs/qtune_{}.log 2>&1'
python3 make_endpoints.py endpoints.jsonl 41 qtune/t_*.jsonl > logs/endpoints.log 2>&1
N=$(wc -l < endpoints.jsonl)
if [ "$N" -gt 0 ]; then
  seq 0 $((N-1)) | xargs -P 41 -I{} sh -c 'timeout 4h python3 onset_dmrg.py endpoints.jsonl {} onset/e_{}.json 0.05 64,128 7 > logs/onset_{}.log 2>&1'
fi
python3 summarize.py . > summary.txt 2>&1
wait $CTRL
python3 summarize.py . > summary.txt 2>&1
python3 fesfit.py ctrl/fes_xx.json ctrl/fes_twoxx.json > ctrl/fesfit.out 2>&1
