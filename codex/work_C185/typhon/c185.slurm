#!/bin/bash
#SBATCH --job-name=c185
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=300G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c185/logs/c185_%j.out
# C185 packed job (Claude Opus 5.5 agent, 2026-11-09). Four streams in parallel:
# (1) 10 runs of the U(1) part of the transverse-DM cusp model (C176 fibdmrg, V2 = -2 + v), 1 thread each;
# (2) 12 runs of the full effective model with its dN = 1 term (blockk.py, 3-site blocks), 1 thread each;
# (3) 4 full-qubit-chain runs at V/Delta = 8 with and without the DM term (C177 fc.py copy), 2 threads each;
# (4) 22 degeneracy-search tasks (degen_opt.py), 1 thread each, 10 at a time.
cd /scratch/lustre/jkudlerflam/c185
mkdir -p logs out outbk oc degen
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
(cat typhon/dmrg_tasks.txt | xargs -P 10 -I CMD sh -c 'CMD') &
(cat typhon/block_tasks.txt | xargs -P 12 -I CMD sh -c 'CMD') &
(export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2; cat typhon/full_tasks.txt | xargs -P 4 -I CMD sh -c 'CMD') &
(cat typhon/degen_tasks.txt | xargs -P 10 -I CMD sh -c 'CMD') &
wait
(cd eff && python3 analyze.py ../out > ../summary_eff.txt 2>&1)
