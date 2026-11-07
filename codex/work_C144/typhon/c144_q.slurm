#!/bin/bash
#SBATCH --job-name=c144_q
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c144/slurm_%j.out
cd /scratch/lustre/jkudlerflam/c144
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# warm the numba cache once, then run 12 tasks in parallel with 4 threads each
/opt/sns/anaconda3/2024.06/bin/python3 -s -B -c "from xzlib import term, ground_energy; print(ground_energy(term('xz'), 12))"
xargs -P 12 -I{} ./run_q.sh {} < qtasks.txt
echo done
