#!/bin/bash
#SBATCH --job-name=cdfock2
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=300G
#SBATCH --time=12:00:00
#SBATCH --output=logs/typhon_job2_%j.out
# Packed job: crossed vacuum energy eps0(beta) at Trotter numbers 9..13 for Richardson extrapolation.
cd /scratch/lustre/jkudlerflam/cd_crossed_fock
mkdir -p data logs
export OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6
B=1,2,3,4,6,8
for n in XX XXZ05 HEIS XXZ2 TFIM XYI XYZ XXZSTAG; do
  (for m in 9 10 11 12 13; do python3 eps0_lanczos.py $n $m $B; done) > logs/t_eps0_$n.log 2>&1 &
done
wait
