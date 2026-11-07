#!/bin/bash
#SBATCH --job-name=cdfock1
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=16:00:00
#SBATCH --output=logs/typhon_job1_%j.out
# Packed job: crossed spectra at m = 8, 9 and insertion analyses at m = 8 (contrib/cd_crossed_fock).
cd /scratch/lustre/jkudlerflam/cd_crossed_fock
mkdir -p data logs
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
MMIN=8 python3 run_spectra.py 0.5 8 > logs/t_D0.5_m8.log 2>&1 &
MMIN=8 python3 run_spectra.py 1 8 > logs/t_D1_m8.log 2>&1 &
MMIN=8 python3 run_spectra.py 2 8 > logs/t_D2_m8.log 2>&1 &
MMIN=7 python3 run_spectra_model.py TFIM 8 > logs/t_TFIM.log 2>&1 &
MMIN=7 python3 run_spectra_model.py XYI 8 > logs/t_XYI.log 2>&1 &
MMIN=7 python3 run_spectra_model.py XYZ 8 > logs/t_XYZ.log 2>&1 &
MMIN=8 python3 run_spectra_model.py XXZSTAG 8 > logs/t_XXZSTAG.log 2>&1 &
MMIN=9 python3 run_spectra.py 1 9 2 4 6 8 9 > logs/t_D1_m9.log 2>&1 &
MMIN=9 python3 run_spectra.py 0 9 2 4 6 8 9 > logs/t_D0_m9.log 2>&1 &
python3 analysis_insertion.py 8 2,4,6,8 XX,XXZ05,HEIS,XXZ2,XXZSTAG > logs/t_ins8.log 2>&1 &
(for b in 4 8; do for n in HEIS XX XXZ05 XXZ2; do python3 two_insertion.py $n $b 8 8; done; done) > logs/t_two8.log 2>&1 &
python3 analysis_insertion.py 7 2,4,6,8 TFIM,XYI,XYZ > logs/t_ins7z2.log 2>&1 &
wait
