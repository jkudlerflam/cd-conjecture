#!/bin/bash
#SBATCH --job-name=cdfock3
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=logs/typhon_job3_%j.out
# Smaller packed job (backfill-friendly): m = 8 crossed spectra and insertion analyses for the U(1) d = 2 chains.
cd /scratch/lustre/jkudlerflam/cd_crossed_fock
mkdir -p data logs
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
MMIN=8 python3 run_spectra.py 0.5 8 0.5 1 2 3 4 6 8 > logs/t3_D0.5_m8.log 2>&1 &
MMIN=8 python3 run_spectra.py 1 8 0.5 1 2 3 4 6 8 > logs/t3_D1_m8.log 2>&1 &
MMIN=8 python3 run_spectra.py 2 8 0.5 1 2 3 4 6 8 > logs/t3_D2_m8.log 2>&1 &
MMIN=8 python3 run_spectra_model.py XXZSTAG 8 0.5 1 2 3 4 6 8 > logs/t3_XXZSTAG_m8.log 2>&1 &
python3 analysis_insertion.py 8 4,8 XX,XXZ05,XXZ2 > logs/t3_ins8.log 2>&1 &
python3 analysis_insertion.py 8 8 HEIS > logs/t3_ins8h.log 2>&1 &
wait
