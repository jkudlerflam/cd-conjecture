#!/bin/bash
#SBATCH -J rt_sel
#SBATCH -p typhon
#SBATCH -c 8
#SBATCH --mem=60G
#SBATCH -t 12:00:00
#SBATCH -a 0-4
#SBATCH -o /scratch/lustre/jkudlerflam/redteam_routeA_cells/sel_%a.out
cd /scratch/lustre/jkudlerflam/redteam_routeA_cells
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
case $SLURM_ARRAY_TASK_ID in
 0) python3 -B -u rt_selector.py full ;;
 1) python3 -B -u rt_selector.py comp 12 0.0 0.05 ;;
 2) python3 -B -u rt_selector.py comp 12 0.0 0.01 ;;
 3) python3 -B -u rt_selector.py comp 12 0.3 0.05 ;;
 4) python3 -B -u rt_selector.py ring ;;
esac
