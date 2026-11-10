#!/bin/bash
#SBATCH --job-name=c193
#SBATCH --partition=apollo
#SBATCH --qos=gpu-medium
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --nodelist=apollo02
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c19x/c193/logs/c193_%j.out
# C193b (Claude Opus 5.5 agent, 2026-11-10).  Correlators on stored C190 states, one GPU, complex128.
# Usage: sbatch c193.slurm [<state list> [<outdir>]]   (defaults: states_c193.txt, out)
# The list has one "<run.json> <ckpt.npz>" per line; finished outputs (<outdir>/<tag>_meas.json, status done)
# are skipped, so a resubmission resumes.  The checkpoint trees are only read.
cd /scratch/lustre/jkudlerflam/c19x/c193
LIST=${1:-states_c193.txt}
OUT=${2:-out}
export CUDA_PATH=/scratch/lustre/jkudlerflam/c18x/c190g/cudahdr
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
echo "GPU $CUDA_VISIBLE_DEVICES on $(hostname) $(date)"
nvidia-smi -L
python3 -c "import cupy; print('cupy', cupy.__version__)"
python3 meas193.py batch "$LIST" "$OUT"
echo "end $(date)"
