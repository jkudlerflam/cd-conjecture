#!/usr/bin/env python3
#SBATCH --job-name=C112-smooth
#SBATCH --array=0-3
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=C112-%A_%a.out

"""Unsubmitted job. Submit from the C112 folder on Typhon.

Partition was not specified because cluster discovery could not connect.
Resource requests are unbenchmarked. Every calculation runs in an allocation.
"""
import os
from pathlib import Path
import subprocess
import sys

if not os.environ.get('SLURM_JOB_ID'):
    raise SystemExit('This script must run as a SLURM job.')
folder = Path(os.environ['SLURM_SUBMIT_DIR'])
suite = ['xx', 'xxz', 'random', 'additional'][int(os.environ['SLURM_ARRAY_TASK_ID'])]
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
subprocess.run([sys.executable, str(folder/'smooth_checks.py'), '--cluster',
                '--suite', suite, '--output', str(folder/f'cluster_{suite}_results.json')],
               cwd=folder, check=True)
