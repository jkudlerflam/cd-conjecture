"""Run a list of ring_boot.py / ed.py tasks in parallel inside one SLURM allocation."""
import os, sys, json, subprocess, itertools
from concurrent.futures import ThreadPoolExecutor
if 'SLURM_JOB_ID' not in os.environ:
    sys.exit('run inside a SLURM allocation')
tasks = [l.strip() for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
ncpu = int(os.environ.get('NPAR', os.environ.get('SLURM_CPUS_PER_TASK', '4')))
def run(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return cmd, r.returncode, r.stdout[-400:], r.stderr[-600:]
with ThreadPoolExecutor(ncpu) as ex:
    for cmd, rc, out, err in ex.map(run, tasks):
        print('RC', rc, cmd, flush=True)
        if rc:
            print(err, flush=True)
