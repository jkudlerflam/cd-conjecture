#!/usr/bin/env python3
"""Bounded, two-worker validation of the registered lead."""
from concurrent.futures import ThreadPoolExecutor
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
from scan import HERE

name='b_r5_2_d1_10_q0_hm7_8'
jobs=[(96,0,192,17),(96,1,192,17),(48,0,192,17),(48,1,192,17)]
parser=argparse.ArgumentParser();parser.add_argument('--job-index',type=int)
args=parser.parse_args()
if args.job_index is not None:jobs=[jobs[args.job_index]]
snapshot=HERE/'validation_snapshots';snapshot.mkdir(exist_ok=True)
for L,P,chi,seed in jobs:
    f=HERE/'dmrg'/f'{name}_L{L}_P{P}_s{seed}.json'
    if f.exists():
        old=json.loads(f.read_text())
        dest=snapshot/f'{f.stem}_chi{old.get("chi","running")}.json'
        if old['status']=='done' and not dest.exists():shutil.copy2(f,dest)

def run(job):
    L,P,chi,seed=job;label=f'{name}_L{L}_P{P}_s{seed}'
    cmd=[sys.executable,str(HERE/'scan.py'),'--single','--config',str(HERE/'stage_b_run.json'),
         '--name',name,'--L',str(L),'--parity',str(P),'--chi',str(chi),'--seed',str(seed)]
    out=HERE/'dmrg'/f'{label}.json'
    if out.exists():
        old=json.loads(out.read_text())
        if old['status']=='done' and old['chi']>=chi:
            print('SKIP',label,chi,flush=True);return
        cmd.append('--resume')
    print('START',label,chi,flush=True);start=time.time()
    with (HERE/'logs'/f'{label}.log').open('a') as f:
        rc=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
    print('DONE',label,chi,'rc',rc,'seconds',round(time.time()-start,2),flush=True)
    if rc:raise RuntimeError(label)

with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(run,jobs))
