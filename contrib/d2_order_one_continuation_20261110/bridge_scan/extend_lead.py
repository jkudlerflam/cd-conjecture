#!/usr/bin/env python3
"""Conditional extension registered in LARGE_SIZE_EXTENSION.md."""
from concurrent.futures import ThreadPoolExecutor
import json
import math
import shutil
import subprocess
import sys
import time
from scan import HERE, D

name='b_r5_2_d1_10_q0_hm7_8'
by_size={L:[json.loads((HERE/'dmrg'/f'{name}_L{L}_P{P}_s17.json').read_text()) for P in (0,1)] for L in (48,96)}
assert all(x['status']=='done' and x['convergence']['sweep_converged'] for rr in by_size.values() for x in rr)
best={L:min(rr,key=lambda x:x['E_physical']) for L,rr in by_size.items()}
slopes={k:6*(best[96][k]-best[48][k])/math.log(2) for k in ['S_fractional','S_central3']}
assert min(slopes.values())>1.05,slopes
selected=best[96]['parity']
jobs=[(96,selected,256),(192,0,192),(192,1,192),(48,0,192),(48,1,192)]
plan=HERE/'large_size_jobs.json'
if not plan.exists():D.atomic_json(plan,dict(name=name,jobs=jobs,selected_L96_parity=selected,trigger_slopes=slopes,registration='LARGE_SIZE_EXTENSION.md'))
snap=HERE/'validation_snapshots';snap.mkdir(exist_ok=True)
for rr in by_size.values():
    for x in rr:
        f=HERE/'dmrg'/f'{name}_L{x["L"]}_P{x["parity"]}_s17.json'
        target=snap/f'{f.stem}_chi{x["chi"]}.json'
        if not target.exists():shutil.copy2(f,target)

def run(job):
    L,P,chi=job;label=f'{name}_L{L}_P{P}_s17';out=HERE/'dmrg'/f'{label}.json'
    cmd=[sys.executable,str(HERE/'scan.py'),'--single','--config',str(HERE/'stage_b_run.json'),
         '--name',name,'--L',str(L),'--parity',str(P),'--chi',str(chi)]
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
