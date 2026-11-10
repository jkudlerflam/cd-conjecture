#!/usr/bin/env python3
"""Preserve the completed largest-size state, then refine its chi."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
from scan import HERE,D

name='b_r5_2_d1_10_q0_hm7_8'
rr=[json.loads((HERE/'dmrg'/f'{name}_L192_P{P}_s17.json').read_text()) for P in (0,1)]
assert all(r['status']=='done' and r['convergence']['sweep_converged'] for r in rr)
r=min(rr,key=lambda x:x['E_physical']);P=r['parity']
assert r['chi']==192,'This preservation step expects the original chi192 state'
label=f'{name}_L192_P{P}_s17'
source=HERE/'dmrg'/f'{label}.json';checkpoint=Path(r['checkpoint'])
destination=HERE/'validation_snapshots'/f'{label}_chi192.json'
ckdest=checkpoint.parent/'validation_snapshots'/f'{label}_chi192_mps.npz'
destination.parent.mkdir(exist_ok=True);ckdest.parent.mkdir(exist_ok=True)
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()
before=[sha(source),sha(checkpoint)]
assert not destination.exists() and not ckdest.exists(),'Do not overwrite an existing preserved state'
shutil.copy2(source,destination);shutil.copy2(checkpoint,ckdest)
after=[sha(source),sha(checkpoint)]
copied=[sha(destination),sha(ckdest)]
assert before==after==copied,'Source changed during preservation'
D.atomic_json(HERE/'largest_size_preservation.json',dict(name=name,L=192,parity=P,chi=192,
    raw_record=str(destination),raw_record_sha256=copied[0],
    original_checkpoint=str(checkpoint),preserved_checkpoint=str(ckdest),checkpoint_sha256=copied[1],
    registration='L192_CHI_VALIDATION.md'))
cmd=[sys.executable,str(HERE/'scan.py'),'--single','--config',str(HERE/'stage_b_run.json'),
     '--name',name,'--L','192','--parity',str(P),'--chi','256','--resume']
print('PRESERVED',label,'START chi256',flush=True)
with (HERE/'logs'/f'{label}.log').open('a') as f:rc=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
print('DONE',label,'chi256 rc',rc,flush=True)
raise SystemExit(rc)
