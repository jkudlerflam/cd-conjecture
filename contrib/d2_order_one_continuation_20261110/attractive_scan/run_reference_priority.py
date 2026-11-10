#!/usr/bin/env python3
"""Begin preferred third size as legacy slots free, with global CPU locks."""
import json,time,subprocess,sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
HERE=Path(__file__).resolve().parent

def old_done(P):
    f=HERE/'results'/f'attr_ref_L192_P{P}_s17.json'
    return f.exists() and json.loads(f.read_text()).get('status')=='done'
def task(P):
    while sum(old_done(p) for p in (0,1))<(1 if P==0 else 2):time.sleep(5)
    lab=f'pref_ref_L192_P{P}_s17';print('START '+lab,flush=True);start=time.time()
    with (HERE/'logs'/f'{lab}.log').open('a') as f:
        rc=subprocess.run([sys.executable,str(HERE/'run.py'),'--name','pref_ref','--L','192','--parity',str(P),'--chi','192'],stdout=f,stderr=subprocess.STDOUT).returncode
    print(f'END {lab} rc={rc} seconds={time.time()-start:.1f}',flush=True)
    if rc:raise RuntimeError(lab)
with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(task,(0,1)))
