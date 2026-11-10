#!/usr/bin/env python3
"""Registered simple preferred reserve, one sequential numerical worker."""
import json,time,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent

def done(n,L,P):
    f=HERE/'results'/f'{n}_L{L}_P{P}_s17.json'
    return f.exists() and json.loads(f.read_text()).get('status')=='done'
while not all(done(n,L,P) for n in ('pref_ref','pref_rat') for L in (48,96) for P in (0,1)):time.sleep(5)
for L,P in ((48,0),(48,1),(96,0),(96,1)):
    lab=f'pref_zrat_L{L}_P{P}_s17';print('START '+lab,flush=True);start=time.time()
    with (HERE/'logs'/f'{lab}.log').open('a') as f:
        rc=subprocess.run([sys.executable,str(HERE/'run.py'),'--name','pref_zrat','--L',str(L),'--parity',str(P),'--chi','128'],stdout=f,stderr=subprocess.STDOUT).returncode
    print(f'END {lab} rc={rc} seconds={time.time()-start:.1f}',flush=True)
    if rc:raise RuntimeError(lab)
