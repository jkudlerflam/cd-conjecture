#!/usr/bin/env python3
"""One priority worker while two registered legacy extensions finish."""
import subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
for n,L,P,c in [('pref_ref',48,0,128),('pref_ref',48,1,128),('pref_ref',96,0,128),('pref_ref',96,1,128),('pref_rat',48,0,128),('pref_rat',48,1,128),('pref_rat',96,0,128),('pref_rat',96,1,128)]:
    lab=f'{n}_L{L}_P{P}_s17';start=time.time();print('START '+lab,flush=True)
    with (HERE/'logs'/f'{lab}.log').open('a') as f:
        rc=subprocess.run([sys.executable,str(HERE/'run.py'),'--name',n,'--L',str(L),'--parity',str(P),'--chi',str(c)],stdout=f,stderr=subprocess.STDOUT).returncode
    print(f'END {lab} rc={rc} seconds={time.time()-start:.1f}',flush=True)
    if rc:raise RuntimeError(lab)
