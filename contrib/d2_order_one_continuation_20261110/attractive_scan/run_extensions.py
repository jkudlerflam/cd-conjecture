#!/usr/bin/env python3
"""Registered conditional promotions, starting after both earlier queues finish."""
from concurrent.futures import ThreadPoolExecutor
import json,math,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
RESULTS=HERE/'results'
if (HERE/'defer_preferred_extensions').exists():
    print('DEFER: '+(HERE/'defer_preferred_extensions').read_text(),flush=True);sys.exit(0)
def read(n,L,P):
    p=RESULTS/f'{n}_L{L}_P{P}_s17.json'
    if not p.exists():return None
    r=json.loads(p.read_text())
    return r if r.get('status')=='done' else None
required=[('attr_ref',192,P) for P in (0,1)]+[(n,L,P) for n in ('pref_ref','pref_rat') for L in (48,96) for P in (0,1)]
while not all(read(*t) is not None for t in required):time.sleep(10)
print('Earlier queues complete; evaluating registered promotion gates.',flush=True)
tasks=[]
for n in ('pref_rat','pref_ref'):
    selected=[min((read(n,L,P) for P in (0,1)),key=lambda r:r['E_physical']) for L in (48,96)]
    a,b=selected
    slopes={k:6*(b[k]-a[k])/math.log(2) for k in ('S_central3','S_fractional')}
    promote=all(v>1 for v in slopes.values()) and all(r['convergence']['sweep_converged'] for r in selected)
    print(json.dumps(dict(name=n,slopes=slopes,promote=promote)),flush=True)
    if promote:tasks +=[(n,192,P,192) for P in (0,1)]+[(n,96,P,192) for P in (0,1)]
def task(t):
    n,L,P,c=t;lab=f'{n}_L{L}_P{P}_s17';start=time.time();print('START '+lab+f' chi{c}',flush=True)
    with (HERE/'logs'/f'{lab}.log').open('a') as f:
        rc=subprocess.run([sys.executable,str(HERE/'run.py'),'--name',n,'--L',str(L),'--parity',str(P),'--chi',str(c)],stdout=f,stderr=subprocess.STDOUT).returncode
    print(f'END {lab} chi{c} rc={rc} seconds={time.time()-start:.1f}',flush=True)
    if rc:raise RuntimeError(lab)
with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(task,tasks))
