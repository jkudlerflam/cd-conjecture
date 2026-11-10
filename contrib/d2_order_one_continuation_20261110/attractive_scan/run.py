#!/usr/bin/env python3
"""Local full-spin attractive-route calculations, preserving old artifacts."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[k]='1'
import argparse
import fcntl
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PREV=ROOT/'contrib/d2_order_one_20261110'
sys.path.insert(0,str(PREV))
from dmrg_screen import run_point
from ed_screen import coupling_metrics

POINTS={
 'attr_ref':dict(Jx=(math.cos(1)+.7)/2,Jy=(math.cos(1)-.7)/2,Dz=math.sin(1)/2,Jz=-.6,hz=0.),
 'attr_rat':dict(Jx=5/8,Jy=-1/12,Dz=5/12,Jz=-3/5,hz=0.),
}


def _one(name,L,P,chi,seed=17):
    if name == "attr_rat" and (HERE/"defer_attr_rat").exists():
        print("DEFER attr_rat",flush=True);return
    config=HERE/'points.json'
    points=json.loads(config.read_text()) if config.exists() else POINTS
    point=dict(name=name,couplings=points[name])
    out=HERE/'results'/f'{name}_L{L}_P{P}_s{seed}.json'
    if out.exists():
        r=json.loads(out.read_text())
        if r.get('status')=='done' and r.get('chi',0)>=chi:
            print('DONE '+out.name,flush=True);return
        resume=True;stages=[chi]
    else:
        resume=False;stages=sorted(set([16,32,64,128,chi]))
        stages=[x for x in stages if x<=chi]
    run_point(point,L,P,chi,seed,resume,stages=stages,max_sweeps=6,
              output_dir=HERE/'results',checkpoint_dir='/private/tmp/d2-order-one-attractive-checkpoints')


def one(name,L,P,chi,seed=17):
    if (HERE/'defer_preferred_extensions').exists() and ((name=='pref_rat' and L==192) or (name in ('pref_ref','pref_rat') and L==96 and chi>128)):
        print('DEFER preferred extension',flush=True);return
    # New launches share two CPU slots. The two pre-existing legacy192 workers
    # predate this lock; parent authorized one extra worker until they finish.
    if name == 'attr_rat' and (HERE/'defer_attr_rat').exists():
        print('DEFER attr_rat',flush=True);return
    if name == 'pref_rat' and L == 192:
        # Parent priority: finish reference third size before rational extensions.
        def ref_done(p):
            f=HERE/'results'/f'pref_ref_L192_P{p}_s17.json'
            if not f.exists():return False
            r=json.loads(f.read_text());return r.get('status')=='done' and r.get('chi',0)>=192
        while not all(ref_done(p) for p in (0,1)):time.sleep(5)
    locks=HERE/'locks';locks.mkdir(exist_ok=True)
    with (locks/f'{name}_L{L}_P{P}_s{seed}.lock').open('a') as job:
        fcntl.flock(job,fcntl.LOCK_EX)
        slot=None
        while slot is None:
            for i in range(2):
                f=(locks/f'cpu{i}.lock').open('a')
                try:fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB);slot=f;break
                except BlockingIOError:f.close()
            if slot is None:time.sleep(1)
        try:_one(name,L,P,chi,seed)
        finally:slot.close()


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--name');ap.add_argument('--L',type=int);ap.add_argument('--parity',type=int,default=0)
    ap.add_argument('--chi',type=int,default=128);ap.add_argument('--seed',type=int,default=17)
    ap.add_argument('--batch',action='store_true')
    a=ap.parse_args()
    if not (HERE/'points.json').exists():
        (HERE/'points.json').write_text(json.dumps(POINTS,indent=2)+'\n')
    if not a.batch:
        one(a.name,a.L,a.parity,a.chi,a.seed);return
    points=json.loads((HERE/'points.json').read_text())
    (HERE/'coupling_metrics.json').write_text(json.dumps({n:coupling_metrics(p) for n,p in points.items()},indent=2)+'\n')
    # Reference extensions start first; rational controls follow in the freed slots.
    tasks=[('attr_ref',192,0,192),('attr_ref',192,1,192),
           ('attr_rat',48,0,128),('attr_rat',48,1,128),
           ('attr_rat',96,0,128),('attr_rat',96,1,128)]
    (HERE/'logs').mkdir(exist_ok=True)
    def task(t):
        n,L,P,c=t;lab=f'{n}_L{L}_P{P}_s17';start=time.time()
        print('START '+lab,flush=True)
        cmd=[sys.executable,str(Path(__file__).resolve()),'--name',n,'--L',str(L),'--parity',str(P),'--chi',str(c)]
        with (HERE/'logs'/f'{lab}.log').open('a') as f:
            rc=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
        print(f'END {lab} rc={rc} seconds={time.time()-start:.1f}',flush=True)
        if rc:raise RuntimeError(lab)
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(task,tasks))


if __name__=='__main__':main()
