#!/usr/bin/env python3
"""Bounded local Stage B runs; two single-BLAS-thread processes, no cluster."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--workers',type=int,default=2)
    ap.add_argument('--config',type=Path,default=HERE/'grid.json')
    ap.add_argument('--names',nargs='*')
    ap.add_argument('--sizes',nargs='*',type=int,default=[24,48])
    ap.add_argument('--chi',type=int)
    args=ap.parse_args()
    names=args.names or ['control_xx','control_two_xx','bridge_r1.5_u1.5_k3',
        'bridge_r2_u1.5_k4','bridge_r1.5_u2_k1','bridge_r4_u1_k4',
        'attractive_Jz-0.5','reference_peotta_exploratory']
    tasks=[]
    for name in names:
        for L in args.sizes:
            for parity in ([0] if name.startswith('control') else [0,1]):
                tasks.append((name,L,parity,args.chi or (64 if L<=24 else 96 if L<=48 else 128)))
    logdir=HERE/'logs';logdir.mkdir(exist_ok=True)
    def run(task):
        name,L,P,chi=task
        label=f'{name}_L{L}_P{P}_s17'
        out=HERE/'dmrg'/f'{label}.json'
        cmd=[sys.executable,str(HERE/'dmrg_screen.py'),'--config',str(args.config),
             '--name',name,'--L',str(L),'--parity',str(P),'--chi',str(chi)]
        if out.exists():
            old=json.loads(out.read_text())
            if old.get('status')=='done' and old.get('chi',0)>=chi:
                print('SKIP '+label,flush=True);return
            cmd.append('--resume')
        print('START '+label+f' chi={chi}',flush=True)
        t=time.time()
        with (logdir/f'{label}.log').open('a') as log:
            rc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
        print(f'END {label} rc={rc} seconds={time.time()-t:.1f}',flush=True)
        if rc: raise RuntimeError(f'Failed {label}; inspect its log')
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(run,tasks))


if __name__=='__main__':main()
