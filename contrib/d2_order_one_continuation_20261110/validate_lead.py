#!/usr/bin/env python3
"""Separate checkpoint continuations and third-size tests of the R=28 lead."""
from pathlib import Path
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
import shutil
import subprocess
import sys

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'d2_order_one_20261110'
sys.path.insert(0,str(OLD))
import dmrg_screen as D

NAME='bridge_r4_u1_k4'
OUT=HERE/'validation'
CK=Path('/private/tmp/d2-order-one-continuation-checkpoints/validation')
POINT=dict(name=NAME,couplings=dict(Jx=.5,Jy=-.5,Jz=1.,Dz=1/16,hz=-7/4))


def run(L,P):
    OUT.mkdir(parents=True,exist_ok=True);CK.mkdir(parents=True,exist_ok=True)
    label=f'{NAME}_L{L}_P{P}_s17'
    target=OUT/f'{label}.json'
    if L==96 and not target.exists():
        old=json.loads((OLD/'dmrg'/f'{label}.json').read_text())
        dest=CK/f'{label}_mps.npz'
        shutil.copy2(old['checkpoint'],dest)
        old['input_record']=str(OLD/'dmrg'/f'{label}.json')
        old['checkpoint']=str(dest)
        old['status']='queued_input'
        D.atomic_json(target,old)
    resume=target.exists()
    chi=256 if L==96 else 192
    if resume:
        old=json.loads(target.read_text())
        if old.get('status')=='done' and old.get('chi',0)>=chi:
            print(f'SKIP {label}',flush=True);return
    D.run_point(POINT,L,P,chi=chi,seed=17,resume=resume,
        stages=[chi] if resume else [16,32,64,128,192],
        max_sweeps=6 if L==96 else 8,output_dir=OUT,checkpoint_dir=CK)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--L',type=int);ap.add_argument('--parity',type=int)
    ap.add_argument('--batch',action='store_true');ap.add_argument('--workers',type=int,default=3)
    args=ap.parse_args()
    if not args.batch:run(args.L,args.parity);return
    logs=HERE/'logs';logs.mkdir(exist_ok=True)
    def task(pair):
        L,P=pair
        print(f'START L{L} P{P}',flush=True)
        with (logs/f'lead_L{L}_P{P}.log').open('a') as log:
            rc=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--L',str(L),'--parity',str(P)],
                              stdout=log,stderr=subprocess.STDOUT).returncode
        print(f'END L{L} P{P} rc={rc}',flush=True)
        if rc:raise RuntimeError((L,P,rc))
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(task,[(96,0),(96,1),(192,0),(192,1)]))


if __name__=='__main__':main()
