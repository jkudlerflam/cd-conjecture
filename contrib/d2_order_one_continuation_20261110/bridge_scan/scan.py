#!/usr/bin/env python3
"""Registered full-qubit scan, importing the unchanged validated wrapper."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PRIOR=ROOT/'contrib/d2_order_one_20261110'
sys.path.insert(0,str(PRIOR))
import dmrg_screen as D
from ed_screen import coupling_metrics


def safe_fraction(v):
    return str(Fraction(v)).replace('-','m').replace('/','_')


def point(r,d,q,h):
    r,d,q,h=map(Fraction,(r,d,q,h))
    p=dict(Jx=float((1+q)/2),Jy=float((q-1)/2),Jz=float(r/4),Dz=float(d),hz=float(h))
    return dict(name=f'b_r{safe_fraction(r)}_d{safe_fraction(d)}_q{safe_fraction(q)}_h{safe_fraction(h)}',
        family=f'r{safe_fraction(r)}_d{safe_fraction(d)}',r=str(r),D=str(d),q=str(q),hz=str(h),
        couplings=p,metrics=coupling_metrics(p),q_cancellation_fraction=float(abs(q)))


def make_grid():
    families=[('2','3/16',['-3/4','-11/16','-5/8','-9/16']),
              ('2','1/4',['-7/8','-13/16','-3/4','-11/16']),
              ('5/2','1/10',['-15/16','-7/8','-27/32','-13/16'])]
    points=[point(r,d,q,h) for r,d,fields in families for q in ['0','1/16','1/8'] for h in fields]
    points += [dict(name='control_xx',couplings=dict(Jx=-.5,Jy=-.5,Jz=0.,Dz=0.,hz=0.)),
               dict(point('2','3/16','0','-1/2'),name='control_prior_seed')]
    D.atomic_json(HERE/'grid.json',dict(points=points,sizes=[24],parities=[0,1],chi=64))


def run_one(args):
    control=HERE/'queue_control.json'
    if control.exists():
        deferred=json.loads(control.read_text()).get('defer',[])
        if any(x['name']==args.name and x['L']==args.L for x in deferred):
            print('DEFERRED by queue_control.json; no numerical output created',args.name,args.L,flush=True)
            return
    cfg=json.loads(args.config.read_text())
    p=next(p for p in cfg['points'] if p['name']==args.name)
    D.run_point(p,args.L,args.parity,args.chi,args.seed,args.resume,
        output_dir=HERE/'dmrg',checkpoint_dir='/private/tmp/d2-order-one-continuation-bridge',
        max_sweeps=args.max_sweeps)


def batch(args):
    cfg=json.loads(args.config.read_text())
    points=[p for p in cfg['points'] if not args.names or p['name'] in args.names]
    tasks=[(p,L,P) for p in points for L in args.sizes for P in ([0] if p['name']=='control_xx' else [0,1])]
    logdir=HERE/'logs';logdir.mkdir(exist_ok=True)
    def launch(task):
        p,L,P=task;name=p['name'];label=f'{name}_L{L}_P{P}_s{args.seed}'
        out=HERE/'dmrg'/f'{label}.json'
        cmd=[sys.executable,str(Path(__file__).resolve()),'--single','--config',str(args.config),
            '--name',name,'--L',str(L),'--parity',str(P),'--chi',str(args.chi),'--seed',str(args.seed),
            '--max-sweeps',str(args.max_sweeps)]
        if out.exists():
            old=json.loads(out.read_text())
            assert old['couplings']==p['couplings']
            if old['status']=='done' and old['chi']>=args.chi:
                print('SKIP',label,flush=True);return
            cmd.append('--resume')
        start=time.time();print('START',label,flush=True)
        with (logdir/(label+'.log')).open('a') as log:
            rc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
        print('DONE',label,'rc',rc,'seconds',round(time.time()-start,2),flush=True)
        if rc:raise RuntimeError(label)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:list(pool.map(launch,tasks))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--make-grid',action='store_true')
    ap.add_argument('--single',action='store_true');ap.add_argument('--config',type=Path,default=HERE/'grid.json')
    ap.add_argument('--names',nargs='*');ap.add_argument('--name');ap.add_argument('--L',type=int,default=24)
    ap.add_argument('--sizes',type=int,nargs='*',default=[24]);ap.add_argument('--parity',type=int,default=0)
    ap.add_argument('--chi',type=int,default=64);ap.add_argument('--seed',type=int,default=17)
    ap.add_argument('--workers',type=int,default=2);ap.add_argument('--resume',action='store_true')
    ap.add_argument('--max-sweeps',type=int,default=6)
    args=ap.parse_args()
    if args.make_grid:make_grid()
    elif args.single:run_one(args)
    else:batch(args)


if __name__=='__main__':main()
