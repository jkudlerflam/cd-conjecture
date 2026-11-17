#!/usr/bin/env python3
"""Authenticate a C195 checkpoint and prepare a separate, fail-closed continuation.

No network access or job submission. Full verification can itself be expensive.
The original run/checkpoint is only read. Use the unmodified C195 driver afterwards.
"""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'codex/work_C195'))
import z2dmrg as Z
import models as Mo


def digest(path):
    h = hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):
            h.update(b)
    return h.hexdigest()


def storage(L, chi, itemsize=8, w=5):
    D=[1]+[min(chi,2**min(j,L-j)) for j in range(1,L)]+[1]
    mps=sum(2*a*b for a,b in zip(D,D[1:]))*itemsize
    env=sum(w*d*d for d in D)*itemsize
    # Dense MPS/env are actually stored despite parity blocks in local solves.
    return dict(MPS_GB=mps/1e9, environments_GB=env/1e9,
                static_GB=(mps+env)/1e9, with_25pct_headroom_GB=1.25*(mps+env)/1e9,
                checkpoint_atomic_pair_GB=2*mps/1e9)


def verify(source, checkpoint):
    r=json.loads(Path(source).read_text()); cfg=r['cfg']; L=cfg['L']
    if r['status']!='done':
        raise ValueError('source run is not done; completion alone is not a convergence certificate')
    with np.load(checkpoint,allow_pickle=False) as z:
        if {k for k in z.files if k.startswith('M')}!={f'M{j}' for j in range(L)}:
            raise ValueError('checkpoint tensor count differs')
        if int(z['P'])!=cfg['P']:
            raise ValueError('checkpoint parity differs')
        if json.loads(str(z['meta'])).get('sweeps')!=len(r['hist']):
            raise ValueError('checkpoint sweep count differs')
        psi=Z.MPS([z[f'M{j}'] for j in range(L)],list(map(int,z['n0'])),int(z['P']))
    if not all(np.isfinite(m).all() for m in psi.M) or not psi.check(1e-12):
        raise ValueError('nonfinite tensors or forbidden parity blocks')
    if psi.dims!=r['dims']:
        raise ValueError('checkpoint bond dimensions differ')
    if cfg['model'] in ('xx1','xx2'):
        mpo,gph=Mo.xx_control(L,cfg['mu'],1 if cfg['model']=='xx1' else 2); shift=0.
    else:
        mpo,gph,shift,_=Mo.candidate(cfg['model'],L,cfg['hz'],cfg.get('gauge'))
    if abs(mpo.const-L*shift-r['E_off'])>1e-10:
        raise ValueError('energy-offset or model convention mismatch')
    env=np.zeros((1,mpo.w,1),dtype=complex); env[0,0,0]=1.
    for j in range(L):
        env=Z.env_left(env,psi.M[j],mpo.W[j],psi.n0[j],psi.n0[j+1],mpo)
    E=float(env[0,2,0].real)+r['E_off']
    S=Z.entropies(psi)
    meas=Z.measure(psi,mpo,gph=gph,gwin=min(4,L-4),cp_rmax=2)
    d=dict(dE=abs(E-r.get('E_phys',r['hist'][-1]['E_phys'])),
           dS_max=float(np.max(np.abs(S-r['S']))),
           ddens_max=float(np.max(np.abs(np.array(meas['dens'])-r['meas']['dens']))),
           parity=Z.parity_expectation(psi),L=L,P=psi.P,Dmax=max(psi.dims),
           checkpoint_sha256=digest(checkpoint),source_sha256=digest(source),
           tensor_bytes=sum(m.nbytes for m in psi.M), status='verified')
    if d['dE']>1e-8 or d['dS_max']>1e-7 or d['ddens_max']>1e-7 or abs(d['parity']-(-1)**psi.P)>1e-10:
        raise ValueError('checkpoint measurements do not match: '+json.dumps(d))
    return r,d


def prepare(source,checkpoint,chi,out,*,test=False):
    out=Path(out).resolve(); checkpoint=Path(checkpoint).resolve(); source=Path(source).resolve()
    if out.exists():
        raise ValueError('refusing to overwrite configuration')
    # Perform provenance checks before producing anything executable.
    r,v=verify(source,checkpoint); c=copy.deepcopy(r['cfg']); old=r['hist'][-1]
    if chi<=old['chi']:
        raise ValueError('new chi must exceed source chi')
    if not test:
        if c['model']!='sym40' or c.get('gauge')!='stag' or c['hz']!=-19.976047250613814:
            raise ValueError('not the registered microscopic model/field/gauge')
        if c['L'] not in (256,512) or (c['L']==512 and chi not in (2048,2560)) or (c['L']==256 and chi!=1536):
            raise ValueError('not a registered stage')
        if c['L']==512 and chi==2560 and old['chi']!=2048:
            raise ValueError('must converge a 2048 source stage first')
        # Reject nonstationary or capped-but-unconverged source stages.
        hs=[x for x in r['hist'] if x['chi']==old['chi'] and x['alpha']==0]
        if len(hs)<2 or abs(hs[-1]['E']-hs[-2]['E'])>4e-9 or abs(hs[-1]['Sbar6']-hs[-2]['Sbar6'])>5e-4:
            raise ValueError('source is not stationary')
    c['tag']=c['tag'].split('_nextchi')[0]+f'_nextchi{chi}'
    c['init']='file:'+str(checkpoint); c['sched']=[]; c.pop('conv',None)
    c['final']=dict(chi=chi,res_tol=1e-9,min_sweeps=3,max_sweeps=8,
                    tolE=4e-9 if c['L']==512 else 2e-9,tolS=5e-4)
    c['maxit']=300; c['mmax']=32
    c['continues']=dict(src=str(source.parent),chi=old['chi'],sweeps=len(r['hist']),
                       **{k:old[k] for k in ('E_phys','Sbar1','Sbar6','Sbar12','trunc')})
    c['provenance']=dict(source=str(source),checkpoint=str(checkpoint),verification=v,
                         registration_sha256=digest(HERE/'REGISTRATION_20261116.md'),test_only=test)
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as f:
        json.dump(c,f,indent=2); f.write('\n')
    return c


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True); p.add_argument('--checkpoint',type=Path,required=True)
    p.add_argument('--chi',type=int,default=2560); p.add_argument('--out',type=Path,required=True)
    p.add_argument('--memory-budget-gb',type=float,required=True)
    a=p.parse_args(); r=json.loads(a.source.read_text())
    need=storage(r['cfg']['L'],a.chi)
    if a.memory_budget_gb<need['with_25pct_headroom_GB']:
        p.error('declared available memory is below estimate: '+json.dumps(need))
    c=prepare(a.source,a.checkpoint,a.chi,a.out)
    print(json.dumps(dict(config=str(a.out),tag=c['tag'],resources=need,
        next='Run unmodified codex/work_C195/run195.py with this config and a NEW output directory; no job submitted.'),indent=2))


if __name__=='__main__': main()
