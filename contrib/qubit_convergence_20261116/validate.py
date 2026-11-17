#!/usr/bin/env python3
"""Small-chain ED, endpoint, parity, offload, interruption and continuation gates.

Runs only this contribution's ignored runtime directory. These are solver
validation results, not large-chi evidence for a conformal field theory.
"""
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits
from continuation import ROOT,HERE,prepare,verify
from controls import exact
import z2dmrg as Z
import z2dmrg_xp as ZX
import models as Mo

RUN=HERE/'runtime'
HZ=-19.976047250613814


def dense_mpo(mpo):
    terms={0:np.ones((1,1))}
    for W in mpo.W:
        nxt={}
        for a,T in terms.items():
            for b in range(mpo.w):
                if np.any(W[a,b]):
                    K=np.kron(T,W[a,b]); nxt[b]=nxt.get(b,0)+K
        terms=nxt
    return terms[2]+mpo.const*np.eye(2**mpo.L)


def literal(L,hz,half_ends=False):
    X=np.array([[0,1],[1,0]]); Y=np.array([[0,-1j],[1j,0]]); zz=np.diag([1,-1])
    bond=(np.kron(X,X)-np.kron(Y,Y))/2+(np.kron(X,Y)-np.kron(Y,X))/160+10*np.kron(zz,zz)
    H=np.zeros((2**L,2**L),complex)
    for j in range(L-1): H+=np.kron(np.eye(2**j),np.kron(bond,np.eye(2**(L-j-2))))
    for j in range(L):
        fac=.5 if half_ends and j in (0,L-1) else 1.
        H+=fac*hz*np.kron(np.eye(2**j),np.kron(zz,np.eye(2**(L-j-1))))
    return H


def run_driver(cfg,dest,extra=None):
    cp=RUN/(cfg['tag']+'.cfg.json'); cp.write_text(json.dumps(cfg,indent=2)+'\n')
    env=dict(os.environ,C190G_BACKEND='numpy',C195_OFFLOAD='1',C195_OFFLOAD_TEST='1',
             C190G_GPU_MIN_CHI='16',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    env.update(extra or {})
    dest=RUN/dest; dest.mkdir(exist_ok=True)
    with (RUN/(dest.name+'.log')).open('a') as log:
        subprocess.run([sys.executable,str(ROOT/'codex/work_C195/run195.py'),str(cp),str(dest)],
                       env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    return dest/(cfg['tag']+'.json'),dest/'ckpt'/(cfg['tag']+'.npz')


def main():
    start=time.time(); RUN.mkdir(exist_ok=True)
    result=dict(kind='small validation only',numpy=np.__version__,python=sys.version,ed=[],resume=[],local_controls=[])
    ZX.USE_GPU=False; ZX.OFFLOAD=True; ZX.OFFLOAD_TEST=True; ZX.GPU_MIN_CHI=16
    with threadpool_limits(limits=1):
        for L in (8,10):
            H=literal(L,HZ)
            mpo,gph,shift,_=Mo.candidate('sym40',L,HZ,gauge='stag')
            phases=np.ones(1,complex)
            for g in gph:
                th=-np.angle(g); phases=np.kron(phases,[np.exp(-.5j*th),np.exp(.5j*th)])
            Hg=phases[:,None]*H*phases.conj()[None,:]
            dM=float(np.max(abs(dense_mpo(mpo)-L*shift*np.eye(2**L)-Hg)))
            half=literal(L,HZ,True)
            endpoint_error=float(np.linalg.norm(H-half,2))
            assert dM<1e-11 and abs(endpoint_error-abs(HZ))<1e-10
            for P in (0,1):
                ix=[i for i in range(2**L) if i.bit_count()%2==P]
                ee=eigh(H[np.ix_(ix,ix)],subset_by_index=[0,0],eigvals_only=True)[0]
                psi=ZX.random_mps(L,P,D0=8,seed=1,dtype=float)
                spec=[dict(chi=32,alpha=a,res_tol=1e-10,maxit=300,mmax=32) for a in (1e-3,0.,0.,0.)]
                psi,h=ZX.dmrg(psi,mpo,spec)
                e=h[-1]['E']+mpo.const-L*shift
                row=dict(L=L,P=P,ED_E=float(ee),DMRG_E=e,error_E=abs(e-ee),
                         matrix_error=dM,missing_endpoint_half_operator_norm=endpoint_error,
                         res_max=h[-1]['res_max'],trunc=h[-1]['trunc'],sweep_seconds=h[-1]['t'])
                assert row['error_E']<1e-10 and row['res_max']<1e-8
                result['ed'].append(row); print('ED',row,flush=True)
    for P in (0,1):
        cfg=dict(tag=f'local_sym_L16_P{P}',model='sym40',gauge='stag',hz=HZ,L=16,P=P,init='rand:1',
            sched=[[8,1e-3,1e-8,2],[16,0.,1e-9,2]],
            final=dict(chi=32,res_tol=1e-9,min_sweeps=3,max_sweeps=5,tolE=2e-9,tolS=5e-4),
            gwin=8,variance=True,maxit=300)
        full,ck=run_driver(cfg,f'full_P{P}',{'C195_OFFLOAD_TEST':'0'})
        interrupted,ckr=run_driver(cfg,f'resume_P{P}',{'C195_EXIT_AFTER':'3'})
        # An old complete output means rerunning is idempotent, otherwise test the stop.
        stopped=json.loads(interrupted.read_text())
        assert stopped['status']=='done' or len(stopped['hist'])==3
        resumed,_=run_driver(cfg,f'resume_P{P}')
        a=json.loads(full.read_text()); b=json.loads(resumed.read_text())
        row=dict(P=P,dE=abs(a['E_phys']-b['E_phys']),dS=float(np.max(abs(np.array(a['S'])-b['S']))),
                 ddens=float(np.max(abs(np.array(a['meas']['dens'])-b['meas']['dens']))),
                 verified_checkpoint=verify(resumed,ckr)[1])
        assert row['dE']<1e-10 and row['dS']<1e-8
        target=RUN/f'continue_P{P}.json'
        if target.exists(): target.unlink() # only this local test configuration
        new=prepare(resumed,ckr,40,target,test=True)
        newf,newck=run_driver(new,f'continued_P{P}')
        row['continued_checkpoint']=verify(newf,newck)[1]
        # A wrong parity JSON must be rejected, not just renamed.
        bad=json.loads(resumed.read_text()); bad['cfg']['P']=1-P
        badfile=RUN/f'wrong_parity_P{P}.json'; badfile.write_text(json.dumps(bad))
        try: verify(badfile,ckr)
        except ValueError as e: row['wrong_parity_rejected']=str(e)
        else: raise AssertionError('wrong parity was accepted')
        result['resume'].append(row); print('resume',row,flush=True)
    for legs in (1,2):
        for L in (32,64):
            n=L//legs; N=round(.1957*n)
            mu=float(-np.cos(np.pi*N/(n+1))-np.cos(np.pi*(N+1)/(n+1)))
            cfg=dict(tag=f'local_xx{legs}_L{L}',model=f'xx{legs}',mu=mu,L=L,P=(N*legs)%2,init='rand:1',
                sched=[[16,1e-3,1e-8,2],[32,0.,1e-9,2],[64,0.,1e-9,2]],
                final=dict(chi=128,res_tol=1e-9,min_sweeps=3,max_sweeps=5,tolE=2e-9,tolS=5e-4),
                gwin=8,variance=False,maxit=300)
            f,_=run_driver(cfg,f'control_xx{legs}_L{L}'); r=json.loads(f.read_text())
            with threadpool_limits(limits=1): ex=exact(L,mu,legs)
            row=dict(model=f'xx{legs}',L=L,mu=mu,N_exact=ex['N'],
                error_E=r['hist'][-1]['E_phys']-ex['E'],
                error_Sbar={str(w):r[f'Sbar{w}']-ex['Sbar'][str(w)] for w in (1,6,12)},
                Sbar={str(w):r[f'Sbar{w}'] for w in (1,6,12)},exact_Sbar=ex['Sbar'],
                res_max=r['hist'][-1]['res_max'],trunc=r['hist'][-1]['trunc'],
                sweep_seconds=r['hist'][-1]['t'])
            result['local_controls'].append(row); print('control',row,flush=True)
            # Small finite chi need not hit the production 0.001 gate; retain misses.
            assert abs(row['error_E'])<1e-5 and max(abs(v) for v in row['error_Sbar'].values())<.01
    result['elapsed_seconds']=time.time()-start; result['validation_status']='PASS'
    (HERE/'validation.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__': main()
