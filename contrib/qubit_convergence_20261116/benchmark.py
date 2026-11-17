#!/usr/bin/env python3
"""Bounded real-gauge C195 benchmark, never a production/convergence result.

Build a full-rank random MPS and host environments, then update three central
bonds with 100 Davidson matvecs, C195 parity blocks and preconditioner. CPU
smoke mode and actual CUDA mode use the same update routines.
"""
import argparse
import gc
import json
import os
from pathlib import Path
import resource
import sys
import time
import numpy as np
from threadpoolctl import threadpool_limits
from continuation import ROOT, storage
sys.path.insert(0,str(ROOT/'codex/work_C195'))
import z2dmrg_xp as Z
import models as Mo


def run(L,chi,xp,nmv=100):
    Z.USE_GPU=xp is not np; Z.SVD_HOST=False
    def host(a): return a.copy() if xp is np else xp.asnumpy(a)
    def sync():
        if xp is not np: xp.cuda.Device().synchronize()
    def now(): sync(); return time.perf_counter()
    t0=now(); j0=L//2-2
    mpo,_,_,_=Mo.candidate('sym40',L,-19.976047250613814,'stag')
    W=[xp.asarray(w) for w in mpo.W]
    D=[1]+[min(chi,2**min(j,L-j)) for j in range(1,L)]+[1]
    n0=[1]+[(d+1)//2 for d in D[1:L]]+[1]
    rng=np.random.default_rng(5)
    def tensor(j):
        out=np.zeros((D[j],2,D[j+1]))
        for qa in (0,1):
            for s in (0,1):
                sl=(Z.bsl(n0[j],D[j],qa),s,Z.bsl(n0[j+1],D[j+1],qa^s))
                out[sl]=rng.standard_normal(out[sl].shape)
        return xp.asarray(out)
    psi=Z.MPS([None]*L,n0,0)
    Le=[None]*(L+1); Re=[None]*(L+1)
    Le[0]=np.zeros((1,mpo.w,1)); Le[0][0,0,0]=1
    Re[L]=np.zeros((1,mpo.w,1)); Re[L][0,2,0]=1
    carry=tensor(0); E=xp.asarray(Le[0])
    for j in range(j0):
        psi.M[j]=carry; psi.M[j+1]=tensor(j+1); Z._qr_move_right(psi,j)
        psi.M[j+1]/=xp.linalg.norm(psi.M[j+1])
        E=Z.env_left(E,psi.M[j],W[j],psi.n0[j],psi.n0[j+1],mpo)
        Le[j+1]=host(E); psi.M[j]=host(psi.M[j]); carry=psi.M[j+1]
    left=carry; carry=tensor(L-1); E=xp.asarray(Re[L])
    for j in range(L-1,j0+1,-1):
        psi.M[j]=carry; psi.M[j-1]=tensor(j-1); Z._rq_move_left(psi,j)
        psi.M[j-1]/=xp.linalg.norm(psi.M[j-1])
        E=Z.env_right(E,psi.M[j],W[j],psi.n0[j],psi.n0[j+1],mpo)
        Re[j]=host(E); psi.M[j]=host(psi.M[j]); carry=psi.M[j-1]
    psi.M[j0]=left; psi.M[j0+1]=carry
    psi.M[j0]/=xp.linalg.norm(xp.tensordot(left,carry,axes=([2],[0])))
    psi.M[j0]=host(psi.M[j0]); psi.M[j0+1]=host(psi.M[j0+1])
    del E,carry,left
    out=dict(L=L,chi=chi,dtype='float64',backend='numpy' if xp is np else 'cupy',build_seconds=now()-t0,
             measured_host_array_GB=sum(a.nbytes for a in psi.M+Le+Re if a is not None)/1e9,
             estimate=storage(L,chi),bonds=[])
    E=xp.asarray(Le[j0]); A=xp.asarray(psi.M[j0])
    for k in range(4):
        j=j0+k; t=now(); D=psi.dims
        R=xp.asarray(Re[j+2]); B=xp.asarray(psi.M[j+1])
        ts=Z.TwoSite(E,W[j],W[j+1],R,psi.n0[j],D[j],psi.n0[j+2],D[j+2],mpo)
        th=xp.tensordot(A,B,axes=([2],[0])); x0=ts.pack_dense(th); ts.build_kron4()
        setup=now()-t; tsolve=now()
        e,v,mv,res=Z.davidson(ts.matvec,x0,None,res_tol=0.,maxit=(nmv-1 if k else 3),mmax=32,
                             keep=6,sigma=1e-3,prec=lambda r,ev:ts.precond4(r,ev,1e-2))
        solve=now()-tsolve; tsplit=now()
        A,B,k0,disc=Z.split(ts.unpack_dense(v),psi.n0[j],D[j],psi.n0[j+2],D[j+2],chi,-1.,'right')
        psi.n0[j+1]=k0; E=Z.env_left(E,A,W[j],psi.n0[j],psi.n0[j+1],mpo)
        psi.M[j]=host(A); psi.M[j+1]=host(B); Le[j+1]=host(E); Re[j+2]=None
        Z.schmidt_probs_center(A,psi.n0[j],D[j],k0,A.shape[2])
        split=now()-tsplit
        if k: out['bonds'].append(dict(j=j,setup_seconds=setup,solve_seconds=solve,
            split_env_entropy_seconds=split,total_seconds=now()-t,matvecs=mv,residual=float(res)))
        A=B
        del ts,th,x0,v,R
    # ru_maxrss is bytes on macOS, KiB on Linux; includes other live CPU arrays.
    out['peak_host_RSS_GB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)/1e9
    if xp is not np:
        free,total=xp.cuda.runtime.memGetInfo()
        out['device_used_GB']=(total-free)/1e9
        out['device_pool_GB']=xp.get_default_memory_pool().total_bytes()/1e9
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--L',type=int,default=512); p.add_argument('--chi',type=int,default=2560)
    p.add_argument('--numpy',action='store_true'); p.add_argument('--memory-budget-gb',type=float,required=True)
    p.add_argument('--out',type=Path,required=True); a=p.parse_args()
    if a.L<16: p.error('L must be at least 16')
    if storage(a.L,a.chi)['with_25pct_headroom_GB']>a.memory_budget_gb:
        p.error('insufficient declared memory; no arrays allocated')
    xp=np
    if not a.numpy:
        import cupy as xp
    with threadpool_limits(limits=1): r=run(a.L,a.chi,xp)
    a.out.write_text(json.dumps(r,indent=2)+'\n'); print(json.dumps(r,indent=2))


if __name__=='__main__':main()
