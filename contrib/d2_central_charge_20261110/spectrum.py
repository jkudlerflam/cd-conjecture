#!/usr/bin/env python3
"""Translation-orbit ED of the literal Pauli ring, no number projection."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
HERE=Path(__file__).resolve().parent

def orbits(L,P):
    states=np.arange(1<<L,dtype=np.int32);par=np.zeros(len(states),np.int8)
    for j in range(L):par^=((states>>j)&1).astype(np.int8)
    members=states[par==P];smallest=members.copy();rot=members.copy()
    for _ in range(1,L):
        rot=(rot>>1)|((rot&1)<<(L-1));np.minimum(smallest,rot,out=smallest)
    reps=np.unique(smallest);ids=np.full(1<<L,-1,np.int32);shifts=np.full(1<<L,-1,np.int8)
    periods=np.zeros(len(reps),np.int8)
    for i,r in enumerate(reps):
        s=int(r);t=0
        while True:
            ids[s]=i;shifts[s]=t;t+=1;s=(s>>1)|((s&1)<<(L-1))
            if s==r:break
        periods[i]=t
    return reps,ids,shifts,periods

def build(L,k,p,data):
    reps,ids,shifts,periods=data
    allowed=(periods.astype(int)*k)%L==0;keep=np.flatnonzero(allowed)
    idx=np.full(len(reps),-1,np.int32);idx[keep]=np.arange(len(keep))
    states=reps[keep];n=periods[keep];N=len(states);source=np.arange(N)
    z=np.array([1-2*((states>>(L-1-j))&1) for j in range(L)])
    diag=p['hz']*z.sum(axis=0).astype(float);rr=[];cc=[];dd=[]
    for j in range(L):
        l=(j+1)%L;diag+=p['Jz']*z[j]*z[l]
        amp=p['Jx']-p['Jy']*z[j]*z[l]+1j*p['Dz']*(z[l]-z[j])
        dest=states^(1<<(L-1-j))^(1<<(L-1-l));target=ids[dest]
        good=allowed[target]&(amp!=0)
        val=amp[good]*np.sqrt(n[good]/periods[target[good]])*np.exp(2j*np.pi*k/L*shifts[dest[good]])
        rr.append(idx[target[good]]);cc.append(source[good]);dd.append(val)
    rr.append(source);cc.append(source);dd.append(diag)
    H=coo_matrix((np.concatenate(dd),(np.concatenate(rr),np.concatenate(cc))),shape=(N,N)).tocsr()
    skew=H-H.getH();err=float(abs(skew).max())
    assert err<1e-10,(L,k,err)
    return H,err

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--model',choices=['candidate','xx'],default='candidate')
    ap.add_argument('--sizes',type=int,nargs='+',default=[8,12,16,20,24]);a=ap.parse_args()
    p=dict(Jx=.5,Jy=-.5,Jz=.625,Dz=.1,hz=-.875) if a.model=='candidate' else dict(Jx=.5,Jy=.5,Jz=0.,Dz=0.,hz=0.)
    out=HERE/f'spectrum_{a.model}.json';rows=json.loads(out.read_text())['rows'] if out.exists() else []
    for L in a.sizes:
        if any(r['L']==L for r in rows):continue
        start=time.time();data=orbits(L,0);sectors=[]
        for k in [0,1,L-1]:
            H,defect=build(L,k,p,data)
            ev,v=eigsh(H,k=3,which='SA',tol=1e-12,maxiter=100000,v0=np.random.default_rng(42).normal(size=H.shape[0]))
            order=np.argsort(ev);ev=ev[order];v=v[:,order]
            residual=[float(np.linalg.norm(H@v[:,i]-ev[i]*v[:,i])) for i in range(len(ev))]
            sectors.append(dict(k=k,dimension=H.shape[0],energies=ev.tolist(),residuals=residual,hermiticity_defect=defect))
        E0=sectors[0]['energies'][0]
        row=dict(L=L,parity=0,sectors=sectors,E0=E0,
            velocities=[L*(r['energies'][0]-E0)/(2*np.pi) for r in sectors[1:]],elapsed_seconds=time.time()-start)
        rows.append(row);tmp=out.with_suffix('.tmp');tmp.write_text(json.dumps(dict(model=a.model,couplings=p,rows=rows),indent=2)+'\n');tmp.replace(out)
        print(json.dumps(row),flush=True)
if __name__=='__main__':main()
