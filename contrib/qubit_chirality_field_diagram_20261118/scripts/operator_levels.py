"""Momentum-projected operator evidence; no Ising dimension imposed in solving."""
import argparse
import json
import signal
import time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from ed import Chain, HERE
from threadpoolctl import threadpool_limits


def momentum_basis(c, k):
    visited=np.zeros(c.n,dtype=bool); rows=[]; cols=[]; vals=[]; count=0
    for i in range(c.n):
        if visited[i]: continue
        orbit=[]; j=i
        while not visited[j]:
            orbit.append(j);visited[j]=True;j=c.trans_dest[j]
        n=len(orbit)
        if k*n % c.L: continue
        for r,j in enumerate(orbit):
            rows.append(j);cols.append(count)
            vals.append(np.exp(-2j*np.pi*k*r/c.L)/np.sqrt(n)/c.phase[j])
        count+=1
    return coo_matrix((vals,(rows,cols)),shape=(c.n,count)).tocsr()


def block(c,H,k,nev=6):
    Q=momentum_basis(c,k); Hq=(Q.conj().T@H@Q).tocsr()
    anti=Hq-Hq.conj().T
    assert max(abs(anti.data),default=0)<1e-11
    vals,V=eigsh(Hq,k=nev,which='SA',tol=1e-12,maxiter=10000,
                 v0=np.random.default_rng(321+k).normal(size=Hq.shape[0]))
    ix=np.argsort(vals);vals,V=vals[ix],V[:,ix]
    residual=max(np.linalg.norm(Hq@V[:,i]-vals[i]*V[:,i]) for i in range(nev))
    assert residual<1e-8
    return vals,V,Q,float(residual)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--sizes',default='12,16')
    ap.add_argument('--fields',default='1.738,1.7466,1.755')
    ap.add_argument('--output',default='operator_levels.json')
    args=ap.parse_args()
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('90-second cap')))
    signal.alarm(90);start=time.monotonic();out=[]
    records=json.loads((HERE/'results'/'coarse.json').read_text())['records']
    records+=json.loads((HERE/'results'/'larger.json').read_text())['records']
    for record in records:
        L,h=record['L'],record['h']
        if h not in list(map(float,args.fields.split(','))) or L not in list(map(int,args.sizes.split(','))): continue
        c=Chain(L,record['parity']); H=c.H(h); kg=record['k']
        vals,V,Q,res=block(c,H,kg)
        E=vals[0];v=c.phase*(Q@V[:,0]); assert abs(E-record['E'])<1e-8
        transverse={}
        inv=np.full(1<<L,-1,dtype=np.int64)
        opposite=np.array([s for s in range(1<<L) if s.bit_count()%2!=c.parity])
        inv[opposite]=np.arange(len(opposite))
        for name in ('X','Y'):
            W=[]
            for j in range(L):
                w=np.zeros(c.n,dtype=complex)
                w[inv[c.states^(1<<j)]]=v*(1 if name=='X' else 1j*c.z[j])
                W.append(w)
            W=np.array(W); raw=W.conj()@W.T
            corr=np.array([np.mean([raw[j,(j+r)%L] for j in range(L)]) for r in range(L)])
            transverse[name]=dict(connected_re=corr.real.tolist(),connected_im=corr.imag.tolist(),
                                  structure=np.fft.fft(corr).real.tolist())
        blockout=[]
        for shift in (0,1,-1,L//2):
            k=(kg+shift)%L
            es,vs,q,rr=(vals,V,Q,res) if shift==0 else block(c,H,k)
            weights={}
            for name in ('Z','bond','kappa','B'):
                op=sum(np.exp(2j*np.pi*shift*j/L)*c.local(name,j,v) for j in range(L))
                # With T|g>=e^{ikg}|g>, coefficient e^{+iqj} produces kg-q.
                # pi has no sign ambiguity; +/-1 block weights below are
                # computed with the coefficient yielding the chosen block.
                if shift not in (0,L//2):
                    op=sum(np.exp(-2j*np.pi*shift*j/L)*c.local(name,j,v) for j in range(L))
                if shift==0: op-=v*np.vdot(v,op)
                b=q.conj().T@(op/c.phase)
                norm=float(np.vdot(op,op).real/L)
                projected=float(np.vdot(b,b).real/L)
                assert abs(norm-projected)<1e-8
                wt=abs(vs.conj().T@b)**2/L
                weights[name]=dict(total_weight=norm,low_level_weights=wt.tolist(),
                                   captured_fraction=float(sum(wt)/norm) if norm>1e-16 else None)
            blockout.append(dict(shift=shift,k=k,gaps=(es-E).tolist(),residual=rr,operators=weights))
        out.append(dict(L=L,h=h,parity=c.parity,k_ground=kg,ground_energy=float(E),
                        ground_residual=res,transverse=transverse,blocks=blockout))
        print('operator blocks',L,h,flush=True)
    signal.alarm(0)
    (HERE/'results'/args.output).write_text(json.dumps(dict(records=out,runtime_seconds=time.monotonic()-start),indent=2)+'\n')


if __name__=='__main__':
    with threadpool_limits(limits=1): main()
