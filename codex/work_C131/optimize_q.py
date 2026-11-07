"""C131 finite-ring search. NumPy/SciPy only. Empirical maxima, not certificates.
Coordinates: Jx,Jy,Jz,Dx,Dy,Dz,bx,by,bz with S=sigma/2,
h=sum J_a S_a S_a + D.(S cross S) + b.(S1+1S)/2.
Eight angular coordinates cover the coefficient sphere; Q is scale invariant.
"""
import argparse, json, time, os
from pathlib import Path
from functools import lru_cache
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from scipy.optimize import minimize
S = np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)/2
I = np.eye(2)
B = np.array([np.kron(s,s) for s in S] +
 [np.kron(S[a],S[b])-np.kron(S[b],S[a]) for a,b in [(1,2),(2,0),(0,1)]] +
 [(np.kron(s,I)+np.kron(I,s))/2 for s in S])
G = np.array([np.kron(s,I)-np.kron(I,s) for s in S])

def gauge_M(h, thorough=False):
    def f(o):
        e=np.linalg.eigvalsh(h+np.einsum('a,aij->ij',o,G))
        return (e[-1]-e[0])/2
    starts=[np.zeros(3)]
    if thorough: starts += list(np.random.default_rng(177).normal(size=(5,3)))
    rs=[minimize(f,o,method='Nelder-Mead',options=dict(xatol=2e-9,fatol=2e-11,maxiter=1200)) for o in starts]
    r=min(rs,key=lambda r:r.fun)
    return float(r.fun),r.x

def gauge_bracket(h):
    upper,o=gauge_M(h,True)
    # Dual: tr X=tr X G_a=0, ||X||_1<=1 implies tr(Xh)<=M.
    # Use all 12 gauge-orthogonal Pauli directions and optimize the convex ball.
    p=np.r_[B, np.array([np.kron(S[a],S[b])+np.kron(S[b],S[a]) for a,b in [(1,2),(2,0),(0,1)]])]
    norms=np.einsum('aij,aji->a',p,p).real
    hh=h+np.einsum('a,aij->ij',o,G)
    e,v=np.linalg.eigh(hh)
    high=v[:,e>e[-1]-1e-6]; low=v[:,e<e[0]+1e-6]
    X=(high@high.conj().T/high.shape[1]-low@low.conj().T/low.shape[1])/2
    x=np.einsum('aij,ji->a',p,X).real/norms
    coeff=np.einsum('aij,ji->a',p,h).real
    def nuc(x): return float(np.abs(np.linalg.eigvalsh(np.einsum('a,aij->ij',x,p))).sum())
    best=[float('-inf'),x.copy()]
    def retain(z):
        zz=z/max(1,nuc(z)); value=float(coeff@zz)
        if value>best[0]: best[:]=[value,zz.copy()]
    retain(x)
    r=minimize(lambda x:-coeff@x,x,jac=lambda x:-coeff,method='SLSQP',constraints=[dict(type='ineq',fun=lambda x:1-nuc(x))],callback=retain,options=dict(ftol=1e-12,maxiter=1500))
    retain(r.x)
    lower=best[0]
    return dict(lower=lower,upper=upper,gap=upper-lower,o=o.tolist(),dual_success=bool(r.success))

@lru_cache(maxsize=4)
def ring_basis(N):
    dim=1<<N; states=np.arange(dim,dtype=np.int64); matrices=[]
    for h in B:
        rows=[]; cols=[]; vals=[]
        for j in range(N):
            k=(j+1)%N
            old=2*((states>>j)&1)+((states>>k)&1)
            base=states & ~((1<<j)|(1<<k))
            for out in range(4):
                val=h[out,old]; keep=np.flatnonzero(val)
                rows.append(base[keep]|((out>>1)<<j)|((out&1)<<k))
                cols.append(states[keep]); vals.append(val[keep])
        A=sp.coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(dim,dim)).tocsr()
        A.eliminate_zeros(); matrices.append(A)
    return matrices

def ring(c,N):
    A=ring_basis(N)
    H=sum((v*a for v,a in zip(c,A) if v!=0),sp.csr_matrix(A[0].shape,dtype=complex))
    if not np.any(H.data.imag): H=H.real
    return H

def energy(c,N,tol=2e-10,verify=False):
    H=ring(c,N)
    if N<=6:
        e,v=np.linalg.eigh(H.toarray()); E=float(e[0]); vec=v[:,0]
    else:
        rng=np.random.default_rng(941+N)
        v0=rng.normal(size=1<<N)
        if np.iscomplexobj(H.data): v0=v0+1j*rng.normal(size=1<<N)
        e,v=eigsh(H,k=1,which='SA',tol=tol,v0=v0,ncv=30,maxiter=10000)
        E=float(e[0]); vec=v[:,0]
    residual=float(np.linalg.norm(H@vec-E*vec)) if verify else None
    return (E,residual) if verify else E

def to_angles(c):
    c=np.asarray(c)/np.linalg.norm(c)
    return np.r_[[np.arctan2(np.linalg.norm(c[i+1:]),c[i]) for i in range(7)],np.arctan2(c[8],c[7])]

def from_angles(a):
    c=np.empty(9); t=1.
    for i in range(8): c[i]=t*np.cos(a[i]); t*=np.sin(a[i])
    c[8]=t
    return c

def evaluate(c,L,verify=False):
    h=np.einsum('a,aij->ij',c,B)
    if verify:
        m=gauge_bracket(h); e,r=energy(c,L,1e-12,True); ee,rr=energy(c,2*L,1e-12,True)
        return dict(L=L,parity='even' if L%2==0 else 'odd',Q=L*(ee-2*e)/m['upper'],M=m,E_L=e,E_2L=ee,residual_L=r,residual_2L=rr,coefficients=(c/m['upper']).tolist())
    m,_=gauge_M(h)
    return L*(energy(c,2*L)-2*energy(c,L))/m

def seed_bond(seed):
    rng=np.random.default_rng(seed+13100)
    if seed==0: c=np.r_[np.ones(3),np.zeros(6)]
    elif seed==1: c=np.r_[[1,1,0],np.zeros(6)]
    elif seed==2: c=np.r_[[1,-1,-1],np.zeros(6)]
    elif seed==3: c=np.r_[[1,1,-.95],np.zeros(6)]
    elif seed==4: c=np.r_[[1.5,.5,0],np.zeros(5),1]
    elif seed%3==0: c=np.r_[np.ones(3),np.zeros(6)]+.2*rng.normal(size=9)
    else: c=rng.normal(size=9)
    return c/np.linalg.norm(c)

def run(L,seed,budget,out):
    started=time.time(); c=seed_bond(seed)
    base=np.ones(9)/3; v=base-c
    rotation=np.eye(9) if np.linalg.norm(v)<1e-12 else np.eye(9)-2*np.outer(v,v)/(v@v)
    a=to_angles(base); best=[-np.inf,c.copy()]; count=0
    def loss(a):
        nonlocal count
        c=rotation@from_angles(a); q=evaluate(c,L); count+=1
        if q>best[0]: best[:]=[q,c.copy()]
        return -q
    # Nonzero simplex edges avoid the tiny default edges at vanishing angles.
    simplex=np.vstack([a,a+0.16*np.eye(8)])
    r=minimize(loss,a,method='Nelder-Mead',options=dict(maxfev=budget,xatol=2e-6,fatol=2e-8,initial_simplex=simplex))
    result=evaluate(best[1],L,True)
    result.update(seed=seed,nfev=count,success=bool(r.success),message=str(r.message),seconds=time.time()-started,method='Nelder-Mead, eight angles',numpy=np.__version__)
    Path(out).write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--L',type=int,default=2); ap.add_argument('--seed',type=int,default=0); ap.add_argument('--budget',type=int,default=1000); ap.add_argument('--out',default='search.json'); ap.add_argument('--task',type=int)
    args=ap.parse_args()
    if args.task is not None:
        args.L=[2,4,6,8,3,5,7][args.task//32]; args.seed=args.task%32
        args.out=f'results/L{args.L}_seed{args.seed:02}.json'; Path('results').mkdir(exist_ok=True)
    run(args.L,args.seed,args.budget,args.out)
