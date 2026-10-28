"""Exact finite-g second-order one-particle T via fixed-charge resolvents.
The inverse is evaluated by a convergent Neumann iteration with a rigorous
exact-arithmetic truncation bound. Floating-point roundoff is not certified.
Only two charge sectors, (0,1) and (1,0), can connect the required sources.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import time, math, itertools
import numpy as np
import scipy.sparse as sp
from structured import tree,visibility

def charge(n,q,D,a,g):
    S=2*n+1;states=[]
    for typ,chg in [(0,0),(1,1),(2,-1)]:
        k=q-chg
        if 0<=k<=n:
            for c in itertools.combinations(range(n),k):
                bits=sum(1<<j for j in c)
                for p in range(S):states.append((typ,bits,p))
    ix={s:i for i,s in enumerate(states)};rows=[];cols=[];vals=[]
    diag=np.array([a+(D if p>=n else 0) for _,_,p in states])
    edges=[(j,j+1) for j in range(n-1)]+[(n+j,n+j+1) for j in range(n)]+[(0,n)]
    adj=[[] for _ in range(S)]
    for i,j in edges:adj[i].append(j);adj[j].append(i)
    for k,(typ,bits,p) in enumerate(states):
        for j in adj[p]:rows.append(ix[typ,bits,j]);cols.append(k);vals.append(-1.)
        if p<n:
            bit=(bits>>p)&1
            nt=(1 if bit else 2) if typ==0 else (0 if (typ==1 and not bit) or (typ==2 and bit) else None)
            if nt is not None:
                rows.append(ix[nt,bits^(1<<p),p]);cols.append(k)
                vals.append(g*(-1)**((bits&((1<<p)-1)).bit_count()))
    off=sp.csr_matrix((vals,(rows,cols)),shape=(len(states),len(states)))
    return diag,off,ix

def exact_T(n,D=10.,a=4.,g=.05,rtol=1e-12,max_seconds=160.):
    start=time.monotonic();d0,H0,ix0=charge(n,0,D,a,g);d1,H1,ix1=charge(n,1,D,a,g)
    den=d0[:,None]+d1[None,:];inv=1/den
    F,_=tree(n,D,a);target=rtol*2*g*g*np.min(F)
    rho=(2+abs(g))/a
    if rho>=1:raise ValueError('Neumann bound requires (2+g)/a<1')
    iters=max(1,math.ceil(math.log(target*(2*a)*(1-rho))/math.log(rho)))
    T=np.zeros((n,n));errs=[]
    target1=[(ix0[2,1<<i,2*n],ix1[1,0,0]) for i in range(n)]
    target2=[(ix0[0,0,0],ix1[0,1<<i,2*n]) for i in range(n)]
    for side in [0,1]:
        for j in range(n):
            src=(ix0[0,0,2*n],ix1[0,1<<j,0]) if side==0 else (ix0[2,1<<j,0],ix1[1,0,2*n])
            X=np.zeros_like(den)
            for k in range(iters):
                Y=-(H0@X+(H1@X.T).T)*inv
                Y[src]+=inv[src]
                X=Y
                if k%30==0 and time.monotonic()-start>max_seconds:
                    raise TimeoutError((n,side,j,k))
            if side==0:T[:,j]=-np.array([X[i,k] for i,k in target1])
            else:T[:,j]-=[X[i,k] for i,k in target2]
    return visibility(n)*T,dict(n=n,D=D,a=a,g=g,sector_dim=int(den.size),iterations=iters,
        inverse_entry_truncation_bound=rho**iters/(2*a*(1-rho)),seconds=time.monotonic()-start)

if __name__=='__main__':
    import json,sys
    ns=list(map(int,sys.argv[1:])) or [1,2,3,4,6,8]
    records=[]
    for n in ns:
        T,r=exact_T(n)
        F,_=tree(n);ref=-2*visibility(n)*.05**2*F
        r.update(T_scale=float(abs(T).max()),leading_relative_error=float(np.linalg.norm(T-ref)/np.linalg.norm(T)),
                 herm_min_over_scale=float(np.linalg.eigvalsh(-(T+T.T)/2)[0]/abs(T).max()))
        np.savez(f'finite_g_n{n}.npz',T=T)
        records.append(r);print(json.dumps(r),flush=True)
        with open('finite_g_results.json','w') as f:json.dump(records,f,indent=2)
