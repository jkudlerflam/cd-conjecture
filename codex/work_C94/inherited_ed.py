"""Small vacancy ED, all outputs local. NumPy and SciPy only.

Vacuum-clock truncation is exact at lambda=0 within the zero-auxiliary
block at K=infinity. Scalar diagnostics are separately identified.
"""
from pathlib import Path
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import itertools
import json
import math
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

HERE = Path(__file__).resolve().parent

def compositions(total, length):
    for separators in itertools.combinations(range(total+length-1), length-1):
        p = (-1,) + separators + (total+length-1,)
        yield tuple(p[j+1]-p[j]-1 for j in range(length))

def selector(b):
    def a(l): return 2*(1-math.cos(math.pi/(l+1)))
    def d(l): return (l+1)*a(l)-l*a(l+1)
    return 2-(d(b)+d(b-1))/2

def cell_energy(l, mu):
    return mu-2*math.cos(math.pi/(l+1))

def bath(N, V):
    gaps = list(compositions(V, N))
    lookup = {g:i for i,g in enumerate(gaps)}
    rows, cols = [], []
    for i,g in enumerate(gaps):
        for j in range(N):
            if g[j]:
                for dest in ((j-1)%N,(j+1)%N):
                    h = list(g); h[j]-=1; h[dest]+=1
                    rows.append(lookup[tuple(h)]); cols.append(i)
    dim = len(gaps)
    H = sp.csr_matrix((-np.ones(len(rows)),(rows,cols)),shape=(dim,dim))
    G = np.asarray(gaps)
    contacts = (G == 0).astype(float)
    L=N+V
    x=np.zeros((dim,N)); x[:,1:]=np.cumsum(G[:,:-1]+1,axis=1)
    logpsi=np.zeros(dim)
    for i in range(N):
        for j in range(i+1,N):
            logpsi+=np.log(2*np.sin(np.pi*(x[:,j]-x[:,i])/L))
    psi=np.exp(logpsi-logpsi.max()); psi/=np.linalg.norm(psi)
    E0=-2*math.sin(math.pi*N/L)/math.sin(math.pi/L)
    alpha=N/L-math.sin(math.pi*N/L)**2/(N*L*math.sin(math.pi/L)**2)
    res=float(np.linalg.norm(H@psi-E0*psi))
    assert res < 1e-10
    assert np.max(abs((psi*psi)@contacts-alpha)) < 1e-11
    return H,contacts,psi,E0,alpha,res

def ground(H, v0=None):
    e,v=sla.eigsh(H,k=1,which='SA',v0=v0,tol=2e-12,maxiter=15000,ncv=40)
    e=float(e[0]); v=v[:,0]
    res=float(np.linalg.norm(H@v-e*v))
    assert res < 2e-9, res
    return e,v,res

def clock_model(lengths, pins, V, J, mu, cache):
    N=sum(lengths)
    if (N,V) not in cache: cache[N,V]=bath(N,V)
    Hgap,A,psi,E0,alpha,free_res=cache[N,V]
    dims=[1 if p else b for b,p in zip(lengths,pins)]
    states=list(itertools.product(*(range(d) for d in dims)))
    idx={x:i for i,x in enumerate(states)}
    dc=len(states); dg=len(psi)
    markers=np.cumsum(lengths)-1
    diag=mu/2*sum((A[:,(m-1)%N]+A[:,m]) for m in markers)
    H=sp.kron(Hgap+sp.diags(J*diag),sp.eye(dc),format='csr')
    start=0
    for c,(b,pin) in enumerate(zip(lengths,pins)):
        if not pin:
            for q in range(b-1):
                rows=[]; cols=[]
                for s,i in idx.items():
                    if s[c]==q:
                        z=list(s); z[c]+=1; k=idx[tuple(z)]
                        rows.extend([i,k]); cols.extend([k,i])
                hop=sp.csr_matrix((-np.ones(len(rows)),(rows,cols)),shape=(dc,dc))
                H+=J*sp.kron(sp.diags(A[:,start+q]),hop,format='csr')
        start+=b
    clockpsi=np.ones(1)
    for d in dims:
        w=np.sin(np.pi*np.arange(1,d+1)/(d+1)); w/=np.linalg.norm(w)
        clockpsi=np.kron(clockpsi,w)
    e,v,res=ground(H,np.kron(psi,clockpsi))
    charge_weight=np.sum(v.reshape(dg,dc)**2,axis=1)
    pin_contacts=[float(charge_weight@A[:,(m-1)%N])
                  for m,p in zip(markers,pins) if p]
    frozen=sum(mu if p else cell_energy(b,mu) for b,p in zip(lengths,pins))
    return dict(lengths=lengths,pins=pins,N=N,V=V,L=N+V,J=J,t=1.0,
                density=N/(N+V),dimension=H.shape[0],energy=e,
                free_energy=E0,alpha=alpha,frozen_energy=frozen,
                dressing_residual=e-E0-J*alpha*frozen,
                pin_terminal_contacts=pin_contacts,
                eigenpair_residual=res,free_state_residual=free_res)

def run():
    cache={}; clock=[]; scalar=[]
    plans=[(3,[3,3],6,[([3,3],[True,False]),([3,3],[True,True]),([2,4],[False,False])],[.04,.02,.01]),
           (15,[15],2,[([15],[True])],[.04,.02,.01]),
           (15,[15,15],1,[([15,15],[True,False]),([14,16],[False,False]),([10,20],[False,False]),([30],[False])],[.02,.01,.005])]
    for b,good,V,defects,Js in plans:
        mu=selector(b)
        for J in Js:
            ref=clock_model(good,[False]*len(good),V,J,mu,cache)
            for lengths,pins in defects:
                r=clock_model(lengths,pins,V,J,mu,cache)
                leading=J*r['alpha']*(r['frozen_energy']-ref['frozen_energy'])
                row=dict(tag='numerical_exact_vacuum_clock_block',b=b,mu_sel=mu,
                    reference=ref,defect=r,cost=r['energy']-ref['energy'],
                    leading_cost=leading,ratio=(r['energy']-ref['energy'])/leading,
                    residual_difference=r['dressing_residual']-ref['dressing_residual'])
                clock.append(row)
                if lengths==good and any(pins):
                    row['pin_variational_lower_bound']=J*sum(r['pin_terminal_contacts'])
                    assert row['cost']+1e-9 >= row['pin_variational_lower_bound']
                print('clock',b,V,J,lengths,pins,row['cost'],row['ratio'],flush=True)
                (HERE/'ed_results.json').write_text(json.dumps(dict(clock=clock,scalar=scalar),indent=2)+'\n')
    mu=selector(15); eb=cell_energy(15,mu)/15
    dp=2*math.cos(math.pi/16)
    dl={l:cell_energy(l,mu)-l*eb for l in [14,16,30,60]}
    for N in [6,8]:
        H,A,psi,E0,alpha,res=bath(N,N)
        for J in [.04,.02,.01]:
            refH=H+sp.diags(J*eb*A.sum(axis=1))
            er,vr,rr=ground(refH,psi)
            for name,profile in [('pin',dp*A[:,0]),('pin_pair',dp*(A[:,0]+A[:,1])),
                                 ('wrong_14_16',dl[14]*A[:,0]+dl[16]*A[:,1]),
                                 ('wrong_60',dl[60]*A[:,0])]:
                hd=refH+sp.diags(J*profile)
                e,v,r=ground(hd,vr)
                lead=J*float((psi*psi)@profile)
                scalar.append(dict(tag='numerical_scalar_contact_surrogate',N=N,L=2*N,
                    dimension=H.shape[0],J=J,t=1,name=name,alpha=alpha,cost=e-er,
                    leading_cost=lead,ratio=(e-er)/lead,eigenpair_residual=r,
                    reference_residual=rr,relaxation_over_J2=(lead-(e-er))/J**2))
                print('scalar',N,J,name,e-er,(e-er)/lead,flush=True)
    result=dict(clock=clock,scalar=scalar)
    (HERE/'ed_results.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__': run()
