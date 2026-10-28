#!/usr/bin/env python3
"""Independent vacancy audit; no campaign implementations imported.

Only writes results beside this file. Requires numpy/scipy.
"""
import itertools as it
import json
import math
from pathlib import Path
import time
import numpy as np
from scipy.linalg import null_space
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh


def squeezing(L, N, q):
    pos = list(it.combinations(range(L), N))
    labs = list(it.product(range(q), repeat=N))
    pi = {p: i for i, p in enumerate(pos)}
    si = {s: i for i, s in enumerate(labs)}
    a, d = len(pos), len(labs)
    rng = np.random.default_rng(280026 + L + N + q)
    z = rng.normal(size=(q*q,q*q))+1j*rng.normal(size=(q*q,q*q))
    h = (z+z.conj().T)/2
    words=[]
    for p in pos:
        for s in labs:
            w=[-1]*L
            for x, lab in zip(p,s): w[x]=lab
            words.append(tuple(w))
    wi={w:i for i,w in enumerate(words)}
    direct_t=np.zeros((a*d,a*d),complex)
    direct_v=np.zeros_like(direct_t)
    for col,w in enumerate(words):
        for x in range(L):
            y=(x+1)%L
            if (w[x]<0)!=(w[y]<0):
                v=list(w);v[x],v[y]=v[y],v[x]
                direct_t[wi[tuple(v)],col]-=1
            elif w[x]>=0 and w[y]>=0:
                for u,v in it.product(range(q),repeat=2):
                    out=list(w);out[x],out[y]=u,v
                    direct_v[wi[tuple(out)],col]+=h[u*q+v,w[x]*q+w[y]]
    R=np.zeros((d,d),complex)
    hs=[]
    for j in range(N):
        hj=np.zeros((d,d),complex)
        k=(j+1)%N
        for col,s in enumerate(labs):
            for u,v in it.product(range(q),repeat=2):
                out=list(s);out[j],out[k]=u,v
                hj[si[tuple(out)],col]+=h[u*q+v,s[j]*q+s[k]]
        hs.append(hj)
    for col,s in enumerate(labs):R[si[(s[-1],)+s[:-1]],col]=1
    bulk=np.zeros((a,a),complex);cut=np.zeros_like(bulk)
    for col,p in enumerate(pos):
        for x in p:
            for dx in (-1,1):
                y=(x+dx)%L
                if y in p:continue
                out=tuple(sorted((set(p)-{x})|{y}))
                row=pi[out]
                if x==L-1 and y==0:cut[row,col]-=1
                elif x==0 and y==L-1:pass
                else:bulk[row,col]-=1
    squeezed_t=np.kron(bulk,np.eye(d))+np.kron(cut,R)+np.kron(cut.conj().T,R.conj().T)
    squeezed_v=np.zeros_like(direct_t)
    for j in range(N):
        active=np.array([(p[(j+1)%N]-p[j])%L==1 for p in pos],float)
        squeezed_v+=np.kron(np.diag(active),hs[j])
    W=[];half_flux_branches=0
    for r in range(N):
        phase=np.exp(2j*np.pi*r/N)
        internal=null_space(R-phase*np.eye(d))
        if not internal.shape[1]:continue
        hc=bulk+phase*cut+phase.conjugate()*cut.conj().T
        ev,U=np.linalg.eigh(hc)
        sea=U[:,ev<ev[0]+1e-10]
        if 2*r==N:half_flux_branches=sea.shape[1]
        W.extend(np.kron(sea[:,k:k+1],internal) for k in range(sea.shape[1]))
    W=np.column_stack(W)
    alpha=N/L-math.sin(math.pi*N/L)**2/(N*L*math.sin(math.pi/L)**2)
    projected=W.conj().T@direct_v@W
    expected=alpha*W.conj().T@np.kron(np.eye(a),sum(hs))@W
    result=dict(L=L,N=N,q=q,dimension=a*d,kinetic_error=float(np.max(abs(direct_t-squeezed_t))),interaction_error=float(np.max(abs(direct_v-squeezed_v))),projection_error=float(np.max(abs(projected-expected))),half_flux_branches=half_flux_branches)
    assert max(result[k] for k in ('kinetic_error','interaction_error','projection_error'))<1e-10
    return result


def compositions(total, length):
    if length==1:
        yield (total,)
    else:
        for a in range(total+1):
            for tail in compositions(total-a,length-1):yield (a,)+tail


def phase_ground(N,H,flux=0):
    states=list(compositions(H,N)); lookup={s:i for i,s in enumerate(states)}
    rows=[];cols=[];vals=[]
    for col,s in enumerate(states):
        for j in range(N):
            if not s[j]:continue
            for dx in (-1,1):
                k=(j+dx)%N
                out=list(s);out[j]-=1;out[k]+=1
                rows.append(lookup[tuple(out)]);cols.append(col);vals.append(-np.exp(1j*dx*flux/N))
    mat=coo_matrix((vals,(rows,cols)),shape=(len(states),len(states))).tocsr()
    return float(eigsh(mat,k=1,which='SA',tol=1e-12,return_eigenvectors=False)[0])


def stiffness(N):
    step=0.01
    e=phase_ground(N,N)
    ef=phase_ground(N,N,step)
    em=phase_ground(N,N-1);ep=phase_ground(N,N+1)
    drude=N*2*(ef-e)/step**2
    inverse_compressibility=N*(em+ep-2*e)
    K=math.pi*math.sqrt(drude/inverse_compressibility)
    v=math.sqrt(drude*inverse_compressibility)
    exact=-2/math.sin(math.pi/(2*N))
    assert abs(e-exact)<1e-9
    return dict(N=N,H=N,dimension=math.comb(2*N-1,N-1),energy_error=e-exact,drude=drude,inverse_compressibility=inverse_compressibility,K=K,v=v)


def exhaustive_threshold(L):
    N=L//2;worst=float('inf');tested=0
    for r in range(N):
        theta=2*np.pi*r/N
        ks=(2*np.pi*np.arange(L)-theta-np.pi*(N-1))/L
        eps=-2*np.cos(ks)
        base=tuple(np.argsort(eps)[:N]);e0=sum(eps[list(base)]);p0=sum(ks[list(base)])
        for occupation in it.combinations(range(L),N):
            de=sum(eps[list(occupation)])-e0
            dp=sum(ks[list(occupation)])-p0
            distance=abs((dp+np.pi/2)%np.pi-np.pi/2)
            slack=de-(4/np.pi)*distance
            worst=min(worst,slack);tested+=1
    assert worst>-1e-10
    return dict(L=L,N=N,all_states_and_twists=tested,minimum_bound_slack=worst)


if __name__=='__main__':
    start=time.monotonic()
    results=dict(squeezing=[squeezing(*args) for args in [(4,2,2),(5,2,3),(5,3,2),(6,4,2),(7,3,2)]],stiffness=[stiffness(N) for N in (3,4,5,6,7)],threshold=[exhaustive_threshold(L) for L in (6,8,10,12,14)])
    results['elapsed_seconds']=time.monotonic()-start
    results['status']='passed'
    Path(__file__).with_name('independent_results.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))
