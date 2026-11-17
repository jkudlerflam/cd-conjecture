"""Exact resonant rank-one parents and full-sector formal elimination.

Uses the previously checked literal-bond and Schur infrastructure. No
particle-number truncation and no inference of conformal central charge.
"""
import sys
from pathlib import Path

import numpy as np
from scipy import sparse

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'qubit_cancellation_classification_20261116'))
from algebra import bond,chain,maxabs,occupations,schur_series,sub


def qvalue(m):
    return float(np.tan(np.pi/m)**2)


def parent(r,q):
    """Return forbidden vector, physical bond, telescoping gauge, transfer."""
    psi=np.array([q*r*r,-r,r,1.],complex)
    h=bond(v=1-2*r*r+q*q*r**4,u=2*r*r-2*q*q*r**4,
           t=-r*r,p=q*r*r,a=r+q*r**3,b=-r-q*r**3)
    G=np.array([[0,q*r**3],[q*r**3,0]],complex)
    T=np.array([[-r,1.],[-q*r*r,-r]],complex)
    return psi,h,G,T


def embed(L,h,periodic=True):
    if periodic:return chain(L,h)
    rows=[];cols=[];data=[]
    for state in range(1<<L):
        for j in range(L-1):
            k=j+1;col=2*((state>>j)&1)+((state>>k)&1)
            base=state&~((1<<j)|(1<<k))
            for row in range(4):
                if h[row,col]!=0:
                    rows.append(base|((row>>1)<<j)|((row&1)<<k));cols.append(state);data.append(h[row,col])
    return sparse.coo_matrix((data,(rows,cols)),shape=(1<<L,)*2).tocsr()


def coefficients(L,q,order=8,periodic=True,q_slope=0.):
    """Sum literal |psi(r)><psi(r)|, retaining all open endpoint terms."""
    dim=1<<L
    vectors=[np.array(v,complex) for v in ([0,0,0,1],[0,-1,1,0],[q,0,0,0])]
    if q_slope!=0:vectors.append(np.array([q_slope,0,0,0],complex))
    length=len(vectors);degree=2*length-2
    hs=[sum((np.outer(vectors[i],vectors[n-i].conj()) for i in range(length) if 0<=n-i<length),
            np.zeros((4,4),complex)) for n in range(min(order,degree)+1)]
    H=[embed(L,h,periodic) for h in hs]
    H += [sparse.csr_matrix((dim,dim),dtype=complex) for _ in range(order-degree)]
    return H


def hierarchy(L,q,order=8,stop=6,periodic=True,q_slope=0.):
    """Eliminate each even positive diagonal constraint up to 'stop'.

    Return the entire remaining basis and every remaining coefficient.
    Off-diagonal or odd coefficients are checked before each elimination.
    """
    K=coefficients(L,q,order,periodic,q_slope);words=np.arange(1<<L)
    stages=[];offset=0
    while offset<=stop:
        leading=K[0];diag=leading.diagonal().real
        off=maxabs(leading-sparse.diags(diag))
        assert off<1e-8, (L,q,offset,off)
        assert min(diag)>-1e-8, (L,q,offset,min(diag))
        keep=np.flatnonzero(abs(diag)<1e-8)
        stages.append({'order':offset,'input_dimension':len(words),'kernel_dimension':len(keep),
                       'offdiagonal_norm':off,'minimum_diagonal':float(min(diag))})
        K=schur_series(K,keep);words=words[keep]
        if offset==stop:break
        assert maxabs(K[1])<1e-8,(L,q,offset+1,maxabs(K[1]))
        K=K[2:];offset+=2
    return words,K,stages
