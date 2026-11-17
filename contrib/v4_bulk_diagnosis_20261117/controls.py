#!/usr/bin/env python3
"""Exact critical-density and ordered-state calibrations of the same metrics."""
import numpy as np
from threadpoolctl import threadpool_limits
from phase_common import HERE,write_json
from observables import spectral_metrics


def free_correlation(L,legs):
    m=L//legs;N=int(np.floor(.296*(m+1)));q=np.arange(1,N+1)
    U=np.sqrt(2/(m+1))*np.sin(np.pi*np.arange(1,m+1)[:,None]*q/(m+1));leg=U@U.T
    C=np.zeros((L,L))
    for k in range(legs):C[k::legs,k::legs]=leg
    return C,N*legs


def main():
    rows=[]
    with threadpool_limits(limits=1):
        # Independent many-body Slater-vector check of the Wick expressions.
        from phase_common import Z
        L=8;N=round(.296*L);q=np.arange(1,N+1)
        U=np.sqrt(2/(L+1))*np.sin(np.pi*np.arange(1,L+1)[:,None]*q/(L+1));C=U@U.T;Q=np.eye(L)-C
        v=np.zeros(2**L);bits=np.array([[(s>>(L-1-j))&1 for s in range(2**L)] for j in range(L)])
        for s in range(2**L):
            ids=np.flatnonzero(bits[:,s])
            if len(ids)==N:v[s]=np.linalg.det(U[ids])
        assert abs(v@v-1)<1e-12
        density_conn=(bits*v**2)@bits.T-np.outer(bits@(v**2),bits@(v**2))
        density_error=float(np.max(abs(density_conn-(np.diag(np.diag(C))-C**2))))
        hop=-np.kron(Z.SP,Z.SM)-np.kron(Z.SM,Z.SP)
        vv=np.array([np.kron(np.eye(2**j),np.kron(hop,np.eye(2**(L-j-2))))@v for j in range(L-1)])
        means=vv.conj()@v;bc=vv.conj()@vv.T-np.outer(means.conj(),means)
        wick=np.array([[sum(C[d,a]*Q[b,c] for a,b in ((i,i+1),(i+1,i))
            for c,d in ((j,j+1),(j+1,j))) for j in range(L-1)] for i in range(L-1)])
        bond_error=float(np.max(abs(bc-wick)));assert max(density_error,bond_error)<1e-12
        for L in (48,96,144):
            ids=list(range(L//4,3*L//4));M=len(ids)
            for legs in (1,2):
                C,N=free_correlation(L,legs);n=np.diag(C);conn=np.diag(n)-abs(C)**2
                row=dict(kind=f'exact_XX_{legs}_leg',L=L,N=N,mu=float(-2*np.cos(np.pi*.296)),
                    density=spectral_metrics(n,conn[np.ix_(ids,ids)],ids))
                entropy=[]
                for cut in range(1,L):
                    ev=np.clip(np.linalg.eigvalsh(C[:cut,:cut]),1e-15,1-1e-15)
                    entropy.append(float(-np.sum(ev*np.log(ev)+(1-ev)*np.log1p(-ev))))
                row['S']=entropy;row['Sbar']={str(w):Z.sbar(entropy,L,w) for w in (1,6,12)}
                assert max(abs(np.array(entropy)-np.array(entropy)[::-1]))<1e-9
                if legs==1:
                    bids=list(range(L//4,3*L//4-1));Q=np.eye(L)-C
                    means=[-2*C[j,j+1] for j in range(L-1)]
                    bc=np.array([[sum(C[d,a]*Q[b,c] for a,b in ((i,i+1),(i+1,i))
                        for c,d in ((j,j+1),(j+1,j))) for j in bids] for i in bids])
                    row['bond']=spectral_metrics(means,bc,bids)
                rows.append(row)
            pattern=np.arange(L)%2
            pinned=np.zeros((M,M));sgn=(-1.)**np.array(ids)
            cat=.25*np.outer(sgn,sgn)
            rows.append(dict(kind='pinned_Neel_product',L=L,
                density=spectral_metrics(pattern,pinned,ids)))
            rows.append(dict(kind='translation_restored_Neel_cat',L=L,
                density=spectral_metrics(np.full(L,.5),cat,ids)))
            assert abs(rows[-2]['density']['profile_peak']['disconnected_over_M']-.25)<1e-12
            assert abs(rows[-1]['density']['connected_peak']['connected_over_M']-.25)<1e-12
    write_json(HERE/'controls.json',dict(rows=rows,slater_validation=dict(L=8,density_error=density_error,bond_error=bond_error),
        scope='Exact free controls at fixed hopping 1 and mu=-2cos(pi*.296), plus density-ordered examples; two legs have range two and are controls only.'))
    for r in rows:print(r['kind'],r['L'],r['density']['RMS'],r['density']['connected_peak']['connected_over_M'])


if __name__=='__main__':main()
