#!/usr/bin/env python3
"""Independent small-ring verification of the one-vacancy squeezed lift.

Uses arbitrary complex Hermitian internal bonds, no cell-model imports.
The local alphabet is vacancy 0 plus internal labels 1,...,q.
"""
from itertools import combinations, product
from pathlib import Path
import json
import numpy as np


def internal_bond(q, M, j, h):
    words = list(product(range(q), repeat=M))
    index = {w:i for i,w in enumerate(words)}
    out = np.zeros((q**M, q**M), complex)
    k = (j+1)%M
    for a,w in enumerate(words):
        col=q*w[j]+w[k]
        for b in range(q*q):
            if h[b,col] != 0:
                v=list(w); v[j],v[k]=divmod(b,q)
                out[index[tuple(v)],a]+=h[b,col]
    return out


def check(L,M,q,seed):
    rng=np.random.default_rng(seed)
    z=rng.normal(size=(q*q,q*q))+1j*rng.normal(size=(q*q,q*q))
    h=(z+z.conj().T)/2
    positions=list(combinations(range(L),M)); pi={x:i for i,x in enumerate(positions)}
    labels=list(product(range(q),repeat=M)); si={x:i for i,x in enumerate(labels)}
    C,S=len(positions),len(labels)
    words=[]
    for x in positions:
        for a in labels:
            w=[0]*L
            for j,i in enumerate(x):w[i]=a[j]+1
            words.append(tuple(w))
    wi={w:i for i,w in enumerate(words)}
    T=np.zeros((C*S,C*S),complex); V=np.zeros_like(T)
    Translation=np.zeros_like(T)
    for col,w in enumerate(words):
        Translation[wi[(w[-1],)+w[:-1]],col]=1
        for i in range(L):
            k=(i+1)%L
            if (w[i]==0) != (w[k]==0):
                a=list(w);a[i],a[k]=a[k],a[i]
                T[wi[tuple(a)],col]-=1
            if w[i] and w[k]:
                hc=q*(w[i]-1)+w[k]-1
                for out in range(q*q):
                    a=list(w);u,v=divmod(out,q);a[i],a[k]=u+1,v+1
                    V[wi[tuple(a)],col]+=h[out,hc]
    R=np.zeros((S,S),complex)
    for i,a in enumerate(labels):R[si[(a[-1],)+a[:-1]],i]=1
    bulk=np.zeros((C,C),complex); boundary=np.zeros_like(bulk)
    B=[np.zeros((C,C),complex) for _ in range(M)]
    for i,x in enumerate(positions):
        occupied=set(x)
        for a in range(L-1):
            if (a in occupied) != (a+1 in occupied):
                y=occupied.copy();y.symmetric_difference_update({a,a+1})
                bulk[pi[tuple(sorted(y))],i]-=1
        if L-1 in occupied and 0 not in occupied:
            y=occupied.copy();y.remove(L-1);y.add(0)
            boundary[pi[tuple(sorted(y))],i]-=1
        for j in range(M-1): B[j][i,i]=int(x[j+1]==x[j]+1)
        B[-1][i,i]=int(x[0]==0 and x[-1]==L-1)
    hs=[internal_bond(q,M,j,h) for j in range(M)]
    Hs=sum(hs)
    Ts=np.kron(bulk,np.eye(S))+np.kron(boundary,R)+np.kron(boundary.conj().T,R.conj().T)
    Vs=sum(np.kron(b,a) for b,a in zip(B,hs))
    exact_T=float(np.max(np.abs(T-Ts)))
    exact_V=float(np.max(np.abs(V-Vs)))
    assert exact_T<1e-12 and exact_V<1e-12
    # Keep every charge-sea ground state in every internal cyclic sector.
    columns=[]; expected_blocks=[]; sectors=[]; momenta=[]
    alpha=(M/L)**2-(np.sin(np.pi*M/L)/(L*np.sin(np.pi/L)))**2
    alpha*=L/M
    for r in range(M):
        theta=2*np.pi*r/M
        projector=sum(np.exp(-1j*theta*k)*np.linalg.matrix_power(R,k) for k in range(M))/M
        vals,vecs=np.linalg.eigh((projector+projector.conj().T)/2)
        W=vecs[:,vals>0.5]
        if W.shape[1]==0:continue
        Hc=bulk+np.exp(1j*theta)*boundary+np.exp(-1j*theta)*boundary.conj().T
        e,u=np.linalg.eigh(Hc); ground=u[:,e<e[0]+1e-9]
        # Resolve charge-ground degeneracy by the full physical translation.
        X=np.column_stack([np.kron(ground[:,a],W[:,b]) for a in range(ground.shape[1]) for b in range(W.shape[1])])
        trans=X.conj().T@Translation@X
        assert np.max(np.abs(trans.conj().T@trans-np.eye(len(trans))))<1e-9
        columns.append(X)
        expected_blocks.append(np.kron(np.eye(ground.shape[1]),alpha*(W.conj().T@Hs@W)))
        sectors.append({"r":r,"spin_dimension":W.shape[1],"charge_ground_degeneracy":ground.shape[1],"energy":float(e[0]),"adjacency_factor":float(alpha),"twist":float(theta)})
        momenta.append(np.linalg.eigvals(trans))
    X=np.column_stack(columns)
    result=X.conj().T@V@X
    expected=np.zeros_like(result);offset=0
    for block in expected_blocks:
        n=len(block);expected[offset:offset+n,offset:offset+n]=block;offset+=n
    error=float(np.max(np.abs(result-expected)))
    assert error<1e-8, (L,M,q,error)
    gram=float(np.max(np.abs(X.conj().T@X-np.eye(X.shape[1]))))
    assert gram<1e-9
    # Show explicitly that individual ordered-bond adjacency weights differ.
    eig,u=np.linalg.eigh(bulk+boundary+boundary.conj().T)
    v=u[:,0]
    raw_weights=[float(np.real(v.conj()@b@v)) for b in B]
    return {"L":L,"M":M,"q":q,"hilbert_dimension":C*S,"exact_hopping_error":exact_T,"exact_interaction_error":exact_V,"projected_first_order_error":error,"low_space_gram_error":gram,"alpha_exact":float(alpha),"ordered_bond_weights_at_theta_0":raw_weights,"sectors":sectors}


if __name__=='__main__':
    cases=[(4,2,2),(5,2,3),(5,3,2),(6,3,2),(6,4,2),(5,3,3)]
    result=[check(*args,101+i) for i,args in enumerate(cases)]
    output={"status":"passed","model":"arbitrary complex Hermitian internal bonds; t=1","cases":result}
    Path(__file__).with_name('results.json').write_text(json.dumps(output,indent=2)+'\n')
    for x in result: print(x['L'],x['M'],x['q'],x['hilbert_dimension'],x['projected_first_order_error'])
