#!/usr/bin/env python3
"""Analytic sine-wave correlation matrices for the finite XX controls."""
import json
import numpy as np
from threadpoolctl import threadpool_limits
from common import HERE,Ana


def main():
    rows=[]
    with threadpool_limits(limits=1):
        for legs in (1,2):
            for L in (32,64,128):
                n=L//legs;N=round(.19*n);q=np.arange(1,N+1)
                mu=-np.cos(np.pi*N/(n+1))-np.cos(np.pi*(N+1)/(n+1))
                U=np.sqrt(2/(n+1))*np.sin(np.pi*np.arange(1,n+1)[:,None]*q/(n+1));C=U@U.T
                ss=[0.]
                for l in range(1,n):
                    v=np.clip(np.linalg.eigvalsh(C[:l,:l]),1e-15,1-1e-15)
                    ss.append(float(-np.sum(v*np.log(v)+(1-v)*np.log1p(-v))))
                ss.append(0.)
                S=ss[1:-1] if legs==1 else [ss[(l+1)//2]+ss[l//2] for l in range(1,L)]
                assert np.max(abs(np.array(S)-Ana.ff_entropies(L,mu,legs)))<1e-10
                levels=-2*np.cos(np.pi*np.arange(1,n+1)/(n+1))-mu
                assert np.count_nonzero(levels<0)==N and min(abs(levels))>1e-6
                rows.append(dict(legs=legs,L=L,N=N*legs,mu=float(mu),E=float(legs*levels[:N].sum()),S=S,
                    Sbar={str(w):Ana.sbar(S,L,w) for w in (1,6,12)},c_profile=Ana.c_prof(S,L)))
    pairs=[]
    for legs in (1,2):
        for L in (32,64):
            u=next(r for r in rows if r['legs']==legs and r['L']==L)
            v=next(r for r in rows if r['legs']==legs and r['L']==2*L)
            pairs.append(dict(legs=legs,lengths=[L,2*L],c_two={str(w):6*(v['Sbar'][str(w)]-u['Sbar'][str(w)])/np.log(2) for w in (1,6,12)}))
    (HERE/'exact_controls.json').write_text(json.dumps(dict(rows=rows,pairs=pairs),indent=2)+'\n')
    for row in pairs:print(row)


if __name__=='__main__':main()
