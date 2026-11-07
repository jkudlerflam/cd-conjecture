"""Bethe doubling defect L(F_2L(2n)-2F_L(n))/M0 at larger L (defect only), Delta>1."""
import json, sys, numpy as np
from circle import Circ, solve
deltas=[1.001,1.01,1.1,1.5,2.,3.,5.,10.,20.]
out=open('scan_large.jsonl','w')
for d in deltas:
    M=Circ(d)
    for L in [127,128,256,512]:
        ns=sorted(set([1,2,3,int(L**(1/3))+1,L//8,L//4,L//2-4,L//2-2,L//2-1,L//2]))
        for n in ns:
            try:
                E,_=solve(L,n,d); E2,_=solve(2*L,2*n,d)
                row=dict(delta=d,L=L,n=n,defect=float(L*(E2-2*E)/M.M0))
            except Exception as ex:
                row=dict(delta=d,L=L,n=n,status='fail',err=str(ex))
            out.write(json.dumps(row)+'\n'); out.flush()
    print('done',d,flush=True)
