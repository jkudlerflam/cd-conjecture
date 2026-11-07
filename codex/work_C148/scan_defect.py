"""Matched doubling defect L(F_2L(2n)-2F_L(n))/M0 and finite-L diagnostics on Delta>1."""
import json, sys, numpy as np
from circle import *
deltas=[1.001,1.01,1.1,1.5,2.,3.,5.,10.,20.]
Ls=[int(a) for a in sys.argv[1].split(',')]
out=open(sys.argv[2],'w')
for d in deltas:
    M=Circ(d)
    for L in Ls:
        ns=sorted(set([1,2,L//4,L//3,(L//2)-2,(L//2)-1,L//2,int(L**(1/3))+1]))
        for n in ns:
            if not 1<=n<=L//2: continue
            try:
                E,_=solve(L,n,d); E2,_=solve(2*L,2*n,d)
                row=dict(delta=d,L=L,n=n,defect=float(L*(E2-2*E)/M.M0))
                row.update(finite_diag(L,n,d,grid=20001))
            except Exception as ex:
                row=dict(delta=d,L=L,n=n,status='fail',err=str(ex))
            out.write(json.dumps(row)+'\n'); out.flush()
    print('done',d,flush=True)
