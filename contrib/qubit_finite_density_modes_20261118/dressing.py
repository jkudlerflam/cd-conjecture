"""Exact integer verification of the leading dressed-density operator."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import sys
from pathlib import Path
import json
import time
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'qubit_cancellation_classification_20261116'))
import algebra as A


def check(L):
    ns=A.occupations(L);energies=np.sum(ns*np.roll(ns,-1,axis=1),axis=1)
    P=np.flatnonzero(energies==0);Q=np.flatnonzero(energies>0)
    pair=A.chain(L,A.bond(p=1))
    B=A.sub(pair,Q,P).real.astype(np.int64)
    reached=np.unique(B.nonzero()[0])
    assert set(energies[Q[reached]]) <= {1,2,3}
    weights=np.zeros(len(Q),dtype=np.int64)
    weights[reached]=36//energies[Q[reached]]**2
    exact36=B.T@sparse.diags(weights,dtype=np.int64)@B
    index={int(s):i for i,s in enumerate(P)}
    rows=[];cols=[];vals=[]
    number=np.sum(ns[P],axis=1)
    for col,s in enumerate(P):
        n=ns[s];N=int(sum(n));D2=int(sum(n*np.roll(n,-2)));D3=int(sum(n*np.roll(n,-3)))
        rows.append(col);cols.append(col);vals.append(36*L-126*N+54*D2+22*D3)
        for j in range(L):
            k=(j+2)%L
            if n[j] and not n[k]:
                dest=int(s)^(1<<j)^(1<<k)
                if dest in index:
                    row=index[dest]
                    rows.extend((row,col));cols.extend((col,row));vals.extend((9,9))
    target=sparse.coo_matrix((np.array(vals,dtype=np.int64),(rows,cols)),shape=exact36.shape).tocsr()
    difference=exact36-target;difference.eliminate_zeros()
    assert difference.nnz==0
    # Independent finite spectra supplement, but do not replace, the bound proof.
    lower=float(np.linalg.eigvalsh(exact36.toarray()/36)[0])
    upper=float(np.linalg.eigvalsh(np.diag(L-8*number/9)-exact36.toarray()/36)[0])
    assert min(lower,upper)>-1e-10
    return dict(L=L,dimension=len(P),integer_identity_nonzeros=difference.nnz,
                minimum_B=lower,minimum_upper_bound_difference=upper)


if __name__=='__main__':
    start=time.monotonic()
    with threadpool_limits(1):out=dict(checks=[check(L) for L in (8,9,10,12)])
    out['runtime_seconds']=time.monotonic()-start
    (HERE/'dressing_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
