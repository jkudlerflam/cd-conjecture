"""Second-order pairing-only strong-repulsion Hamiltonian, in units Delta^2/V.
H_eff diagonal = -sum_{00 bonds} 1/(1+n_left+n_right) - m N.
Allowed two-site particle hops have amplitude -1/2.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json
import numpy as np
from scipy.sparse import coo_matrix
from checks import basis,low

def heff(L,N,Q=None):
    b=basis(L,N=N)
    nb=((b << 1)&((1 << L)-1)) | (b >> (L-1))
    b=b[(b&nb)==0]
    if Q is not None:
        q=sum((1 if j%2==0 else -1)*((b>>j)&1) for j in range(L))
        b=b[q==Q]
    if not len(b):return None,None
    diag=np.zeros(len(b));rows=[];cols=[];vals=[]
    index={int(x):i for i,x in enumerate(b)}
    for j in range(L):
        k=(j+1)%L
        mask=(((b>>j)&1)+((b>>k)&1))==0
        denom=1+((b>>((j-1)%L))&1)+((b>>((j+2)%L))&1)
        diag-=mask/denom
        k=(j+2)%L
        mask=((b>>j)&1)!=((b>>k)&1)
        for col,tar in zip(np.where(mask)[0],b[mask]^((1 << j)|(1 << k))):
            row=index.get(int(tar))
            if row is not None:rows.append(row);cols.append(int(col));vals.append(-.5)
    rows.extend(range(len(b)));cols.extend(range(len(b)));vals.extend(diag)
    return coo_matrix((vals,(rows,cols)),shape=(len(b),len(b))).tocsr(),b

if __name__=='__main__':
    results=[]
    for N in range(1,7):
        for Q in range(N%2,N+1,2):
            H,b=heff(24,N,Q);e,u=low(H,k=1)
            r=dict(L=24,N=N,Q=Q,dim=len(b),E_above_vacuum=float(e[0]+24),
                   E_per_particle=float((e[0]+24)/N))
            results.append(r);print(r,flush=True)
    open('singular_corner.json','w').write(json.dumps(results,indent=2))
