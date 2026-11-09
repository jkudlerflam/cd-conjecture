"""Second-order audit of a transverse DM/facilitated-spin-flip repair.

F = sum_j (n_j X_(j+1)-X_j n_(j+1)) = sum_j (n_(j-1)-n_(j+1))*X_j.
Periodic spin basis. Existing check_reduction builds T and adjacency R.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(key,'2')
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix,diags
from check_reduction import construct


def facilitated(L):
    dim=1<<L;r=[];c=[];v=[]
    for s in range(dim):
        for j in range(L):
            a=((s>>((j-1)%L))&1)-((s>>((j+1)%L))&1)
            if a:
                r.append(s^(1<<j));c.append(s);v.append(a)
    return coo_matrix((v,(r,c)),shape=(dim,dim),dtype=complex).tocsr()


def formulas(L,p):
    idx={int(s):i for i,s in enumerate(p)};fr=[];fc=[];fv=[];mr=[];mc=[];mv=[]
    for col,s0 in enumerate(p):
        s=int(s0);N=s.bit_count()
        n2=sum(((s>>j)&1)*((s>>((j+2)%L))&1) for j in range(L))
        fr.append(col);fc.append(col);fv.append(-2*N+2*n2)
        for j in range(L):
            k=(j+1)%L
            if ((s>>j)&1)!=((s>>k)&1):
                f=s^((1<<j)|(1<<k))
                if f in idx:
                    fr.append(idx[f]);fc.append(col);fv.append(1.)
            f=s^(1<<j)
            if f in idx:
                a=((s>>((j+2)%L))&1)-((s>>((j-2)%L))&1)
                if a:
                    mr.append(idx[f]);mc.append(col);mv.append(a)
    size=(len(p),len(p))
    return (coo_matrix((fv,(fr,fc)),shape=size,dtype=complex).tocsr(),
            coo_matrix((mv,(mr,mc)),shape=size,dtype=complex).tocsr())


def maximum(H):
    return float(np.max(np.abs(H.data))) if H.nnz else 0.


def main():
    a=argparse.ArgumentParser();a.add_argument('--max-L',type=int,default=12)
    a.add_argument('--output',default='facilitated_results.json');args=a.parse_args()
    t0=time.time();out={'checks':[]}
    for L in range(7,args.max_L+1):
        p,num,adj,T,W,K2,K3=construct(L,0.5j,1.9,periodic=True)
        F=facilitated(L);inv=np.zeros(len(adj));inv[adj>0]=1/adj[adj>0]
        R=diags(inv);f2=-F[p]@R@F[:,p]
        mixed=-T[p]@R@F[:,p]-F[p]@R@T[:,p]
        f2expected,mexpected=formulas(L,p)
        rec=dict(L=L,P_dimension=len(p),PFP_max=maximum(F[p][:,p]),
                 F_squared_error=maximum(f2-f2expected),mixed_error=maximum(mixed-mexpected),
                 mixed_max=maximum(mixed))
        out['checks'].append(rec);print(rec,flush=True)
        assert rec['PFP_max']<1e-12 and rec['F_squared_error']<1e-12 and rec['mixed_error']<1e-12
    out['elapsed_seconds']=time.time()-t0
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
