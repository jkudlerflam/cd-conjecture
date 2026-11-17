"""Exact small open-interval recursions for the first resonant hop."""
import json
import time
from pathlib import Path
import numpy as np
import sympy as sp

from model import embed


def exact(A):
    A=A.tocoo()
    return sp.MutableSparseMatrix(A.shape[0],A.shape[1],
                                  {(int(i),int(j)):sp.Integer(round(v.real))
                                   for i,j,v in zip(A.row,A.col,A.data) if v!=0})


def canon(x):
    return sp.expand(sp.radsimp(sp.expand(x)))


def take(A,rows,cols=None):
    return A.extract(rows,rows if cols is None else cols)


def schur(H,keep):
    remove=sorted(set(range(H[0].rows))-set(keep))
    R=sp.diag(*[canon(1/H[0][j,j]) for j in remove])
    B=[take(h,remove,keep) for h in H];C=[take(h,remove) for h in H]
    X=[sp.zeros(len(remove),len(keep))];out=[take(H[0],keep)]
    for n in range(1,len(H)):
        a=B[n]
        for k in range(1,n):a-=C[k]*X[n-k]
        X.append((R*a).applyfunc(canon))
        a=take(H[n],keep)
        for k in range(1,n):a-=B[k].T*X[n-k]
        out.append(a.applyfunc(canon))
    return out


def run(m):
    q={3:sp.Integer(3),4:sp.Integer(1),5:5-2*sp.sqrt(5),6:sp.Rational(1,3)}[m]
    L=m+1;order=6 if m<=4 else 8;dim=1<<L
    v0=np.array([0,0,0,1]);v1=np.array([0,-1,1,0]);v2=np.array([1,0,0,0])
    outer=np.outer
    H=[exact(embed(L,outer(v0,v0),False)),
       exact(embed(L,outer(v0,v1)+outer(v1,v0),False)),
       exact(embed(L,outer(v1,v1),False))+q*exact(embed(L,outer(v0,v2)+outer(v2,v0),False)),
       q*exact(embed(L,outer(v1,v2)+outer(v2,v1),False)),
       q*q*exact(embed(L,outer(v2,v2),False))]+[sp.zeros(dim) for _ in range(order-4)]
    words=list(range(dim));stages=[]
    for n in range(0,order,2):
        H=[h.applyfunc(canon) for h in H]
        keep=[i for i in range(len(words)) if H[0][i,i]==0]
        assert all(i==j or val==0 for (i,j),val in H[0].todok().items())
        H=schur(H,keep);words=[words[i] for i in keep]
        assert all(v==0 for v in H[1].todok().values())
        stages.append({'order':n,'kernel_dimension':len(words)})
        if n<order-2:H=H[2:]
    K=H[2];t=-canon(K[words.index(1<<m),words.index(1)])
    N=[s.bit_count() for s in words]
    mixing=[v for (i,j),v in K.todok().items() if N[i]!=N[j]]
    assert all(v==0 for v in mixing)
    eigen_dimension=len(words)-K.rank()
    assert eigen_dimension==L+1
    return {'m':m,'q':str(q),'open_interval_length':L,'first_nonzero_order':order,
            'hop_m_magnitude':str(t),'hop_m_decimal':float(t.evalf()),
            'all_lower_kept_orders_zero':True,'leading_number_conservation':True,
            'leading_kernel_dimension':eigen_dimension,'stages':stages}


def main():
    start=time.perf_counter();rows=[]
    for m in (3,4,5,6):
        rows.append(run(m));print(rows[-1],flush=True)
    out={'sympy_version':sp.__version__,'results':rows,'runtime_seconds':time.perf_counter()-start}
    Path(__file__).with_name('exact_local_results.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
