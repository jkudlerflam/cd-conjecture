"""Rational scalar two-path resolvent, no floating point inference."""
import json
from fractions import Fraction as Q
import numpy as np


def inverse(a):
    n=len(a);b=np.array([[Q(int(i==j)) for j in range(n)] for i in range(n)],dtype=object)
    a=np.array([[Q(int(x)) for x in row] for row in a],dtype=object)
    for j in range(n):
        pivot=next(i for i in range(j,n) if a[i,j])
        if pivot!=j:a[[j,pivot]]=a[[pivot,j]];b[[j,pivot]]=b[[pivot,j]]
        val=a[j,j];a[j]=a[j]/val;b[j]=b[j]/val
        for i in range(n):
            if i!=j and a[i,j]:
                val=a[i,j];a[i]-=val*a[j];b[i]-=val*b[j]
    return b


def weights(n,d):
    m=n+1;A=d*np.eye(m,dtype=int)-np.diag(np.ones(n,dtype=int),1)-np.diag(np.ones(n,dtype=int),-1)
    R=inverse(np.kron(A,np.eye(m,dtype=int))+np.kron(np.eye(m,dtype=int),A))
    v=R[:,n*m];F=[]
    for i in range(n):
        row=[];left=np.arange(i*m,(i+1)*m)
        for j in range(n):
            right=np.arange(j,m*m,m)
            row.append(2*sum(v[x]*R[x,y]*v[y] for x in left for y in right))
        F.append(row)
    return F

if __name__=='__main__':
    rows=[]
    for n in (2,3):
        for d in (3,5):
            F=weights(n,d)
            assert F[1][0]>F[0][1]
            rows.append(dict(n=n,d=d,F=[[str(x) for x in row] for row in F],F10_minus_F01=str(F[1][0]-F[0][1])))
    print(json.dumps(rows,indent=2))
    with open('exact_weights_results.json','w') as f:json.dump(rows,f,indent=2)
