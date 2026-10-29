"""Exact rational certificate for the SU(2)-averaged stiffness operator.
6 B = 9 T - 16 S.S, in units t=1. Certify sum_{j=1}^3 B_j < 9/2.
All arithmetic in the LDL certificate uses Python exact fractions.
"""
from fractions import Fraction
from itertools import product
import json
import numpy as np

def certificate():
    groups={}
    for s in product(range(3),repeat=4):
        key=(s.count(1),s.count(2))
        groups.setdefault(key,[]).append(s)
    results=[]
    for key,states in groups.items():
        ix={s:i for i,s in enumerate(states)}; n=len(states)
        C=np.zeros((n,n),dtype=int)
        for col,s in enumerate(states):
            for j in range(3):
                a,b=s[j:j+2]
                if a and b:
                    C[col,col]+= -4 if a==b else 4
                    if a!=b:
                        r=list(s);r[j],r[j+1]=b,a; C[ix[tuple(r)],col]-=8
                elif a!=b:
                    r=list(s);r[j],r[j+1]=b,a; C[ix[tuple(r)],col]+=9
        A=27*np.eye(n,dtype=int)-C
        L=[[Fraction(int(i==j)) for j in range(n)] for i in range(n)]
        piv=[]
        for j in range(n):
            p=Fraction(int(A[j,j]))-sum(L[j][k]**2*piv[k] for k in range(j))
            assert p>0,(key,j,p)
            piv.append(p)
            for i in range(j+1,n):
                L[i][j]=(Fraction(int(A[i,j]))-sum(L[i][k]*L[j][k]*piv[k] for k in range(j)))/p
        for i in range(n):
            for j in range(n):
                assert sum(L[i][k]*piv[k]*L[j][k] for k in range(n))==int(A[i,j])
        results.append(dict(counts=key,dimension=n,pivots=[str(p) for p in piv],lambda_max_B=float(np.linalg.eigvalsh(C)[-1]/6)))
    return results
if __name__=='__main__':
    r=certificate()
    with open('tj_certificate.json','w') as f:json.dump(r,f,indent=2)
    print('EXACT PASS',len(r),'charge blocks, total dimension',sum(x['dimension'] for x in r))
    print('Numerical maximum B4/3:',max(x['lambda_max_B'] for x in r)/3)
    print('Certified A/(pi t) <= 3/2, ratio <=',1.5/np.log2(3))
