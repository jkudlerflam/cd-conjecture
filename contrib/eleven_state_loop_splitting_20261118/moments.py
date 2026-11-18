"""Exact rational endpoint-moment criterion; no unbounded-n extrapolation."""
from fractions import Fraction as F
from math import factorial


def inverse(A):
    n=len(A)
    M=[list(row)+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for j in range(n):
        pivot=M[j][j]
        M[j]=[x/pivot for x in M[j]]
        for i in range(n):
            if i!=j:
                factor=M[i][j]
                M[i]=[x-factor*y for x,y in zip(M[i],M[j])]
    return [row[n:] for row in M]


def bilinear(x,A,y):
    return sum(x[i]*A[i][j]*y[j] for i in range(len(x)) for j in range(len(y)))


def criterion(n):
    # Time rescaled so 2a=1. Other positive a give the same nonzero test.
    C=[[F(factorial(2*i+2*j+2),factorial(2*i+1)*factorial(2*j+1))
        for j in range(n)] for i in range(n)]
    Q=[[C[i][j]*F(i-j,2*(i+j+1)) for j in range(n)] for i in range(n)]
    J=[[F(0) for j in range(n)] for i in range(n)]
    if n>=2:J[-2][-1]=F(1);J[-1][-2]=F(-1)
    X=[];norms=[]
    for j in range(n):
        tail=inverse([row[j:] for row in C[j:]])
        x=[F(0)]*j+[row[0]/tail[0][0] for row in tail]
        X.append(x);norms.append(bilinear(x,C,x))
    for j in range(n):
        assert norms[j]>0
        for l in range(j+1,n):assert bilinear(X[j],C,X[l])==0
    L=[sum(bilinear(X[j],Q,X[l])*bilinear(X[j],J,X[l])/norms[l]
           for l in range(j+1,n)) for j in range(n-1)]
    if n>=2:assert L[-1]==-F(1,16*n-12)
    return C,Q,J,X,norms,L


if __name__=='__main__':
    for n in (2,3,4):
        print(n,[str(x) for x in criterion(n)[-1]])
