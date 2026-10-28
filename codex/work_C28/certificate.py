"""Exact rational certificate of the four n=2 leading Fermi crossings.

Fractions is Python's standard library, not an additional numerical package.
"""
from fractions import Fraction as F
import json

def inverse(A):
    n=len(A);B=[[F(x) for x in row]+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for i in range(n):
        k=next(k for k in range(i,n) if B[k][i]);B[i],B[k]=B[k],B[i]
        q=B[i][i];B[i]=[x/q for x in B[i]]
        for k in range(n):
            if k!=i:
                q=B[k][i];B[k]=[a-q*b for a,b in zip(B[k],B[i])]
    return [row[n:] for row in B]

def main():
    m=3;A=[[0]*9 for _ in range(9)]
    for i in range(m):
        for j in range(m):
            q=i*m+j;A[q][q]=8
            for u,v in [(i-1,j),(i+1,j),(i,j-1),(i,j+1)]:
                if 0<=u<m and 0<=v<m:A[q][u*m+v]=-1
    R=inverse(A);s=6;v=[r[s] for r in R];Fs=[]
    for i,j in [(0,0),(0,1),(1,0)]:
        z=[v[q] if q%m==j else F(0) for q in range(9)]
        rz=[sum(a*b for a,b in zip(row,z)) for row in R]
        Fs.append(2*sum(v[q]*rz[q] for q in range(9) if q//m==i))
    a,b,c=Fs;eps=F(1,200)
    # cos(1.4)<.18 and >0 follow from alternating Taylor bounds.
    x=F(7,5);cos_upper=1-x*x/2+x**4/24-x**6/720+x**8/40320
    cos_lower=cos_upper-x**10/3628800
    assert 0<cos_lower<cos_upper<F(18,100)
    upper_mu_a=8*a*(b+c)*eps
    lower_bc_margin=F(82,100)*4*b*c*(1-eps*eps)
    assert 0<upper_mu_a<lower_bc_margin
    out=dict(F00=str(a),F01=str(b),F10=str(c),F_decimal=[float(x) for x in Fs],
             upper_mu_abs_a=str(upper_mu_a),lower_bc_one_minus_cos=str(lower_bc_margin),
             strict_margin_ratio=float(lower_bc_margin/upper_mu_a),
             cos14_lower=str(cos_lower),cos14_upper=str(cos_upper))
    with open('certificate_results.json','w') as f:json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
