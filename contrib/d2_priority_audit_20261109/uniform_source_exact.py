"""Exact Gaussian-integer check of the sixth-order canonical SW coefficient.

X_n = A_n / D**n for invariant-graph coefficients, with D=lcm(1,...,L).
Store real and imaginary parts as Python integer matrices; no floating point
is used after the integer microscopic matrices are constructed.
"""
from pathlib import Path
import sys, json, math
from fractions import Fraction
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'codex/work_C186'))
import sw_general as sw

def z(A): return (A.real.astype('int64').astype(object), A.imag.astype('int64').astype(object))
def zero(shape): return (np.zeros(shape,dtype=object),np.zeros(shape,dtype=object))
def add(A,B): return (A[0]+B[0],A[1]+B[1])
def sub(A,B): return (A[0]-B[0],A[1]-B[1])
def scale(a,A): return (a*A[0],a*A[1])
def dot(A,B): return (A[0]@B[0]-A[1]@B[1],A[0]@B[1]+A[1]@B[0])
def adj(A): return (A[0].T,-A[1].T)
def comm(A,B): return sub(dot(A,B),dot(B,A))

def exact(L,p,f,a,b,m):
    P=sw.pspace(L,True); Q=[s for s in range(2**L) if s not in P]
    B=sw.build_full(L,sw.bond_matrix(0,0,p,-f,f),sw.site_matrix(0,0),True)
    C=sw.build_full(L,sw.bond_matrix(0,a+1j*b,0,0,0),sw.site_matrix(m,0),True)
    bpq,bqp,bqq=(z(B[np.ix_(i,j)]) for i,j in ((P,Q),(Q,P),(Q,Q)))
    cpp,cpq,cqp,cqq=(z(C[np.ix_(i,j)]) for i,j in ((P,P),(P,Q),(Q,P),(Q,Q)))
    D=math.lcm(*range(1,L+1)); inv=np.array([D//sw.adj_count(s,L,True) for s in Q],dtype=object)[:,None]
    A=[zero((len(Q),len(P))) for _ in range(7)]
    for n in range(1,7):
        R=bqp if n==1 else scale(D,cqp) if n==2 else zero((len(Q),len(P)))
        R=add(R,dot(bqq,A[n-1]))
        if n>=2: R=add(R,scale(D,sub(dot(cqq,A[n-2]),dot(A[n-2],cpp))))
        for i in range(1,n):
            j=n-1-i
            if j>=1: R=sub(R,dot(dot(A[i],bpq),A[j]))
            j=n-2-i
            if j>=1: R=sub(R,scale(D,dot(dot(A[i],cpq),A[j])))
        A[n]=scale(-inv,R)
    # H_B,n = B_n / D**(n-1), U_n = U_numerator / D**n.
    HB=[zero((len(P),len(P))) for _ in range(7)]
    U=[zero((len(P),len(P))) for _ in range(7)]
    for n in range(2,7):
        HB[n]=add(dot(bpq,A[n-1]),scale(D,dot(cpq,A[n-2])))
        if n==2: HB[n]=add(HB[n],scale(D,cpp))
        for i in range(1,n): U[n]=add(U[n],dot(adj(A[i]),A[n-i]))
    pos={s:i for i,s in enumerate(P)}; row,col=pos[4],pos[0]
    coeffs=[]; hermiticity=[]
    for n in range(2,7):
        # 8 D**(n-1) times the canonical coefficient.
        H=scale(8,HB[n])
        for i in range(2,n-1): H=add(H,scale(4,comm(U[i],HB[n-i])))
        if n==6:
            UU=dot(U[2],U[2])
            H=add(H,add(sub(scale(3,dot(HB[2],UU)),dot(UU,HB[2])),scale(-2,dot(dot(U[2],HB[2]),U[2]))))
        coeffs.append([str(Fraction(H[k][row,col],8*D**(n-1))) for k in (0,1)])
        defect=sub(H,adj(H)); hermiticity.append(all(not q for part in defect for q in part.flat))
    expected=2*p**3*f*b
    assert coeffs[:-1]==[['0','0']]*4, coeffs
    assert coeffs[-1]==['0',str(expected)], (coeffs,expected)
    assert all(hermiticity)
    return dict(L=L,p=p,f=f,a=a,b=b,m=m,source_scaled_orders_2_to_6=coeffs,expected_order6_imag=expected,exact_hermiticity=hermiticity)

if __name__=='__main__':
    cases=[(L,p,f,a,1,m) for L in (6,7,8) for p,f,a,m in
           [(1,1,0,0),(2,1,0,0),(1,2,0,0),(1,1,1,0),(1,1,0,1)]]
    rows=[]
    for case in cases:
        row=exact(*case); rows.append(row); print(json.dumps(row),flush=True)
    Path(__file__).with_suffix('.json').write_text(json.dumps({'rows':rows},indent=2)+'\n')
