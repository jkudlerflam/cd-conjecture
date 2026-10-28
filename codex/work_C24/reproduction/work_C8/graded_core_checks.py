"""Finite-cell checks of the graded-swap core, with no CFT inference."""
import itertools
import json
import math
import numpy as np
from scipy.sparse import csr_matrix, eye, kron
from boundary_ff_checks import chain


def code(n,N):
    q=n+1; D=2**n
    V=np.zeros((q**N,D),complex)
    for mask in range(D):
        colors=[j+1 for j in range(n) if mask>>j&1]
        r=len(colors); norm=math.sqrt(math.factorial(N)/math.factorial(N-r))
        for occupied in itertools.combinations(range(N),r):
            for perm in itertools.permutations(colors):
                letters=[0]*N
                for i,a in zip(occupied,perm): letters[i]=a
                inversions=sum(perm[i]>perm[j] for i in range(r) for j in range(i+1,r))
                V[np.ravel_multi_index(tuple(letters),(q,)*N),mask]=(-1)**inversions/norm
    return V


def one_case(n,N):
    q=n+1; D=2**n
    Pg=np.zeros((q*q,q*q))
    for a,b in itertools.product(range(q),repeat=2): Pg[b*q+a,a*q+b]=-1 if a*b else 1
    h=(np.eye(q*q)-Pg)/2
    bond_rank=int(np.linalg.matrix_rank(h))
    assert bond_rank==n*(n+3)//2
    H=chain(h,q,N,sparse=True); V=code(n,N)
    kernel_error=np.linalg.norm(H@V)
    isometry_error=np.linalg.norm(V.conj().T@V-np.eye(D))
    R=np.array([mask.bit_count() for mask in range(D)])
    P=np.diag((-1.)**R); g=np.diag(np.sqrt(N-R)/N)
    port_error=0.; clifford_approx=0.
    for a in range(n):
        f=np.zeros((D,D))
        for mask in range(D):
            if mask>>a&1: f[mask^(1<<a),mask]=(-1)**((mask&((1<<a)-1)).bit_count())
        ann=g@f
        Ax=ann+ann.T; Ay=-1j*(ann-ann.T)
        X=np.zeros((q,q),complex); X[0,a+1]=X[a+1,0]=1
        Y=np.zeros((q,q),complex); Y[0,a+1]=-1j; Y[a+1,0]=1j
        for F,A,B,Gam in [(X,Ax,-1j*Ay@P,f+f.T),(Y,Ay,1j*Ax@P,-1j*(f-f.T))]:
            for site,target in [(0,A),(N-1,B)]:
                Fphys=kron(csr_matrix(F),eye(q**(N-1))) if site==0 else kron(eye(q**(N-1)),csr_matrix(F))
                actual=V.conj().T@(Fphys@V)
                port_error=max(port_error,np.linalg.norm(actual-target))
            clifford_approx=max(clifford_approx,np.linalg.norm(np.sqrt(N)*A-Gam,2))
    predicted=1-np.sqrt(1-(n-1)/N)
    assert max(kernel_error,isometry_error,port_error,abs(clifford_approx-predicted))<1e-10
    out=dict(n=n,N=N,q=q,physical_dimension=q**N,code_dimension=D,bond_rank=bond_rank,
             kernel_error=float(kernel_error),isometry_error=float(isometry_error),port_error=float(port_error),
             normalized_clifford_error=float(clifford_approx),predicted_error=float(predicted))
    if q**N<=256:
        es=np.linalg.eigvalsh(H.toarray())
        ground_dimension=int(np.sum(es<1e-9)); assert ground_dimension==D
        out.update(ground_dimension=ground_dimension,gap=float(es[D]))
    return out


def main():
    out=[one_case(n,N) for n,N in [(1,3),(2,3),(2,5),(3,4),(3,5)]]
    with open('graded_core_checks.json','w') as f: json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
