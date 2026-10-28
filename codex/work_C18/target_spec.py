"""Uniform route-A target specification, not a local history-machine compiler.
The physical period b is the sole input. Write r=floor(log2(b-1)/10).
The r=0 case is a trivial one-dimensional logical target.
"""
import json
from pathlib import Path
import numpy as np

I=np.eye(2); X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]]); Z=np.diag([1,-1])


def kron(xs):
    out=np.array([[1.]])
    for x in xs: out=np.kron(out,x)
    return out


def primes(count):
    ps=[]; x=2
    while len(ps)<count:
        if all(x%p for p in ps if p*p<=x): ps.append(x)
        x+=1
    return ps


def target_from_b(b):
    assert b>=3
    n=b-2; r=((n+1).bit_length()-1)//10
    if r==0: return dict(b=b,n=n,r=0,N=1,K=np.zeros((1,1)),A=[np.ones((1,1))],B=[np.ones((1,1))])
    N=2**r; ps=primes(2*r); den=10*r*ps[-1]
    theta=np.pi*np.sqrt(ps[::2])/den
    kF=np.pi/3+np.pi*np.sqrt(ps[1::2])/den
    mu=-2*np.cos(kF)
    gam=[]; zz=[]
    P=kron([Z]*r)
    for f in range(r):
        gam.extend([kron([Z]*f+[v]+[I]*(r-f-1)) for v in (X,Y)])
        zz.append(kron([I]*f+[Z]+[I]*(r-f-1)))
    eta=[1j*g@P for g in gam]
    A=eta+zz+[np.eye(N)]; B=gam+zz+[np.eye(N)]
    m=3*r+1; K=np.zeros((m,m))
    for f in range(r):
        o,e=2*f,2*f+1
        K[o,o]=K[e,e]=-np.sin(theta[f])/2
        K[o,e]=-np.cos(theta[f])/2; K[e,o]=np.cos(theta[f])/2
        K[2*r+f,-1]=K[-1,2*r+f]=mu[f]/4
    D=4*((((18*r+19)*N)+3)//4)
    return dict(b=b,n=n,r=r,N=N,K=K,A=A,B=B,theta=theta,kF=kF,mu=mu,
                C9_operator_slots=m,C9_energy_coordinates=m+2,
                abstract_cell_dimension_bound=D)


def verify(b):
    spec=target_from_b(b); r,N=spec['r'],spec['N']
    K,A,B=spec['K'],spec['A'],spec['B']
    H=sum(K[i,j]*np.kron(A[i],B[j]) for i,j in zip(*np.nonzero(K)))
    # Independently construct Jordan-Wigner fermions on two ordered cells.
    ferm=[]
    for f in range(2*r):
        ferm.append(kron([Z]*f+[(X+1j*Y)/2]+[I]*(2*r-f-1)))
    Href=np.zeros((N*N,N*N),complex)
    for f in range(r):
        a,c=ferm[f],ferm[f+r]; theta=spec['theta'][f]; mu=spec['mu'][f]
        Href-=np.exp(1j*theta)*a.conj().T@c+np.exp(-1j*theta)*c.conj().T@a
        # Half of each onsite term belongs to this single bond.
        Href-=mu/2*(a.conj().T@a+c.conj().T@c-np.eye(N*N))
    err=float(np.linalg.norm(H-Href))
    assert err<1e-12
    # C9 gives a fixed algorithm for coefficients. Verify its dimension
    # requirement here, not its previously proved spectral dilation again.
    assert spec['abstract_cell_dimension_bound']>=max(16*N,(6*(3*r+3)+1)*N)
    return dict(b=b,n=spec['n'],r=r,N=N,operator_slots=3*r+1,
                abstract_dimension_bound=spec['abstract_cell_dimension_bound'],
                work_qubits_for_dimension=(spec['abstract_cell_dimension_bound']//4-1).bit_length(),
                two_cell_JW_error=err,theta=spec['theta'].tolist(),kF=spec['kF'].tolist())

if __name__=='__main__':
    # Large b values are input integers only. No physical chain is allocated.
    out=[verify(2**(10*r)+1) for r in (1,2,3)]
    Path('target_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
