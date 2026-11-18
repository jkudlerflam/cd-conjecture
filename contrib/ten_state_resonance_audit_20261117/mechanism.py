"""Nine-state single rail: its scalar reduction and literal generated vertices."""
import numpy as np
from scipy import sparse
from scipy.linalg import solve


def fermions(m):
    states=np.arange(1<<m);out=[]
    for i in range(m):
        col=states[(states>>i)&1==1];row=col^(1<<i)
        val=np.array([(-1)**int(x&((1<<i)-1)).bit_count() for x in col])
        out.append(sparse.csr_matrix((val,(row,col)),shape=(1<<m,)*2))
    return out


def scalar(n,mass=2):
    s=n+2;H=mass*np.eye(s)-np.eye(s,k=1)-np.eye(s,k=-1)
    R=np.linalg.inv(np.kron(H,np.eye(s))+np.kron(np.eye(s),H))
    source=np.zeros(s*s);source[(s-1)*s]=1
    x=(R@source).reshape(s,s);U=np.zeros((s*s,n));V=U.copy()
    for i in range(n):U[i*s:(i+1)*s,i]=x[i,:];V[i::s,i]=x[:,i]
    C=U.T@R@V;A=U.T@R@U+V.T@R@V
    return H,R,-2*C,A


def pair_coefficients(n,mass=2,max_order=4):
    """K_(2r)=-S^T R (C R)^(2r) S, excluding lambda² alpha⁴ g^(2r).

    Carrier species order AA,AB,BA,BB, then pair positions, then Fock.
    Source BA - AB*P_left; conversion B<-A annihilates a local-cell bit.
    This is derived from the literal single-rail signs, not a chosen interaction.
    """
    H,R,T,A=scalar(n,mass);s=len(H);d=1<<(2*n);f=fermions(2*n)
    parity=sparse.diags([(-1)**(x&((1<<n)-1)).bit_count() for x in range(d)],format='csr',dtype=float)
    local_right=[parity@f[n+i] for i in range(n)]
    C=sparse.csr_matrix((4*s*s*d,)*2)
    for i in range(n):
        pi=np.zeros((s,s));pi[i,i]=1
        posleft=sparse.kron(pi,np.eye(s));posright=sparse.kron(np.eye(s),pi)
        for right in (0,1):
            E=np.zeros((4,4));E[2+right,right]=1
            term=sparse.kron(E,sparse.kron(posleft,f[i]),format='csr');C+=term+term.T
        for left in (0,1):
            E=np.zeros((4,4));E[2*left+1,2*left]=1
            term=sparse.kron(E,sparse.kron(posright,local_right[i]),format='csr');C+=term+term.T
    S=np.zeros((4*s*s*d,d));idx=(s-1)*s
    S[(2*s*s+idx)*d:(2*s*s+idx+1)*d]=np.eye(d)
    S[(s*s+idx)*d:(s*s+idx+1)*d]=-parity.toarray()
    def applyR(X):return np.einsum('ab,sbcm->sacm',R,X.reshape(4,s*s,d,d)).reshape(4*s*s*d,d)
    X=applyR(S);out={0:-S.T@X}
    for order in range(1,max_order+1):
        X=applyR(C@X)
        if order%2==0:out[order]=-S.T@X
        else:assert np.max(abs(S.T@X))<1e-12
    return out


def normal_order(K,m):
    """Constant, quadratic, and antisymmetrized two-particle vertex."""
    f=fermions(m);d=1<<m;k=K[0,0];inds=[1<<i for i in range(m)]
    h=K[np.ix_(inds,inds)]-k*np.eye(m);G=k*np.eye(d)
    for i in range(m):
        for j in range(m):G+=h[i,j]*(f[i].T@f[j]).toarray()
    pairs=[(i,j) for i in range(m) for j in range(i+1,m)]
    pairinds=[(1<<i)+(1<<j) for i,j in pairs]
    W=(K-G)[np.ix_(pairinds,pairinds)]
    quartic=np.zeros_like(K)
    for a,(i,j) in enumerate(pairs):
        for b,(k1,l) in enumerate(pairs):
            quartic+=W[a,b]*(f[i].T@f[j].T@f[l]@f[k1]).toarray()
    return dict(constant=k,h=h,W=W,pairs=pairs,gaussian_residual=np.linalg.norm(K-G,2),higher_residual=np.linalg.norm(K-G-quartic,2))


def projected_controls(n):
    """Projection of literal N, sum P_L1, and LL/RR bit swaps."""
    s=n+2;v=np.sqrt(2/(s+1))*np.sin(np.pi*np.arange(1,s+1)/(s+1))
    left=np.array([sum(v[i+1:]**2) for i in range(n)])
    swap=np.diag(1-v[1:n]**2,k=1)+np.diag(1-v[1:n]**2,k=-1)
    return left,swap


def resonant_seed(T,center=1/5,radius=2/3):
    """An explicit REAL onsite matrix with four roots and reciprocal resonance.

    A congruence is used only to prove the zeros; physical eigenvectors and
    velocities must be calculated from the unwhitened Bloch symbol.
    """
    assert T.shape==(2,2)
    S=(T+T.T)/2;A=(T-T.T)/2;det=np.linalg.det(S)
    assert det<0 and np.linalg.det(T)<0
    absS=(np.trace(S)*S-2*det*np.eye(2))/np.sqrt(np.trace(S)**2-4*det)
    rho2=A[0,1]**2/(-det);assert 0<rho2<1
    x0=(1-rho2)*center
    epsilon=2*np.sqrt(rho2+(1-rho2)*(radius**2-rho2*center**2))
    M=epsilon*absS-2*x0*S
    left,swap=projected_controls(2)
    eta=(M[0,0]-M[1,1])/(left[0]-left[1]);mu=M[0,0]-eta*left[0];kappa=M[0,1]/swap[0,1]
    assert np.max(abs(M-(mu*np.eye(2)+eta*np.diag(left)+kappa*swap)))<1e-12
    return dict(M=M,rho2=rho2,epsilon=epsilon,x0=x0,mu=mu,eta=eta,kappa=kappa,
                cos_roots=[center-radius,center+radius])


def bloch(T,M,k):return M+T*np.exp(1j*k)+T.T*np.exp(-1j*k)


def fermi_modes(T,M,cos_roots):
    modes=[]
    for k in sorted([sign*np.arccos(x) for x in cos_roots for sign in (-1,1)]):
        E,U=np.linalg.eigh(bloch(T,M,k));i=np.argmin(abs(E));u=U[:,i]
        v=float((u.conj()@(1j*T*np.exp(1j*k)-1j*T.T*np.exp(-1j*k))@u).real)
        modes.append(dict(k=float(k),v=v,u=u,residual=float(abs(E[i]))))
    return modes


def vertex(W,pairs,modes,a,b,c,d):
    """Coefficient of psi_a† psi_b† psi_d psi_c from the literal bond vertex.

    Each local one-particle wave has components (u,e^(ik)u) on its two cells.
    Antisymmetric wedge factors include direct and exchange contributions.
    """
    def wave(x):return np.r_[modes[x]['u'],np.exp(1j*modes[x]['k'])*modes[x]['u']]
    def wedge(x,y):
        u=wave(x);v=wave(y);return np.array([u[i]*v[j]-u[j]*v[i] for i,j in pairs])
    return np.vdot(wedge(a,b),W@wedge(c,d))
