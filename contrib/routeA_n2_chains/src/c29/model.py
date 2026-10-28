"""C29 balanced multi-bus model. NumPy/SciPy only.
Fixed legal cell sector, no projection of vacuum-walker excitations.
"""
import numpy as np
from scipy.linalg import eigh, block_diag
from scipy.integrate import quad


def kron_all(xs):
    ans=np.ones((1,1),complex)
    for x in xs: ans=np.kron(ans,x)
    return ans


def data_ops(n):
    I=np.eye(2); Z=np.diag([1.,-1.]); sm=np.array([[0.,1.],[0.,0.]])
    cs=[kron_all([Z if i<j else sm if i==j else I for i in range(n)]) for j in range(n)]
    return cs,kron_all([Z]*n)


def excitation(n,delta,phases,g=.2,t=1.):
    N=2**n; m=n+1; cs,P=data_ops(n)
    H=delta*np.eye(3*m*N,dtype=complex)
    def sl(j,a): return slice((3*j+a)*N,(3*j+a+1)*N)
    phases=np.asarray(phases)
    for j in range(n):
        for a in range(3):
            val=-t*np.exp(1j*phases[j,a])*np.eye(N)
            H[sl(j+1,a),sl(j,a)]+=val
            H[sl(j,a),sl(j+1,a)]+=val.conj().T
        H[sl(j,1),sl(j,0)]+=g*cs[j]
        H[sl(j,0),sl(j,1)]+=g*cs[j].conj().T
        H[sl(j,2),sl(j,0)]+=g*cs[j].conj().T
        H[sl(j,0),sl(j,2)]+=g*cs[j]
    sources=[]
    for j in (0,n):
        row=[]
        for a in range(3):
            S=np.zeros((len(H),N),complex)
            S[sl(j,a)]=P if j==n and a==0 else np.eye(N)
            row.append(S)
        sources.append(row)
    return H,sources


def bond(n,delta,phases,phi,g=.2,coefficient=False,inverse_power=1):
    H,(L,R)=excitation(n,delta,phases,g=0. if coefficient else g)
    e,U=eigh(H); L=[U.conj().T@s for s in L];R=[U.conj().T@s for s in R]
    W=np.kron(R[0],L[0])+np.exp(1j*phi)*np.kron(R[1],L[2])+np.exp(-1j*phi)*np.kron(R[2],L[1])
    den=e[:,None]+e[None,:]
    if not coefficient: return -W.conj().T@(W/den.reshape(-1,1)**inverse_power)
    C=U.conj().T@(excitation(n,delta,phases,g=1.)[0]-H)@U
    X=(W/den.reshape(-1,1)).reshape(len(H),len(H),-1)
    Y=(np.einsum('ab,bck->ack',C,X)+np.einsum('ab,cbk->cak',C,X)).reshape(len(H)**2,-1)
    return -Y.conj().T@(Y/den.reshape(-1,1))


def heat_weights(n,delta,t=1.):
    m=n+1
    A=delta*np.eye(m)-t*(np.diag(np.ones(n),1)+np.diag(np.ones(n),-1))
    e,U=eigh(A)
    def f(s,i,tau):
        w=U[s]*U[i]; ea=e[:,None]; eb=e[None,:]; diff=ea-eb
        same=abs(diff)<1e-12
        arr=(np.exp(-eb*tau)-np.exp(-ea*tau))/np.where(same,1.,diff)
        arr[same]=np.broadcast_to(tau*np.exp(-ea*tau),diff.shape)[same]
        return np.sum(w[:,None]*w[None,:]*arr)
    return np.array([[quad(lambda tau:f(n,i,tau)*f(0,j,tau),0,np.inf,epsabs=1e-16,epsrel=1e-10)[0] for j in range(n)] for i in range(n)])


def formula(n,delta,phases,phi):
    phases=np.asarray(phases)
    kap=np.r_[0.,np.cumsum(phases[:,0]-(phases[:,1]+phases[:,2])/2)]
    chi=np.r_[0.,np.cumsum((phases[:,2]-phases[:,1])/2)]
    F=heat_weights(n,delta)
    return np.array([[2j*np.exp(1j*phi)*F[i,j]*np.exp(1j*(chi[n]-chi[i]+chi[j]))*np.sin(kap[n]-kap[i]-kap[j]) for j in range(n)] for i in range(n)])


def quadratic_fit(K,nmodes):
    cs,_=data_ops(nmodes); D=len(K)
    gam=[x for c in cs for x in (c+c.conj().T,-1j*(c-c.conj().T))]
    scalar=np.trace(K).real/D; fit=scalar*np.eye(D,dtype=complex)
    for a in range(2*nmodes):
        for b in range(a+1,2*nmodes):
            op=1j*gam[a]@gam[b]
            fit+=np.trace(op@K).real/D*op
    inds=[1<<(nmodes-1-j) for j in range(nmodes)]
    one=K[np.ix_(inds,inds)]-K[0,0]*np.eye(nmodes)
    return one,float(np.linalg.norm(K-fit,2))


def parameters(n):
    # Bulk flux pi/2, terminal link flux zero; chi zero.
    ph=np.zeros((n,3));ph[:-1,0]=np.pi/2
    # n=1 has only a terminal link, so give it flux pi/3 for the control.
    if n==1: ph[0,0]=np.pi/3
    return [(3.,ph,np.pi/2,1.),(5.,ph,0.,8.)]


def cell(n,g=.2):
    N=2**n;m=n+1; pars=parameters(n)
    A=-np.diag(np.ones(n),1)-np.diag(np.ones(n),-1)
    ev,uv=eigh(A)
    if uv[0,0]<0: uv=-uv
    parts=[np.kron(A-ev[0]*np.eye(m),np.eye(N))]
    ss=[]
    for d,ph,phi,w in pars:
        h,s=excitation(n,d,ph,g);parts.append(h);ss.append(s)
    H=block_diag(*parts);D=len(H)
    G=np.zeros((D,N),complex);G[:m*N]=np.kron(uv[:,0,None],np.eye(N))
    ports=[[],[]]
    for end,j in enumerate((0,n)):
        V=np.zeros((m*N,N));V[j*N:(j+1)*N]=np.eye(N)
        for r in range(len(pars)):
            fs=[];start=m*N+r*3*m*N
            for a in range(3):
                F=np.zeros((D,D),complex)
                F[start:start+3*m*N,:m*N]=ss[r][end][a]@V.T
                fs.append(F)
            ports[end].append((fs[0]+fs[0].conj().T,fs[1]+fs[2].conj().T,fs[2]+fs[1].conj().T))
    E=np.diag(np.r_[np.zeros(m*N),np.ones(D-m*N)])
    q=np.array([i.bit_count() for i in range(N)])
    charges=np.r_[np.tile(q,m),np.tile(np.concatenate([q+a for j in range(m) for a in (0,1,-1)]),len(pars))]
    return H,G,ports,float(uv[0,0]),E,charges


def combined_bond(n,g=.2,coefficient=False):
    return sum(w*w*bond(n,d,ph,phi,g,coefficient) for d,ph,phi,w in parameters(n))



def target_data(n,g=.2):
    T=sum(w*w*formula(n,d,ph,phi) for d,ph,phi,w in parameters(n))
    mu=np.array([.5,-.7])[:n]*np.abs(np.diag(T)) if n<=2 else np.zeros(n)
    cs,_=data_ops(n)
    W=sum(mu[i]*(cs[i].conj().T@cs[i]) for i in range(n))
    return T,mu,W


def cubic_bond(n,g=.2,eta=10.):
    return -2*eta*sum(w*w*bond(n,d,ph,phi,g,inverse_power=2) for d,ph,phi,w in parameters(n))
