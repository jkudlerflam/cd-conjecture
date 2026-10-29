# C82 local copy of C75 scalar routines for the NEW n=16 certificate.
# C75 reuse of C33 scalar-certificate implementation. See PROOF.md.
"""Polynomial-size exact coefficients of lambda^2 g^2, NumPy/SciPy only.
The finite-g interaction is distinguished from this Taylor coefficient.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from scipy.linalg import eigvalsh

def carrier(n,D=10.,a=4.,td=0.,flux=None,loopflux=.9,chargedloopflux=0.):
    if flux is None: flux=np.pi/(2*n)
    S=2*n+1
    H=[]
    for typ in range(3):
        A=np.diag(np.r_[np.full(n,a),np.full(n+1,a+D)]).astype(complex)
        edges=[(j,j+1,-1.) for j in range(n-1)]
        edges += [(n+j,n+j+1,-np.exp(1j*(flux if typ==0 else 0.))) for j in range(n)]
        edges += [(0,n,-1.)]
        if td:edges += [(n-1,2*n,-td*np.exp(1j*(loopflux if typ==0 else chargedloopflux)))]
        for i,j,v in edges:A[j,i]+=v;A[i,j]+=np.conj(v)
        H.append(sp.csc_matrix(A))
    return H

def visibility(n):
    return (2/(n+2))**2*np.sin(np.pi/(n+2))**4

def tree(n,D=10.,a=4.,normalized=False):
    """Return F_ij=int f_M,i f_0,j. Exactly specified g^2 coefficient.
    Sparse positive solves avoid exponentially small spectral subtractions.
    """
    H=carrier(n,D,a,flux=0)[0].real
    S=H.shape[0];A=sp.kronsum(H,H,format='csc');lu=splu(A)
    b=np.zeros(S*S);b[(S-1)*S]=1
    C=H[n:,n:].toarray()
    logbeta=-np.linalg.slogdet(C)[1]
    x=lu.solve(b).reshape(S,S)
    if normalized:x/=np.exp(logbeta)
    V=np.zeros((S*S,n))
    for j in range(n):V[j::S,j]=x[:,j]
    RV=lu.solve(V).reshape(S,S,n)
    F=np.array([2*x[i,:]@RV[i,:,:] for i in range(n)])
    C=H[n:,n:].toarray()
    logbeta=-np.linalg.slogdet(C)[1]
    return F,logbeta

class Leading:
    """Two-carrier resolvent and one conversion on each source vector.
    No logical Fock-space matrix is built. State keys only store reached bits.
    """
    def __init__(self,n,D=10.,a=4.,td=0.,flux=None,loopflux=.9,bphi=np.pi/2,chargedloopflux=0.):
        self.n=n;self.S=2*n+1;self.bphi=bphi
        if flux is None:flux=np.pi/(2*n)
        self.source_phase=np.exp(-1j*n*flux)
        if abs(self.source_phase.real)<1e-15:self.source_phase=1j*self.source_phase.imag
        self.H=carrier(n,D,a,td,0.,loopflux-n*flux,chargedloopflux)
        self.lu={}
        for x in range(3):
            for y in range(3):
                self.lu[x,y]=splu((sp.kron(self.H[x],sp.eye(self.S))+sp.kron(sp.eye(self.S),self.H[y])).tocsc())
        self.rsource={}
        b=np.zeros(self.S*self.S,complex);b[(self.S-1)*self.S]=1
        for xy in [(0,0),(1,2),(2,1)]:self.rsource[xy]=self.lu[xy].solve(b).reshape(self.S,self.S)
    def Y(self,bits):
        out={};n=self.n
        for x,y,c in [(0,0,self.source_phase*(-1)**((bits&((1<<n)-1)).bit_count())),(1,2,np.exp(1j*self.bphi)),(2,1,np.exp(-1j*self.bphi))]:
            R=c*self.rsource[x,y]
            for cell,s in enumerate([x,y]):
                local=(bits>>(cell*n))&((1<<n)-1)
                for j in range(n):
                    bit=(local>>j)&1
                    if s==0:ns=1 if bit else 2
                    elif s==1:
                        if bit:continue
                        ns=0
                    else:
                        if not bit:continue
                        ns=0
                    sign=(-1)**((local&((1<<j)-1)).bit_count())
                    nb=bits^(1<<(cell*n+j))
                    key=(ns,y,nb) if cell==0 else (x,ns,nb)
                    if key not in out:out[key]=np.zeros((self.S,self.S),complex)
                    if cell==0:out[key][j,:]+=sign*R[j,:]
                    else:out[key][:,j]+=sign*R[:,j]
        return out
    def calculate(self):
        n=self.n;bits=[0]+[1<<i for i in range(2*n)]
        Ys=[self.Y(b) for b in bits]
        # Solve each distinct source block, then contract all compatible keys.
        out=np.zeros((2*n+1,2*n+1),complex)
        for j,Y in enumerate(Ys):
            for key,v in Y.items():
                w=self.lu[key[:2]].solve(v.ravel()).reshape(self.S,self.S)
                for i,X in enumerate(Ys):
                    if key in X:out[i,j]-=np.vdot(X[key],w)
        h=out[1:,1:]-out[0,0]*np.eye(2*n)
        return visibility(n)*h[:n,n:],visibility(n)*(h[:n,:n]+h[n:,n:]),float(np.max(abs(out-out.conj().T)))

if __name__=='__main__':
    for n in [1,2,3]:
        F,lb=tree(n)
        T,M,err=Leading(n).calculate()
        ref=-2*visibility(n)*F
        print(n,'tree_error',np.max(abs(T-ref))/np.max(abs(ref)),'M/T',np.max(abs(M))/np.max(abs(T)),'hermiticity',err,flush=True)
