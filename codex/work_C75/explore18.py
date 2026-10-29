"""Two-species leading reduction, independent of the logical Fock size."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import sys
sys.dont_write_bytecode=True
sys.path.insert(0,os.path.abspath('../work_C33'))
from structured import carrier,visibility
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

class Leading18:
    def __init__(self,n,D=10.,a=4.,td=0.,flux=None,loopflux=np.pi/2,theta=0.):
        self.n=n;self.S=2*n+1;self.theta=theta
        if flux is None:flux=np.pi/(2*n)
        self.source_phase=np.exp(-1j*n*flux)
        self.H=carrier(n,D,a,td,0.,loopflux,loopflux)[:2]
        self.lu={(x,y):splu(sp.kronsum(self.H[y],self.H[x],format='csc')) for x in range(2) for y in range(2)}
        b=np.zeros(self.S*self.S,complex);b[(self.S-1)*self.S]=1
        self.rs={xy:self.lu[xy].solve(b).reshape(self.S,self.S) for xy in [(0,1),(1,0)]}
    def Y(self,bits):
        out={};n=self.n
        for x,y,c in [(0,1,self.source_phase*np.exp(1j*self.theta)*(-1)**((bits&((1<<n)-1)).bit_count())),(1,0,np.exp(-1j*self.theta))]:
            R=c*self.rs[x,y]
            for cell,s in enumerate([x,y]):
                local=(bits>>(cell*n))&((1<<n)-1)
                for j in range(n):
                    bit=(local>>j)&1
                    if (s==0 and bit!=1) or (s==1 and bit!=0):continue
                    ns=1-s
                    sign=(-1)**((local&((1<<j)-1)).bit_count())
                    nb=bits^(1<<(cell*n+j))
                    key=(ns,y,nb) if cell==0 else (x,ns,nb)
                    if key not in out:out[key]=np.zeros((self.S,self.S),complex)
                    if cell==0:out[key][j,:]+=sign*R[j,:]
                    else:out[key][:,j]+=sign*R[:,j]
        return out
    def matrix(self,states):
        ys=[self.Y(b) for b in states];out=np.zeros((len(states),len(states)),complex)
        for j,Y in enumerate(ys):
            for key,v in Y.items():
                w=self.lu[key[:2]].solve(v.ravel()).reshape(self.S,self.S)
                for i,X in enumerate(ys):
                    if key in X:out[i,j]-=np.vdot(X[key],w)
        return visibility(self.n)*out
    def calculate(self):
        n=self.n;out=self.matrix([0]+[1<<i for i in range(2*n)])
        h=out[1:,1:]-out[0,0]*np.eye(2*n)
        return h[:n,n:],h[:n,:n]+h[n:,n:]

if __name__=='__main__':
    from fast_loop import fast
    for n in [1,2,3]:
        for theta in [0,np.pi/4]:
            T,M=Leading18(n,td=1e-5,theta=theta).calculate()
            R,N=fast(n,D=10,a=4,td=1e-5)
            print(n,theta,'T ratios',T[0,0]/R[0,0],'M ratio',np.linalg.norm(M)/max(np.linalg.norm(N),1e-100),'T shape error',np.max(abs(T-T[0,0]/R[0,0]*R))/np.max(abs(T)),'Mreal',np.max(abs(M.real)),flush=True)
            print('T',T,'M',M,flush=True)
