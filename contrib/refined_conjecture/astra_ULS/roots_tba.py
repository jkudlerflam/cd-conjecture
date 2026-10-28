from tba_check import TBA,e0
from scipy.optimize import brentq
import numpy as np,sys
N=int(sys.argv[1]);M=int(sys.argv[2]);lo=float(sys.argv[3]);hi=float(sys.argv[4]);target=float(sys.argv[5])
K=int(sys.argv[6]) if len(sys.argv)>6 else 1024
obj=TBA(N,M,K,R=max(24,6*N));x=None

def calc(T):
 global x
 out,x=obj.solve(T,x,tol=2e-9)
 print(out,flush=True)
 return out['ratio']-target
root=brentq(calc,lo,hi,xtol=2e-6)
s=target*N*(N-1)*root/6
print('RESULT',dict(N=N,M=M,K=K,v=2*np.pi/N,e0=e0(N),T=root,t=N*N*root,ell=2*np.pi/(N*root),sigma=s/np.log(N)),flush=True)
