"""C80 numerical helpers, copied locally for the C84 checks."""
import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import splu
from scipy.linalg import svd

def gaussian_residual(H,n):
 m=2*n;inds=[1<<j for j in range(m)]
 h=H[np.ix_(inds,inds)]-H[0,0]*np.eye(m)
 cs=[]
 for j in range(m):
  c=np.zeros((2**m,2**m))
  for x in range(2**m):
   if (x>>j)&1:c[x^(1<<j),x]=(-1)**((x&((1<<j)-1)).bit_count())
  cs.append(c)
 ref=H[0,0]*np.eye(len(H),dtype=complex)
 for i in range(m):
  for j in range(m):ref+=h[i,j]*cs[i].T@cs[j]
 return float(np.linalg.norm(H-ref,2)/np.linalg.norm(H,2))

def effective(z,ideal,n):
 H,V=z['HH'],z['V'];nv=4**n
 vac=np.flatnonzero(z['exc']==0);src=np.flatnonzero(z['exc']==2)
 ev,U=np.linalg.eigh(H[vac][:,vac].toarray());gap=ev[nv]-ev[nv-1]
 E=ev[:nv];G=np.zeros((len(z['states']),nv),complex);G[vac]=U[:,:nv]
 S=(V@G)[src];A=H[src][:,src];XS=np.empty_like(S);cache={}
 for j,e in enumerate(E):
  key=round(float(e),11)
  if key not in cache:cache[key]=splu(A-e*eye(len(src),format='csc'))
  XS[:,j]=cache[key].solve(S[:,j])
 B=S.conj().T@XS;H2=-(B+B.conj().T)/2
 Gi=np.zeros_like(G);Gi[:len(ideal)]=ideal
 uu,ss,vh=svd(G.conj().T@Gi);C=uu@vh
 return E,H2,C.conj().T@np.diag(E)@C,C.conj().T@H2@C,float(gap),float(min(ss))

