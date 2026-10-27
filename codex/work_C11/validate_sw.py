from c11_calculation import *
from scipy.sparse import csr_matrix,eye,kron as skron
from scipy.sparse.linalg import eigsh

def mixed(H,G,L,R):
 ev,U=eigh(H);e0=np.trace(G.conj().T@H@G).real/G.shape[1]; gap=ev-e0;w=np.zeros_like(gap);w[gap>1e-9]=1/gap[gap>1e-9];Res=(U*w)@U.conj().T
 return [[G.conj().T@(a@Res@b+b@Res@a)@G for b in R]for a in L]
def spk(*a):return reduce(skron,a)
H,G,L,R,gam,eta,P,ev=double(1,1,.3);N=len(H);h=len(P);e0=ev[0];I=np.eye(h)
L0=[G.conj().T@x@G for x in L];R0=[G.conj().T@x@G for x in R]
H1=sum(np.kron(a,b) for a,b in zip(R0,L0));H2=bond_sw(sectors(H,G,R),sectors(H,G,L));M=mixed(H,G,L,R)
for nc in [2,3]:
 Id=eye(N,format='csr');Hs=csr_matrix(H-e0*np.eye(N));Ls=list(map(csr_matrix,L));Rs=list(map(csr_matrix,R))
 Hfree=sum(spk(*([Id]*j+[Hs]+[Id]*(nc-j-1)))for j in range(nc))
 V=sum(spk(*([Id]*j+[a,b]+[Id]*(nc-j-2)))for j in range(nc-1)for a,b in zip(Rs,Ls))
 if nc==2: he1,he2=H1,H2
 else:
  he1=kron(H1,I)+kron(I,H1);he2=kron(H2,I)+kron(I,H2)-sum(kron(R0[a],M[a][b],L0[b])for a in range(len(L))for b in range(len(L)))
 for lam in [.004,.002,.001]:
  eig=np.sort(eigsh(Hfree+lam*V,k=h**nc,which='SA',tol=1e-12,return_eigenvectors=False))
  pred=np.linalg.eigvalsh(lam*he1+lam**2*he2)
  err=np.max(np.abs(eig-pred))
  print('cells',nc,'lambda',lam,'max error',err,'error/lambda^3',err/lam**3,flush=True)
