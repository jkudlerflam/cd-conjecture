# C44 snapshot of ../work_C39/model.py, unchanged implementation.
"""C39 two independent folded buses. Legal-cell matrices in a tree gauge."""
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as sla

def fermions(n):
 N=2**n; out=[]
 for j in range(n):
  c=np.zeros((N,N))
  for b in range(N):
   if b>>j&1: c[b^(1<<j),b]=(-1)**((b&((1<<j)-1)).bit_count())
  out.append(c)
 return out

def scalar(n,w,D,a=4.):
 m=2*n+2; H=np.diag([a+D]+[a]*n+[a+D]*(n+1))
 edges=[(0,1,1.),(1,n+1,1.),(2*n,2*n+1,1.)]
 edges += [(1+j,2+j,w) for j in range(n-1)]
 edges += [(n+1+j,n+2+j,1.) for j in range(n-1)]
 for i,j,t in edges: H[i,j]=H[j,i]=-t
 return H

def excited(n,w,D,g,a=4.):
 N=2**n;m=2*n+2; C=fermions(n)
 H=sp.kron(sp.eye(3),sp.kron(sp.csc_matrix(scalar(n,w,D,a)),sp.eye(N)),format='csc')
 conv=sp.lil_matrix(H.shape,dtype=complex)
 for j,c in enumerate(C):
  pos=1+j
  for s,A in [(1,c),(2,c.T)]:
   r=s*m*N+pos*N; t=pos*N
   conv[r:r+N,t:t+N]=A;conv[t:t+N,r:r+N]=A.T
 return H+g*conv.tocsc(),conv.tocsc()

def sources(n):
 N=2**n;m=2*n+2; P=np.diag([(-1)**b.bit_count() for b in range(N)])
 out={}
 for side,pos in [('L',0),('R',m-1)]:
  out[side]={}
  for s in range(3):
   A=np.zeros((3*m*N,N),complex)
   A[(s*m+pos)*N:(s*m+pos+1)*N]=(-1j*P if side=='R' and s==0 else np.eye(N))
   out[side][s]=A
 return out

def source_pair(n,phi):
 S=sources(n)
 return (np.kron(S['R'][0],S['L'][0])+np.exp(1j*phi)*np.kron(S['R'][1],S['L'][2])+np.exp(-1j*phi)*np.kron(S['R'][2],S['L'][1]))

def k2(n,w,D,g,phi,leading=False):
 H,C=excited(n,w,D,0 if leading else g)
 d=H.shape[0]; A=sp.kronsum(H,H,format='csc'); B=source_pair(n,phi)
 lu=sla.splu(A)
 X=lu.solve(B)
 if leading:
  CC=sp.kronsum(C,C,format='csc'); Y=CC@X
  return -Y.conj().T@lu.solve(Y)
 return -B.conj().T@X

def hopping(H,n):
 # two-cell basis index left*N+right, little-endian within each cell
 N=2**n
 return H[np.ix_([2**j*N for j in range(n)],[2**j for j in range(n)])]

if __name__=='__main__':
 for phi in [np.pi/2,np.pi,0]:
  H=k2(2,.3,20,0,phi,True)
  print(phi,hopping(H,2))


def cell(n,w=.8,v=.65,D=.2,g=.4,q=.8,a=4.):
 N=2**n;k=n+2;m=2*n+2
 Lv=np.zeros((k,k))
 for j in range(k-1):
  Lv[j,j]+=q*q;Lv[j+1,j+1]+=1;Lv[j,j+1]-=q;Lv[j+1,j]-=q
 Hv=sp.kron(sp.csc_matrix(Lv),sp.eye(N),format='csc')
 HA,_=excited(n,w,D,g,a);HB,_=excited(n,w*v,D,g,a)
 H=sp.block_diag([Hv,HA,HB],format='csc');d=H.shape[0]
 weights=q**np.arange(k);weights=weights/np.linalg.norm(weights)
 G=np.zeros((d,N));G[:k*N]=np.kron(weights[:,None],np.eye(N))
 charges=np.array([b.bit_count() for b in range(N)])
 Q=np.r_[np.tile(charges,k),*[np.concatenate([np.tile(charges+s,m) for s in [0,1,-1]]) for _ in range(2)]]
 typ=np.r_[np.zeros(k*N,int),np.ones(3*m*N,int),np.full(3*m*N,2,int)]
 F={};src=sources(n)
 for bus in [1,2]:
  offset=k*N+(bus-1)*3*m*N
  for side,pos in [('L',0),('R',k-1)]:
   vi=np.arange(pos*N,(pos+1)*N)
   maps=[]
   for s in range(3):
    Y=sp.lil_matrix((d,d),dtype=complex)
    Y[offset:offset+3*m*N,vi]=src[side][s]
    maps.append(Y.tocsc())
   F[bus,side,0]=maps[0]+maps[0].getH()
   F[bus,side,1]=maps[1]+maps[2].getH()
   F[bus,side,2]=F[bus,side,1].getH()
 # Conditional vacuum-token data density.
 R=np.zeros(d)
 for j in range(n):
  R[(j+1)*N:(j+2)*N]=[(b>>j)&1 for b in range(N)]
 return dict(n=n,d=d,N=N,H=H,G=G,Q=Q,typ=typ,F=F,R=R,weights=weights,
             params=dict(w=w,v=v,D=D,g=g,q=q,a=a))

def exact_bond(c,eta=.6):
 n=c['n'];p=c['params'];N=c['N'];Hs=[]
 for wr,phi,weight in [(p['w'],np.pi/2,1.),(p['w']*p['v'],np.pi,eta)]:
  H,_=excited(n,wr,p['D'],p['g'],p['a']);e,U=la.eigh(H.toarray())
  S=sources(n);Z={side:{s:U.conj().T@S[side][s] for s in range(3)} for side in ['L','R']}
  B=(np.kron(Z['R'][0],Z['L'][0])+np.exp(1j*phi)*np.kron(Z['R'][1],Z['L'][2])+np.exp(-1j*phi)*np.kron(Z['R'][2],Z['L'][1]))
  denom=(e[:,None]+e[None,:]).reshape(-1)
  Hs.append(-weight*(B.conj().T@(B/denom[:,None])))
 vis=(c['weights'][0]*c['weights'][-1])**2
 return vis*sum(Hs)
