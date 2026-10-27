import numpy as np
from scipy.linalg import eigh
from functools import reduce
I2=np.eye(2); X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]]); Z=np.diag([1,-1])
def kron(*a): return reduce(np.kron,a)
def ingredients(qp,qm):
 q=qp+qm; r=2*qp*qm; n=r//2; h=2**n
 gam=[]
 for j in range(n):
  for m in (X,Y): gam.append(kron(*([Z]*j+[m]+[I2]*(n-j-1))))
 P=kron(*([Z]*n)); eta=[1j*a@P for a in gam]
 F=[]
 for i in range(qp):
  for j in range(qp,q):
   a=np.zeros((q,q),complex); a[i,j]=a[j,i]=np.sqrt(q/2); F.append(a)
   a=np.zeros((q,q),complex); a[i,j]=-1j*np.sqrt(q/2); a[j,i]=1j*np.sqrt(q/2); F.append(a)
 om=np.eye(q).reshape(-1)/np.sqrt(q); wa=[np.kron(f,np.eye(q))@om for f in F]
 ops=[np.outer(a,om.conj()) for a in wa]
 return q,r,h,gam,eta,P,F,om,ops

def single(qp,qm,g,Delta=1):
 q,r,h,gam,eta,P,F,om,ops=ingredients(qp,qm)
 H=Delta*np.kron(np.eye(q*q)-np.outer(om,om),np.eye(h))
 B=sum(np.kron(a,b) for a,b in zip(ops,gam)); H=H+g*(B+B.conj().T)
 D=np.sqrt(Delta**2+4*r*g*g); theta=.5*np.arctan2(2*g*np.sqrt(r),Delta)
 V=sum(np.kron((a@om)[:,None],b) for a,b in zip(ops,gam))/np.sqrt(r)
 G=np.cos(theta)*np.kron(om[:,None],np.eye(h))-np.sin(theta)*V
 ports=[kron(f,np.eye(q),np.eye(h)) for f in F]
 return H,G,ports,gam,eta,P

def double(qp,qm,g,Delta=1):
 q,r,h,gam,eta,P,F,om,ops=ingredients(qp,qm); A=q*q; IA=np.eye(A); Ih=np.eye(h)
 Q=IA-np.outer(om,om)
 H=Delta*(kron(Q,IA,Ih)+kron(IA,Q,Ih))
 B=sum(kron(a,IA,b)+kron(IA,a,c) for a,b,c in zip(ops,gam,eta)); H=H+g*(B+B.conj().T)
 ev,U=eigh(H); G0=kron(om[:,None],om[:,None],Ih)
 # Polar align eigenvectors with the bare logical basis.
 Ug=U[:,:h]; ov=Ug.conj().T@G0; u,s,vh=np.linalg.svd(ov); G=Ug@u@vh
 left=[kron(f,np.eye(q),IA,Ih) for f in F]; right=[kron(IA,f,np.eye(q),Ih) for f in F]
 return H,G,left,right,gam,eta,P,ev

def sectors(H,G,ports,tol=1e-9):
 h=G.shape[1]; e0=np.trace(G.conj().T@H@G).real/h
 ev,U=eigh(H); groups=[]
 for i,e in enumerate(ev-e0):
  if e<tol: continue
  if not groups or abs(e-groups[-1][0])>tol: groups.append([e,[]])
  groups[-1][1].append(i)
 # A sector is the r x r array of matrices G^dag Fa Psector Fb G.
 Fg=np.stack([f@G for f in ports]); out=[]
 out.append((0., np.einsum('aij,bjk->abik',np.stack([G.conj().T@x for x in Fg]),np.stack([G.conj().T@x for x in Fg]))))
 for e,ix in groups:
  V=U[:,ix]; M=np.einsum('ki,ail->akl',V.conj().T,Fg)
  E=np.einsum('aki,bkj->abij',M.conj(),M)
  out.append((e,E))
 return out

def bond_sw(seca,secb):
 h=seca[0][1].shape[-1]; H=np.zeros((h,h,h,h),complex)
 for e,A in seca:
  for f,B in secb:
   if e+f<1e-10: continue
   H-=np.einsum('abij,abkl->ikjl',A,B)/(e+f)
 return H.reshape(h*h,h*h)

def pauli_coeff(H):
 # Return Clifford bilinear products on a two-cell bond separately below.
 pass

def run():
 for qp,qm in [(1,1),(1,2)]:
  for g in [.03,.1,.3]:
   H,G,F,gam,eta,P=single(qp,qm,g); r=len(F); h=len(P)
   S=sectors(H,G,F); H2=bond_sw(S,S)
   coeff=[]
   for a in range(r):
    for b in range(a+1,r):
     row=[]
     for c in range(r):
      for d in range(c+1,r):
       O=np.kron(gam[a]@gam[b],gam[c]@gam[d]); row.append(np.trace(O.conj().T@H2).real/h**2)
     coeff.append(row)
   print('SINGLE',qp,qm,g,'Hermitian error',np.linalg.norm(H2-H2.conj().T),'quartic matrix',np.array(coeff),flush=True)
   if r==2:
    u=1/np.sqrt(1+4*r*g*g); D=1/u
    k=u*(u-1)**2*(3*u**3+18*u**2+38*u+24)/(D*r*r*(u+1)*(u+3))
    print('analytic kappa',k,'difference',coeff[0][0]-k,flush=True)
  H,G,L,R,gam,eta,P,ev=double(qp,qm,.1)
  h=len(P); print('DOUBLE',qp,qm,'lowest',ev[:h+4], 'ground spread',np.ptp(ev[:h]),flush=True)
  for name,F,target in [('left',L,gam),('right',R,eta)]:
   M=[G.conj().T@f@G for f in F]; alpha=np.mean([np.trace(a@b).real/h for a,b in zip(M,target)])
   print(name,'alpha',alpha,'residual',max(np.linalg.norm(a-alpha*b) for a,b in zip(M,target)),flush=True)
  if np.ptp(ev[:h])<1e-8:
   H2=bond_sw(sectors(H,G,R),sectors(H,G,L))
   print('double bond quartics',[[round(np.trace(np.kron(gam[a]@gam[b],gam[c]@gam[d]).conj().T@H2).real/h**2,12) for c in range(len(gam)) for d in range(c+1,len(gam))] for a in range(len(gam)) for b in range(a+1,len(gam))],flush=True)
if __name__=='__main__': run()
