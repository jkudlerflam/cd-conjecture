import numpy as np
from itertools import combinations,product

def leading(qp,qm):
 q=qp+qm;r=2*qp*qm;F=[]
 for i in range(qp):
  for j in range(qp,q):
   a=np.zeros((q,q),complex);a[i,j]=a[j,i]=np.sqrt(q/2);F.append(a)
   a=np.zeros((q,q),complex);a[i,j]=-1j*np.sqrt(q/2);a[j,i]=1j*np.sqrt(q/2);F.append(a)
 pairs=list(combinations(range(r),2)); n=len(pairs);J=np.zeros((r,r,n),complex);L=J.copy()
 for a,b in product(range(r),repeat=2):
  for z,(c,d) in enumerate(pairs):
   J[a,b,z]=int(a==c and b==d)-int(a==d and b==c)
   L[a,b,z]=np.trace(F[c]@F[a]@F[b]@F[d]-F[d]@F[a]@F[b]@F[c])/q
 J=J.reshape(r*r,n);L=L.reshape(r*r,n)
 K=np.real_if_close(119/4*J.T@J-7/4*(J.T@L+L.T@J)-.5*L.T@L).real
 M=np.zeros_like(K);chiral={}; pd={p:i for i,p in enumerate(pairs)}
 for i,A in enumerate(pairs):
  for j,B in enumerate(pairs):
   if abs(K[i,j])<1e-8:continue
   for lr in product(range(2),repeat=4):
    if sum(lr)%2:continue
    c=4*K[i,j]*(-1)**(lr[2]+lr[3]);ix=[fl+r*ch for fl,ch in zip(A+B,lr)]
    if len(set(ix))<4:continue
    c*=(-1)**sum(ix[a]>ix[b] for a in range(4) for b in range(a+1,4));key=tuple(sorted(ix))
    if sum(lr)==2:M[pd[key[:2]],pd[tuple(t-r for t in key[2:])]]+=c
    else:chiral[key]=chiral.get(key,0)+c
 return K,M,chiral,pairs

if __name__=='__main__':
 for m in range(1,7):
  K,M,C,pairs=leading(1,m); pd={p:i for i,p in enumerate(pairs)}; charge=np.zeros(len(pairs))
  for i in range(m):charge[pd[(2*i,2*i+1)]]=1/np.sqrt(m)
  gc=charge@M@charge
  if m>1:
   spin=np.zeros(len(pairs));spin[pd[(0,1)]]=1/np.sqrt(2);spin[pd[(2,3)]]=-1/np.sqrt(2)
   pair=np.zeros(len(pairs));pair[pd[(0,2)]]=1/np.sqrt(2);pair[pd[(1,3)]]=-1/np.sqrt(2)
   gs=spin@M@spin;gp=pair@M@pair
   expected=[gc]+[gs]*(m*m-1)+[gp]*(m*(m-1))
   err=np.max(np.abs(np.linalg.eigvalsh(M)-np.sort(expected)))
   print('m',m,'charge',gc,'SU(m)',gs,'pair',gp,'spectrum residual',err,'chiral sample',C.get((0,1,2,3)),flush=True)
  else:print('m',m,'charge',gc,flush=True)
