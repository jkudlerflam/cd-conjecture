"""C44 code-transport, helix spectra and physical-ring gluing checks."""
import itertools,json,time
from pathlib import Path
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
from audit_numeric import core,norm,tl
import core_model
import wrapper_tools as wt

def cell_code_words(n,p,w=.2):
 k=n+2;N=k+p;ch=wt.charge(k,p,w);q=.8
 ww=q**np.arange(k);ww/=la.norm(ww);out=[]
 for bits in range(2**n):
  vec={}
  for ci,holes in enumerate(ch['cfg']):
   for tok in range(k):
    occupied=[('A','v' if tok==0 else None,None)]
    for j in range(n):occupied.append(('d','v',(bits>>j)&1) if tok==j+1 else ('L' if j+1<tok else 'R',None,(bits>>j)&1))
    occupied.append(('B','v' if tok==k-1 else None,None))
    word=[];j=0
    for pos in range(N):
     if pos in holes:word.append(('H',0 if j==0 else (2 if j==k else 1)))
     else:word.append(('occ',0 if j==0 else (2 if j==k-1 else 1),occupied[j]));j+=1
    vec[tuple(word)]=ch['chi'][ci]*ww[tok]
  out.append(vec)
 return out,ch

def transport():
 out=[]
 for n in (1,2):
  a,ch1=cell_code_words(n,1);b,ch2=cell_code_words(n,2);D=2**n;nb=n+3
  bra={}
  for i,j in itertools.product(range(D),repeat=2):
   for x,c in b[i].items():
    for y,d in a[j].items():bra[(x+y,i*D+j)]=c*d
  P=np.zeros((D*D,D*D));wrong=P.copy()
  for i,j in itertools.product(range(D),repeat=2):
   for x,c in a[i].items():
    for y,d in b[j].items():
     word=x+y
     for target in range(D*D):wrong[target,i*D+j]+=bra.get((word,target),0)*c*d
     if word[nb]!=('H',0):continue
     word=word[:nb]+(('H',2),)+word[nb+1:]
     for target in range(D*D):P[target,i*D+j]+=bra.get((word,target),0)*c*d
  beta=wt.transport(ch1,ch2)
  err=norm(P-beta*np.eye(D*D));assert err<1e-12 and norm(wrong)==0
  # Explicit single-site helix tensor factors before/after the local rotation.
  phi=1/nb;plus=lambda t:np.array([1,np.exp(1j*t)])/np.sqrt(2)
  def state(length):return [(plus(j*phi),plus((j-length+1)*phi)) for j in range(length)]
  ket=state(nb)+state(nb+1);target=state(nb+1)+state(nb)
  Ub=np.diag([1,np.exp(1j*nb*phi)]);z=np.exp(-1j*nb*phi)
  for j,((ka,kb),(ba,bb)) in enumerate(zip(ket,target)):
   z*=np.vdot(ba,Ub@ka if j==nb else ka)*np.vdot(bb,Ub@kb if j==nb else kb)
  expected=np.cos(phi/2)**(2*nb)
  assert abs(z-expected)<1e-12
  out.append(dict(n=n,logical_pair_dimension=D*D,charge_core_supports=[sum(map(len,a)),sum(map(len,b))],
    converted_overlap_error=err,unconverted_norm=norm(wrong),beta_vacancy=beta,helix_overlap_error=float(abs(z-expected))))
 return out

def helix_ed():
 out=[]
 for N in (4,5,6):
  sing=np.array([0,1,-1,0])/np.sqrt(2);h=np.outer(sing,sing)
  H=np.kron(np.diag([0,1]),np.eye(2**(N-1)))
  for j in range(N-1):H+=np.kron(np.eye(2**j),np.kron(h,np.eye(2**(N-j-2))))
  e=la.eigvalsh(H);assert abs(e[0])<1e-12 and e[1]>=1/(N*(2*N-1))-1e-12
  out.append(dict(N=N,dimension=len(H),kernel_dimension=int(sum(abs(e)<1e-10)),gap=float(e[1]),lower_bound=1/(N*(2*N-1))))
 return out

def marker_states(b,M,r):
 for longs in itertools.combinations(range(M),r):
  gaps=np.full(M,b,dtype=int);gaps[list(longs)]+=1;rel=np.r_[0,np.cumsum(gaps)[:-1]]
  for first in range(gaps[-1]):yield tuple(map(int,rel+first))

def payload(K,D,M):
 zs=list(itertools.product(range(D),repeat=M));ix={z:i for i,z in enumerate(zs)}
 H=np.zeros((len(zs),len(zs)),complex)
 for i,z in enumerate(zs):
  for j in range(M):
   jj=(j+1)%M
   for x,y in itertools.product(range(D),repeat=2):
    zz=list(z);zz[j]=x;zz[jj]=y
    H[ix[tuple(zz)],i]+=K[x*D+y,z[j]*D+z[jj]]
 return H,zs,ix

def ring_check(n,M,r):
 b=n+3;L=b*M+r;D=2**n;c=core(n);K=core_model.exact_bond(c)
 # Rescale the tiny nonscalar diagnostic bond to visible size. No inference about velocities.
 K=(K+K.conj().T)/2
 K=tl(K);K=K/norm(K)*.3
 Hp,zs,zi=payload(K,D,M);patterns=sorted(marker_states(b,M,r));states=list(itertools.product(patterns,zs));ix={s:i for i,s in enumerate(states)}
 H=np.kron(np.eye(len(patterns)),Hp)
 for i,(p,z) in enumerate(states):
  for j,x in enumerate(p):
   if (x-p[j-1])%L==b and (p[(j+1)%M]-x)%L==b+1:
    pp=list(p);pp[j]=(x+1)%L;pp=tuple(sorted(pp));zz=(z[-1],)+z[:-1] if x==L-1 else z
    a=ix[pp,zz];H[a,i]-=1;H[i,a]-=1
 # Translation eigenfunctions chi(Sz)=exp(2 pi i m/M)chi(z).
 seen=set();orbits=[]
 for z in zs:
  if z in seen:continue
  orb=[];zz=z
  while zz not in seen:seen.add(zz);orb.append(zz);zz=zz[1:]+zz[:1]
  orbits.append(orb)
 pay=[]
 for m in range(M):
  cols=[]
  for orb in orbits:
   if (m*len(orb))%M:continue
   vec=np.zeros(len(zs),complex)
   for j,z in enumerate(orb):vec[zi[z]]=np.exp(2j*np.pi*m*j/M)/np.sqrt(len(orb))
   cols.append(vec)
  if cols:
   B=np.array(cols).T
   for e in la.eigvalsh(B.conj().T@Hp@B):pay.append((e,m))
 levels=[]
 for occ in itertools.combinations(range(M),r):
  q=2*np.pi*(np.asarray(occ)+(r-1)/2)/M;mg=(sum(occ)+r*(r-1)//2)%M
  for ep,mp in pay:
   for ell in range((-mg-mp)%M,L,M):levels.append(ep-2*np.cos(q+2*np.pi*ell/(L*M)).sum())
 ev=la.eigvalsh(H);error=float(max(abs(ev-np.sort(levels))))
 assert len(ev)==len(levels) and error<1e-10
 return dict(n=n,M=M,r=r,L=L,dimension=len(ev),spectrum_error=error,scope='Projected physical marker graph with rescaled exact bare C39 second-order bond, all payload charges. Both marker and payload cyclic ordering retained.')

if __name__=='__main__':
 start=time.monotonic();out=dict(transport=transport(),helix=helix_ed(),rings=[ring_check(1,3,1),ring_check(1,4,2),ring_check(2,3,1)])
 out['seconds']=time.monotonic()-start
 Path('transport_ring_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
