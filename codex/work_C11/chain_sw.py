from c11_calculation import *
from itertools import combinations,product

def herm_basis(gam):
 out=[((),np.eye(len(gam[0])))]
 for k in range(1,len(gam)+1):
  for ix in combinations(range(len(gam)),k):
   out.append((ix,1j**(k*(k-1)//2)*reduce(np.matmul,[gam[a] for a in ix])))
 return out

def mixed_response(H,G,L,R):
 ev,U=eigh(H); h=G.shape[1]; e0=np.trace(G.conj().T@H@G).real/h
 gap=ev-e0; w=np.zeros_like(gap); w[gap>1e-8]=1/gap[gap>1e-8]; Res=(U*w)@U.conj().T
 return np.array([[G.conj().T@(a@Res@b+b@Res@a)@G for b in R] for a in L])

def expand_continuum(terms,r):
 # terms = coefficient * ordered product of four canonical lattice Majoranas;
 # use Gamma_{j,a}=sqrt(2a)[R_a+(-1)^j L_a], a=1.
 # Only tree-level normal-ordered dimension-2 terms are returned.
 out={}
 for coeff,ops in terms:
  for lr in product(range(2),repeat=4):
   if sum(lr)%2: continue # rapidly oscillating
   c=4*coeff; inds=[]
   for (site,flav),chir in zip(ops,lr):
    c*=(-1)**(site*chir); inds.append(chir*r+flav)
   if len(set(inds))<4: continue
   inversions=sum(inds[i]>inds[j] for i in range(4) for j in range(i+1,4))
   key=tuple(sorted(inds)); out[key]=out.get(key,0)+c*(-1)**inversions
 return {k:v for k,v in out.items() if abs(v)>1e-10}

for qp,qm in [(1,1),(1,2)]:
 for g in [.03,.1]:
  H,G,L,R,gam,eta,P,ev=double(qp,qm,g); r=len(gam); h=len(P); bas=herm_basis(gam)
  alpha=np.trace(G.conj().T@L[0]@G@gam[0]).real/h
  H2=bond_sw(sectors(H,G,R),sectors(H,G,L)); terms=[]
  print('\nMODEL',qp,qm,'g',g,'alpha',alpha,flush=True)
  for (ix,A),(jx,B) in product(bas,repeat=2):
   coeff=np.trace(np.kron(A,B)@H2).real/h**2
   if abs(coeff)<1e-9: continue
   if len(ix)==len(jx)==2:
    terms.append((-coeff,[(0,a) for a in ix]+[(1,b) for b in jx]))
   elif not ix and not jx: pass
   else: print('bond other',ix,jx,coeff,flush=True)
  bondcont=expand_continuum(terms,r)
  M=mixed_response(H,G,L,R); adj=[]
  for a,b in product(range(r),repeat=2):
   K=P@M[a,b]
   for ix,O in bas:
    c=np.trace(O@K)/h
    if abs(c)<1e-9: continue
    print('adjacent K',a,b,ix,c,'Hamiltonian prefactor',alpha**2*c,flush=True)
    if len(ix)==2:
     adj.append((alpha**2*c.real,[(0,a)]+[(1,c) for c in ix]+[(2,b)]))
  adjcont=expand_continuum(adj,r); total=expand_continuum(terms+adj,r)
  print('bond continuum',bondcont,flush=True)
  print('adj continuum',adjcont,flush=True)
  print('total continuum',total,flush=True)
