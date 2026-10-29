"""Bounded alternative tests. These do not prove a global alphabet lower bound."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import itertools,json
import numpy as np
import bond_shift as b
import bond18 as parent

def quotient_parts(pair,p):
 ls=b.alphabet(12);a,c=pair
 def f(l):return a if l==c else l
 new=list(dict.fromkeys(map(f,ls)));ix={l:i for i,l in enumerate(new)};d=len(new)
 h,_,o,v=b.parts(p,12);out=[];conflicts=0
 for M,arity in [(h,2),(o,1),(v,2)]:
  Q=np.zeros((d**arity,d**arity),complex);seen={}
  for j,i in zip(*np.nonzero(M)):
   def index(k):
    if arity==1:return ix[f(ls[k])]
    return ix[f(ls[k//12])]*d+ix[f(ls[k%12])]
   key=index(j),index(i)
   if key in seen and abs(seen[key]-M[j,i])>1e-12:conflicts+=1
   else:seen[key]=M[j,i]
  for key,val in seen.items():Q[key]=val
  out.append(Q)
 def allowed(x,y):return any(b.allowed(l,r,12) for l in ls for r in ls if f(l)==x and f(r)==y)
 return new,f,out,conflicts,allowed

def main():
 out={'single_letter_merges':[]};p=parent.params(2,D=4,Delta=.5,g=.04,td=1e-4,K=100)
 ls=b.alphabet(12);ws=b.words(2,12)
 for a,c in itertools.combinations(ls,2):
  if b.charge(a)!=b.charge(c):continue
  new,f,parts,conflicts,allowed=quotient_parts((a,c),p)
  images=[tuple(map(f,w)) for w in ws]
  mark=lambda l:l[0] in ['e','m']
  rec=dict(letters=[a,c],legal_cell_images=len(set(images)),parent_legal_dimension=len(ws),rule_coefficient_conflicts=conflicts,merges_marker_with_data=mark(a)!=mark(c))
  if a[0] in ['L','R'] and c[0] in ['L','R']:
   # Full zero-excitation n=2 cell, including zero-penalty extra words.
   d=len(new);ix={l:i for i,l in enumerate(new)};h,o,v=parts
   data=[l for l in new if l[0] in ['L','R']];marks=[l for l in new if l[0]=='e']
   words=[x+(m,) for x in itertools.product(data,repeat=2) for m in marks];idx={w:i for i,w in enumerate(words)}
   fixed=f(('e',0));H=np.zeros((len(words),len(words)),complex)
   for col,w in enumerate(words):
    H[col,col]+=sum(o[ix[l],ix[l]] for l in w)
    full=(fixed,)+w
    for pos in range(3):
     H[col,col]+=p['K']*(not allowed(full[pos],full[pos+1]))
     ci=ix[full[pos]]*d+ix[full[pos+1]]
     for ri in np.flatnonzero(h[:,ci]):
      z=list(full);z[pos]=new[ri//d];z[pos+1]=new[ri%d]
      if z[0]==fixed and tuple(z[1:]) in idx:H[idx[tuple(z[1:])],col]+=h[ri,ci]
   eig=np.linalg.eigvalsh(H)
   rec['zero_excitation_cell_dimension']=len(words)
   rec['lowest_six_energies_at_K100']=eig[:6].tolist()
  out['single_letter_merges'].append(rec)
 # Fold deletion: identical real carriers on two sites.
 H=np.array([[4.,-1],[-1,4.]])
 A=np.kron(H,np.eye(2))+np.kron(np.eye(2),H);R=np.linalg.inv(A)
 x=R[:,2].reshape(2,2);U=[];V=[]
 for j in range(2):
  u=np.zeros((2,2));v=u.copy();u[j,:]=x[j,:];v[:,j]=x[:,j]
  U.append(u.ravel());V.append(v.ravel())
 F=np.array([[u@R@v for v in V] for u in U])
 eig=np.array([np.linalg.eigvalsh(F*np.exp(1j*k)+F.T*np.exp(-1j*k)) for k in np.linspace(-np.pi,np.pi,4001)])
 out['different_fold']=dict(n=2,F=F.tolist(),bc_minus_a_squared=float(F[0,1]*F[1,0]-F[0,0]**2),band_intervals=[[float(min(eig[:,j])),float(max(eig[:,j]))] for j in range(2)],common_field_interval_overlap=float(min(eig.max(axis=0))-max(eig.min(axis=0))))
 # Two noncrossing heads carrying A and B. Their spatial order is invariant.
 L=6;states=[(i,j,s) for i in range(L) for j in range(i+1,L) for s in range(2)];ix={w:i for i,w in enumerate(states)}
 H=np.zeros((len(states),len(states)))
 for col,(i,j,s) in enumerate(states):
  for a,c in [(i-1,j),(i+1,j),(i,j-1),(i,j+1)]:
   if 0<=a<c<L:H[ix[a,c,s],col]=-1
 order=np.diag([(-1)**s for i,j,s in states])
 out['two_token_order']=dict(sites=L,dimension=len(H),order_commutator=float(np.max(abs(H@order-order@H))),cross_order_matrix_element=float(np.max(abs(H[::2,1::2]))),two_order_total_charges=[0,0])
 json.dump(out,open('alternative_checks.json','w'),indent=2)
 print(json.dumps(out),flush=True)
if __name__=='__main__':main()
