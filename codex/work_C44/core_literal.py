# C44 snapshot of ../work_C39/literal.py. Modified to return individual
# gauged bond contributions and endpoint terms for the wrapper calculation.
"""Independent physical-letter implementation of the 46-letter core."""
import itertools,json
import numpy as np
import scipy.sparse as sp
import scipy.linalg as la
import core_model as model

def make(n,w=.8,v=.65,D=.2,g=.4,q=.8,a=3.5,K=10.):
 species=['v']+[(r,s) for r in [1,2] for s in [0,1,-1]]
 EX=species[1:]
 letters=[(z,None,b) for z in ['L','R'] for b in [0,1]]
 letters += [('d',s,b) for s in species for b in [0,1]]
 letters += [('x',s,b) for s in EX for b in [0,1]]
 letters += [(end,s,None) for end in ['A','B'] for s in [None]+species]
 assert len(letters)==46
 ix={z:i for i,z in enumerate(letters)};d=46
 Q=np.array([int(s in EX) for z,s,b in letters]);high=np.array([int(s in EX and z in ['A','B','x']) for z,s,b in letters])
 charge=np.array([(b or 0)+(s[1] if s in EX else 0) for z,s,b in letters])
 onsite=sp.diags((a-3)*Q+D*high,format='lil',dtype=complex)
 for r in [1,2]:
  for z1,z2 in [(('d',(r,0),1),('d',(r,1),0)),(('d',(r,0),0),('d',(r,-1),1))]:
   i,j=ix[z1],ix[z2];onsite[i,j]=onsite[j,i]=g
 onsite=onsite.tocsr()
 bond=sp.diags(np.repeat(Q,d)+2*np.tile(Q,d),format='lil',dtype=complex)
 def edge(aa,bb,z,vac=False):
  i=ix[aa[0]]*d+ix[aa[1]];j=ix[bb[0]]*d+ix[bb[1]]
  if vac:bond[i,i]+=q*q;bond[j,j]+=1;z=-q
  bond[j,i]+=z;bond[i,j]+=np.conj(z)
 for b,c in itertools.product([0,1],repeat=2):
  edge((('d','v',b),('R',None,c)),(('L',None,b),('d','v',c)),-q,True)
  for r,s in EX:
   edge((('d',(r,s),b),('R',None,c)),(('L',None,b),('d',(r,s),c)),-w*(v if r==2 else 1)*((-1)**b if s==0 else 1))
   edge((('x',(r,s),b),('R',None,c)),(('L',None,b),('x',(r,s),c)),-np.exp(1j*np.pi/(2*n)*(s==0))*((-1)**b if s==0 else 1))
 for b in [0,1]:
  edge((('A','v',None),('R',None,b)),(('A',None,None),('d','v',b)),-q,True)
  edge((('d','v',b),('B',None,None)),(('L',None,b),('B','v',None)),-q,True)
  for r,s in EX:
   edge((('A',(r,s),None),('R',None,b)),(('A',None,None),('d',(r,s),b)),-1)
   edge((('x',(r,s),b),('B',None,None)),(('L',None,b),('B',(r,s),None)),-np.exp(1j*np.pi/(2*n)*(s==0))*((-1)**b if s==0 else 1))
   for am in [None]+species:
    edge((('A',am,None),('d',(r,s),b)),(('A',am,None),('x',(r,s),b)),-1)
 kinetic=bond.tocsr()
 bond+=sp.kron(onsite,sp.eye(d),format='lil')/2+sp.kron(sp.eye(d),onsite,format='lil')/2
 def incoming(x):
  z,s,b=x
  if z=='B':return 'R' if s is None else 'L'
  return 'R' if z=='R' else 'L'
 def outgoing(x):
  z,s,b=x
  if z=='A':return 'L' if s is None else 'R'
  return 'L' if z=='L' else 'R'
 bad=np.array([r[0]=='A' or l[0]=='B' or outgoing(l)!=incoming(r) for l in letters for r in letters])
 bond+=sp.diags(K*bad,format='lil')
 left=sp.diags(2*Q+K*np.array([z!='A' for z,s,b in letters]),format='csr')+onsite/2
 right=sp.diags(Q+K*np.array([z!='B' for z,s,b in letters]),format='csr')+onsite/2
 # Legal words in the independent model's ordering.
 words=[];gauges=[];N=2**n;k=n+2;m=2*n+2
 sectors=[('v',p) for p in range(k)]
 sectors += [((r,s),p) for r in [1,2] for s in [0,1,-1] for p in range(m)]
 for s,p in sectors:
  if s=='v':kind='A' if p==0 else ('B' if p==k-1 else 'd');site=p
  elif p==0:kind='A';site=0
  elif p==m-1:kind='B';site=k-1
  elif p<=n:kind='d';site=p
  else:kind='x';site=p-n
  for bits in range(N):
   ww=[('A',s if kind=='A' else None,None)]
   for j in range(n):
    b=bits>>j&1
    ww.append((kind,s,b) if site==j+1 else ('L' if j+1<site else 'R',None,b))
   ww.append(('B',s if kind=='B' else None,None));words.append(tuple(ww))
   phase=1
   if s in EX and s[1]==0:
    prefix=0 if kind=='A' else (n if kind=='B' else site-1)
    phase=(-1)**((bits&((1<<prefix)-1)).bit_count())
    if kind in ['x','B']:phase*=np.exp(1j*np.pi/(2*n)*prefix)
   gauges.append(phase)
 wi={word:i for i,word in enumerate(words)};dc=len(words);H=sp.lil_matrix((dc,dc),dtype=complex)
 parts=[sp.lil_matrix((dc,dc),dtype=complex) for _ in range(k-1)]
 endpart=sp.lil_matrix((dc,dc),dtype=complex)
 pair=bond.tocsc();ends=[left.tocsc(),right.tocsc()]
 for col,word in enumerate(words):
  for pos in range(k-1):
   idx=ix[word[pos]]*d+ix[word[pos+1]]
   for z in range(pair.indptr[idx],pair.indptr[idx+1]):
    rr=pair.indices[z];ww=word[:pos]+(letters[rr//d],letters[rr%d])+word[pos+2:]
    assert ww in wi,'legal leakage'
    H[wi[ww],col]+=pair.data[z]
    parts[pos][wi[ww],col]+=pair.data[z]
  for pos,op in [(0,ends[0]),(k-1,ends[1])]:
   idx=ix[word[pos]]
   for z in range(op.indptr[idx],op.indptr[idx+1]):
    ww=word[:pos]+(letters[op.indices[z]],)+word[pos+1:];assert ww in wi
    H[wi[ww],col]+=op.data[z]
    endpart[wi[ww],col]+=op.data[z]
 gauge=np.array(gauges);Hg=sp.diags(gauge.conj())@H.tocsr()@sp.diags(gauge)
 ref=model.cell(n,w,v,D,g,q,a)
 err=la.norm((Hg-ref['H']).toarray(),np.inf)
 b=bond.tocoo();comm=max(abs((charge[b.row//d]+charge[b.row%d]-charge[b.col//d]-charge[b.col%d])*b.data),default=0)
 # Gershgorin bound of the full positive kinetic piece.
 kin=kinetic.tocsr();lower=np.min((kin.diagonal().real-(abs(kin).sum(axis=1).A.ravel()-abs(kin.diagonal())))[np.repeat(Q,d)+np.tile(Q,d)>0])
 out=dict(n=n,alphabet=d,bond_dimension=d*d,core_dimension=dc,entry_error=err,charge_commutator=float(comm),kinetic_gershgorin_lower=float(lower),code_residual=float(la.norm(ref['H']@ref['G'])))
 assert err<1e-12 and comm==0 and lower>=-1e-12
 ref['parts']=np.array([(sp.diags(gauge.conj())@p.tocsr()@sp.diags(gauge)).toarray() for p in parts])
 ref['ends']=(sp.diags(gauge.conj())@endpart.tocsr()@sp.diags(gauge)).toarray()
 ref['literal_check']=out
 return ref,bond.tocsr(),left,right
if __name__=='__main__':
 result=[]
 for n in [1,2,3]:
  ref,b,l,r=make(n);out=ref['literal_check'];print(out);result.append(out)
 with open('literal_checks.json','w') as f:json.dump(result,f,indent=2)
