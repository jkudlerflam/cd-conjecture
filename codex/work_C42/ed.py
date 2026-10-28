"""Untruncated legal-sector two/three-cell ED, exact symmetry restriction."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import itertools,json,time,signal,argparse
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import scipy.linalg as la
import model

def sector(c,M,Q):
 d=c['d'];groups={(t,q):np.flatnonzero((c['typ']==t)&(c['Q']==q)) for t in range(3) for q in range(-1,c['n']+2)}
 chunks=[]
 for ts in itertools.product(range(3),repeat=M):
  if ts.count(1)%2 or ts.count(2)%2:continue
  for qs in itertools.product(range(-1,c['n']+2),repeat=M):
   if sum(qs)!=Q:continue
   gg=[groups[t,q] for t,q in zip(ts,qs)]
   if any(len(x)==0 for x in gg):continue
   ar=np.zeros(tuple(len(x) for x in gg),dtype=np.int64)
   for j,g in enumerate(gg):
    shape=[1]*M;shape[j]=len(g);ar+=g.reshape(shape)*d**(M-1-j)
   chunks.append(ar.ravel())
 ids=np.sort(np.concatenate(chunks));coord=np.array([(ids//d**(M-1-j))%d for j in range(M)]).T
 return ids,coord

def build(c,M,Q,eta=.6):
 d=c['d'];ids,co=sector(c,M,Q);dim=len(ids);rows=[];cols=[];vals=[]
 def add(old,new,value):
  jj=np.searchsorted(ids,new)
  assert np.all(jj<len(ids)) and np.all(ids[jj]==new)
  rows.append(jj);cols.append(old);vals.append(np.broadcast_to(value,old.shape))
 H=c['H'].tocoo()
 for site in range(M):
  stride=d**(M-1-site)
  index=[np.flatnonzero(co[:,site]==k) for k in range(d)]
  for i,j,z in zip(H.row,H.col,H.data):
   old=index[j]
   if len(old):add(old,ids[old]+(i-j)*stride,z)
 H0=sp.coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(dim,dim)).tocsr()
 rows=[];cols=[];vals=[]
 for site in range(M-1):
  for bus,phi,weight in [(1,np.pi/2,1),(2,np.pi,np.sqrt(eta))]:
   for sr,sl,amp in [(0,0,1),(1,2,np.exp(1j*phi)),(2,1,np.exp(-1j*phi))]:
    maps=[]
    for side,s in [('R',sr),('L',sl)]:
     f=c['F'][bus,side,s].tocoo();dest=np.full(d,-1);value=np.zeros(d,complex)
     assert len(set(f.col))==len(f.col)
     dest[f.col]=f.row;value[f.col]=f.data;maps.append((dest,value))
    a,za=maps[0];b,zb=maps[1];x=co[:,site];y=co[:,site+1]
    old=np.flatnonzero((a[x]>=0)&(b[y]>=0));new=ids[old]+(a[x[old]]-x[old])*d**(M-1-site)+(b[y[old]]-y[old])*d**(M-2-site)
    add(old,new,weight*amp*za[x[old]]*zb[y[old]])
 V=sp.coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(dim,dim)).tocsr()
 N=c['N'];codebits=[b for b in range(N**M) if b.bit_count()==Q];lookup={b:j for j,b in enumerate(codebits)}
 P=np.zeros((dim,len(codebits)),complex)
 for r in np.flatnonzero(np.all(c['typ'][co]==0,axis=1)):
  bits=0;weight=1
  for j in range(M):
   local=co[r,j];bits=bits*N+local%N;weight*=c['weights'][local//N]
  P[r,lookup[bits]]=weight
 R=c['R'][co].sum(axis=1)
 return H0,V,R,P,codebits

def prediction(c,K,M,bits,field):
 N=c['N'];full=np.zeros((N**M,N**M),complex)
 for j in range(M-1):full+=np.kron(np.kron(np.eye(N**j),K),np.eye(N**(M-2-j)))
 R=c['G'].T@(c['R'][:,None]*c['G'])
 for j in range(M):full+=field*np.kron(np.kron(np.eye(N**j),R),np.eye(N**(M-1-j)))
 return full[np.ix_(bits,bits)]

def run(n,M,Q,lams):
 start=time.monotonic();c=model.cell(n,a=3.5);K=model.exact_bond(c);H0,V,R,P,bits=build(c,M,Q)
 pred=prediction(c,K,M,bits,.03);out=dict(n=n,cells=M,Q=Q,dimension=H0.shape[0],code_dimension=P.shape[1],source_dimension=3*(2*n+2)*2**n,results=[])
 out['code_residual']=float(la.norm(H0@P));out['first_order']=float(la.norm(P.conj().T@(V@P)))
 out['hermiticity']=float(sla.norm(V-V.getH()))
 warm=P.copy()
 for lam in lams:
  H=H0+lam*V+sp.diags(lam**2*.03*R)
  k=P.shape[1]+1
  if H.shape[0]>100000:
   pre=sp.diags(1/np.maximum(H.diagonal().real,.2))
   ev,U=sla.lobpcg(H,warm,M=pre,largest=False,tol=8e-11,maxiter=100)
   k=P.shape[1]+1
  else:
   ev,U=sla.eigsh(H,k=k,which='SA',tol=2e-12,ncv=max(30,3*k),maxiter=20000)
  order=np.argsort(ev);ev=ev[order];U=U[:,order];low=U[:,:k-1];warm=low.copy()
  residual=np.max(la.norm(H@U-U*ev,axis=0))
  # Rayleigh-Ritz on the computed subspace.
  low,_=la.qr(low,mode='economic');e,rot=la.eigh(low.conj().T@(H@low));low=low@rot
  X,s,Yh=la.svd(P.conj().T@low);polar=X@Yh;eff=(polar*e)@polar.conj().T
  error=la.norm(eff-lam**2*pred,2)
  off=eff-lam**2*pred;np.fill_diagonal(off,0)
  out['results'].append(dict(offdiag_error=float(la.norm(off,2)),lam=lam,error=float(error),error_over_lam4=float(error/lam**4),residual=float(residual),gap=float(ev[-1]-ev[k-2]) if len(ev)==k else None,minimum_overlap=float(min(s))))
 out['seconds']=time.monotonic()-start
 with open(f'ed_n{n}_m{M}_q{Q}.json','w') as f:json.dump(out,f,indent=2)
 print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('n',type=int);p.add_argument('M',type=int);p.add_argument('Q',type=int);p.add_argument('--seconds',type=int,default=180);args=p.parse_args()
 signal.signal(signal.SIGALRM,lambda *x:(_ for _ in ()).throw(TimeoutError('per-run time budget')));signal.alarm(args.seconds)
 run(args.n,args.M,args.Q,[.1,.05])
