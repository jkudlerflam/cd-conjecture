"""C80 literal-bond, full-source resolvent, and small ED checks."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,time
import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import eigsh,splu
from scipy.linalg import svd
import bond18 as parent
import bond_shift as b

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

def main():
 out={'bonds':[],'projection':[],'finite_penalty':[],'ed':[]};start=time.monotonic()
 def save():json.dump(out,open('small_checks.json','w'),indent=2)
 for size in [14,13,12]:
  p=parent.params(2,g=.071,td=.003,flux=.34,theta=.23,field=.17,musel=.4,K=100)
  h=b.bond(p,size);q=np.array([b.charge(l) for l in b.alphabet(size)]);qq=(q[:,None]+q[None,:]).ravel()
  rec=dict(alphabet=size,dimension=len(h),nonzeros=int(np.count_nonzero(h)),hermiticity=float(np.max(abs(h-h.conj().T))),charge_commutator=float(np.max(abs(h*(qq[:,None]-qq[None,:])))))
  assert rec['hermiticity']<1e-13 and rec['charge_commutator']==0
  out['bonds'].append(rec)
  for n in [1,2]:
   p=parent.params(n,g=.04,D=4,Delta=.5,td=1e-4,K=100)
   ref=parent.pair(p);z=b.cluster(p,size,closure=False)
   dh=z['HH']-ref['HH'];dv=z['V']-ref['V']
   rec=dict(alphabet=size,n=n,dimension=len(z['states']),H0_difference=float(max(abs(dh.data),default=0)),V_difference=float(max(abs(dv.data),default=0)),leading_gaussian_residual=gaussian_residual(ref['H22'],n))
   assert rec['H0_difference']==0 and rec['V_difference']==0 and rec['leading_gaussian_residual']<1e-12
   out['projection'].append(rec)
 save()
 for n in [1,2]:
  p=parent.params(n,g=.04,D=4,Delta=.5,td=1e-4)
  ref=parent.pair(p)
  for K in [100.,1000.]:
   z=b.cluster(dict(p,K=K),12)
   E,H2,H0aligned,H2aligned,gap,ov=effective(z,ref['G'],n)
   rec=dict(n=n,K=K,dimension=len(z['states']),legal_dimension=z['legal_dim'],zero_port_splitting=float(np.ptp(E)),H0_reference_error=float(np.linalg.norm(H0aligned-2*ref['E0']*np.eye(4**n),2)),H2_reference_error=float(np.linalg.norm(H2aligned-ref['H2'],2)),code_gap=gap,min_overlap_singular_value=ov)
   out['finite_penalty'].append(rec);print(rec,flush=True);save()
   if K!=100:continue
   for lam in [.2,.1,.05]:
    H=z['HH']+lam*z['V'];count=4**n
    vals,vec=eigsh(H,k=count,sigma=float(min(E))-.25,which='LM',ncv=min(H.shape[0]-1,max(40,5*count)),tol=2e-12,maxiter=8000)
    pred=np.linalg.eigvalsh(np.diag(E)+lam*lam*H2)
    err=float(max(abs(np.sort(vals)-pred)))
    resid=float(max(np.linalg.norm(H@vec-vec*vals,axis=0)))
    rec=dict(n=n,K=K,lambda_=lam,dimension=len(z['states']),retained=count,max_energy_error=err,error_over_lambda4=err/lam**4,eigenpair_residual=resid)
    assert resid<1e-9
    out['ed'].append(rec);print(rec,flush=True);save()
 out['seconds']=time.monotonic()-start;save()
if __name__=='__main__':main()
