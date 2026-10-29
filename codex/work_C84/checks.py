"""C84 literal bond, independent bit-carrying reference, scalar reduction and ED."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu,eigsh
import bond18 as parent
import bond11 as b
from parent_checks import gaussian_residual,effective
from structured import carrier,visibility

def reference(p):
 sites=p['n'];n=sites-1
 Hp,Gp,fp,exp,wp=parent.cell(p)
 sel=[i for i,w in enumerate(wp) if [l[2] for l in w if l[2] is not None][-1]==0]
 H=Hp[np.ix_(sel,sel)];G=Gp[sel,:2**n];ex=exp[sel]
 fs={key:val[sel][:,sel] for key,val in fp.items()}
 V=np.exp(1j*p['theta'])*sp.kron(fs['right','A'],fs['left','B'])+np.exp(-1j*p['theta'])*sp.kron(fs['right','B'],fs['left','A'])
 V=(V+V.conj().T).tocsc();HH=sp.kronsum(sp.csc_matrix(H),sp.csc_matrix(H),format='csc')
 GG=np.kron(G,G);perm=[a*2**n+c for c in range(2**n) for a in range(2**n)];GG=GG[:,perm]
 E0=p['musel']-2*np.cos(np.pi/(sites+2));src=np.flatnonzero(np.kron(ex,ex));S=(V@GG)[src]
 A=HH[src][:,src]-2*E0*sp.eye(len(src),format='csc')
 H2=-S.conj().T@splu(A).solve(S)
 H0p=parent.cell(dict(p,g=0))[0];H0=H0p[np.ix_(sel,sel)]
 CC=sp.kronsum(sp.csc_matrix((H-H0)/p['g']),sp.csc_matrix((H-H0)/p['g']),format='csc')[src][:,src]
 A0=sp.kronsum(sp.csc_matrix(H0),sp.csc_matrix(H0),format='csc')[src][:,src]-2*E0*sp.eye(len(src),format='csc')
 lu=splu(A0);X=lu.solve(S);Y=CC@X
 H22=-Y.conj().T@lu.solve(Y);H20=-S.conj().T@X
 return dict(HH=HH,V=V,G=GG,H2=H2,H22=H22,H20=H20,E0=E0,A0=A0,CC=CC,S=S)

def scalar(p):
 sites=p['n'];n=sites-1;S=2*sites+1
 a=p['Delta']+2*np.cos(np.pi/(sites+2))
 H=carrier(sites,D=p['D'],a=a,td=p['td'],flux=0,loopflux=np.pi/2)[0]
 lu=splu(sp.kronsum(H,H,format='csc'));src=np.zeros(S*S);src[(S-1)*S]=1
 x=lu.solve(src).reshape(S,S);U=np.zeros((S*S,n),complex);V=U.copy()
 for j in range(n):
  U[j*S:(j+1)*S,j]=x[j,:];V[j::S,j]=x[:,j]
 C=U.conj().T@lu.solve(V);A=U.conj().T@lu.solve(U)+V.conj().T@lu.solve(V)
 return -2*visibility(sites)*C.real,-2j*visibility(sites)*A.imag

def main():
 out=dict(bond={},projection=[],conversion=[],finite_penalty=[],ed=[]);start=time.monotonic()
 def save():json.dump(out,open('small_checks.json','w'),indent=2)
 p=parent.params(3,g=.071,td=.003,flux=.34,theta=.23,field=.17,musel=.4,K=100)
 h=b.bond(p);ls=b.alphabet();q=np.array([b.charge(l) for l in ls]);qq=(q[:,None]+q[None,:]).ravel()
 marks=np.array([b.marker(l) for l in ls]);ml=np.repeat(marks,11);mr=np.tile(marks,11)
 rec=dict(alphabet=11,dimension=len(h),nonzeros=int(np.count_nonzero(h)),hermiticity=float(np.max(abs(h-h.conj().T))),charge_commutator=float(np.max(abs(h*(qq[:,None]-qq[None,:])))),left_marker_commutator=float(np.max(abs(h*(ml[:,None]!=ml[None,:])))),right_marker_commutator=float(np.max(abs(h*(mr[:,None]!=mr[None,:])))))
 assert max(rec[k] for k in ['hermiticity','charge_commutator','left_marker_commutator','right_marker_commutator'])==0
 out['bond']=rec;print(rec,flush=True)
 for n in [1,2]:
  p=parent.params(n+1,g=.04,D=4,Delta=.5,td=1e-4,K=100)
  ref=reference(p);z=b.cluster(p,closure=False)
  dh=z['HH']-ref['HH'];dv=z['V']-ref['V']
  H22=ref['H22'];inds=[1<<j for j in range(2*n)];hm=H22[np.ix_(inds,inds)]-H22[0,0]*np.eye(2*n)
  T,M=scalar(p);sc=max(np.linalg.norm(T),np.linalg.norm(M))
  rec=dict(n=n,data_sites=n+1,dimension=len(z['states']),H0_difference=float(max(abs(dh.data),default=0)),V_difference=float(max(abs(dv.data),default=0)),leading_gaussian_residual=gaussian_residual(H22,n),scalar_T_relative_error=float(np.linalg.norm(hm[:n,n:]-T)/sc),scalar_mu_relative_error=float(np.linalg.norm(hm[:n,:n]+hm[n:,n:]-M)/sc))
  assert rec['H0_difference']==0 and rec['V_difference']==0
  assert max(rec[k] for k in ['leading_gaussian_residual','scalar_T_relative_error','scalar_mu_relative_error'])<1e-10
  out['projection'].append(rec);print(rec,flush=True)
  for g in [.08,.04,.02]:
   hg=-ref['S'].conj().T@splu(ref['A0']+g*ref['CC']).solve(ref['S'])
   err=float(np.linalg.norm(hg-ref['H20']-g*g*H22,2))
   out['conversion'].append(dict(n=n,g=g,remainder_norm=err,remainder_over_g4=err/g**4))
  for K in [100.,1000.]:
   z=b.cluster(dict(p,K=K))
   E,H2,H0a,H2a,gap,ov=effective(z,ref['G'],n)
   rec=dict(n=n,K=K,dimension=len(z['states']),legal_dimension=z['legal_dim'],zero_port_splitting=float(np.ptp(E)),H0_reference_error=float(np.linalg.norm(H0a-2*ref['E0']*np.eye(4**n),2)),H2_reference_error=float(np.linalg.norm(H2a-ref['H2'],2)),code_gap=gap,min_overlap_singular_value=ov)
   out['finite_penalty'].append(rec);print(rec,flush=True);save()
   if K!=100:continue
   for lam in [.2,.1,.05]:
    H=z['HH']+lam*z['V'];count=4**n
    vals,vec=eigsh(H,k=count,sigma=float(min(E))-.25,which='LM',ncv=min(H.shape[0]-1,max(40,5*count)),tol=2e-12,maxiter=8000)
    pred=np.linalg.eigvalsh(np.diag(E)+lam*lam*H2)
    err=float(max(abs(np.sort(vals)-pred)));resid=float(max(np.linalg.norm(H@vec-vec*vals,axis=0)))
    rec=dict(n=n,K=K,lambda_=lam,dimension=len(z['states']),retained=count,max_energy_error=err,error_over_lambda4=err/lam**4,eigenpair_residual=resid)
    assert resid<1e-9
    out['ed'].append(rec);print(rec,flush=True);save()
 out['seconds']=time.monotonic()-start;save()
if __name__=='__main__':main()
