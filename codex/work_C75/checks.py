"""Small exact diagonalization and literal-bond checks. No large jobs."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,time
import numpy as np
from scipy.sparse.linalg import eigsh
import bond18,bond17,bond15
from explore18 import Leading18

def fermions(m):
    out=[]
    for j in range(m):
        c=np.zeros((2**m,2**m))
        for b in range(2**m):
            if (b>>j)&1:c[b^(1<<j),b]=(-1)**((b&((1<<j)-1)).bit_count())
        out.append(c)
    return out
def gaussian_residual(H,n):
    inds=[1<<i for i in range(2*n)];h=H[np.ix_(inds,inds)]-H[0,0]*np.eye(2*n)
    cs=fermions(2*n);ref=H[0,0]*np.eye(len(H),dtype=complex)
    for i in range(2*n):
        for j in range(2*n):ref+=h[i,j]*cs[i].T@cs[j]
    return np.linalg.norm(H-ref,2)/np.linalg.norm(H,2)

def main():
    out={'bond':[],'leading':[],'penalty':[],'ed':[]};start=time.monotonic()
    def save():json.dump(out,open('small_checks.json','w'),indent=2)
    for mod in [bond18,bond17,bond15]:
        p=bond18.params(2,K=100.,g=.071,td=.003,flux=.34,theta=.23,field=.17,musel=.4)
        h=mod.bond(p);q=np.array([bond18.charge(l) for l in mod.LET]);q=(q[:,None]+q[None,:]).ravel()
        r=dict(d=mod.d,hermiticity=float(np.max(abs(h-h.conj().T))),charge_commutator=float(np.max(abs(h*(q[:,None]-q[None,:])))),nonzeros=int(np.count_nonzero(h)))
        assert r['hermiticity']<1e-13 and r['charge_commutator']==0
        out['bond'].append(r);print(r,flush=True)
    for n in [1,2]:
        p=bond18.params(n);z=bond18.pair(p)
        ref=Leading18(n,td=p['td'],theta=p['theta']).matrix(list(range(4**n)))
        r=dict(n=n,cell_dim=len(z['cell']),pair_dim=z['HH'].shape[0],source_relative_error=float(np.linalg.norm(ref-z['H22'])/np.linalg.norm(ref)),gaussian_relative_residual=gaussian_residual(z['H22'],n),code_residual=float(np.max(abs(z['cell']@z['code']-z['E0']*z['code']))))
        assert r['source_relative_error']<1e-11 and r['gaussian_relative_residual']<1e-11
        out['leading'].append(r);print(r,flush=True)
        # Finite-g remainder, before extracting very small hopping entries.
        for g in [.08,.04,.02]:
            y=bond18.pair(dict(p,g=g));rem=np.linalg.norm(y['H2']-y['H20']-g*g*y['H22'],2)
            r=dict(n=n,g=g,remainder_norm=float(rem),remainder_over_g4=float(rem/g**4))
            out.setdefault('conversion',[]).append(r)
    for n in [1,2]:
        p=bond18.params(n,g=.04,D=4.,Delta=.5,td=1e-4)
        ref=bond18.pair(p)
        for K in [100.,1000.]:
            z=bond15.cluster(dict(p,K=K));error=np.linalg.norm(z['H2']-ref['H2'],2)
            legal=np.flatnonzero(z['bad']==0)
            # The first legal states are in original18 order. Exact projection
            # equality checks the quotient prescription, including every term.
            nr=len(ref['HH'].toarray()) if n==1 else ref['HH'].shape[0]
            diff=z['HH'][:nr,:nr]-ref['HH']
            maxdiff=float(max(abs(diff.data),default=0.))
            assert maxdiff<1e-12
            r=dict(n=n,K=K,full_invariant_dim=z['HH'].shape[0],legal_dim=len(legal),source_dim=z['source_dim'],H2_difference_norm=float(error),K_times_error=float(K*error),projected_H0_difference=maxdiff,vac_gap=float(z['vac_gap']))
            out['penalty'].append(r);print(r,flush=True);save()
            if K!=100.:continue
            for lam in [.2,.1,.05]:
                H=z['HH']+lam*z['V'];count=4**n
                vals,vec=eigsh(H,k=count,sigma=2*z['E0']-.25,which='LM',ncv=min(H.shape[0]-1,max(40,5*count)),tol=2e-12,maxiter=8000)
                vals.sort();pred=2*z['E0']+lam*lam*np.linalg.eigvalsh(z['H2'])
                err=float(np.max(abs(vals-pred)))
                # eigsh vector ordering is independent of sorted vals above.
                hv=H@vec;ray=np.sum(vec.conj()*hv,axis=0)
                resid=float(np.max(np.linalg.norm(hv-vec*ray,axis=0)))
                r=dict(n=n,K=K,lambda_=lam,dimension=H.shape[0],retained=count,max_energy_error=err,error_over_lambda4=err/lam**4,eigenpair_residual=resid)
                out['ed'].append(r);print(r,flush=True);save()
    out['seconds']=time.monotonic()-start;save()
if __name__=='__main__':main()
