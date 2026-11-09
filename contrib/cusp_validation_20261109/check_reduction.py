"""Independent spin-basis SW check. Run with numpy and scipy, no repo imports.

H=V*D+T+W0/V. D counts adjacent occupied pairs, T flips 00<->11
on NN bonds, and W0 contains NN hopping and a chemical potential.
The low subspace P is ker D. R=Q/D. Canonical SW through physical V^-2:
 K2=P W0 P-P T R T P,
 K3=P T R T R T P-P T R W0 P-P W0 R T P.
All basis signs are spin signs, not Jordan-Wigner conventions.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(key,'2')
import argparse, json, time
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh


def construct(L,w0,mu0,periodic=False):
    dim=1<<L
    states=np.arange(dim,dtype=np.int64)
    nums=np.array([int(s).bit_count() for s in states])
    adj=np.array([(int(s)&(int(s)>>1)).bit_count() + (int(s)&1)*((int(s)>>(L-1))&1)*periodic for s in states])
    p=np.flatnonzero(adj==0)
    tr,tc,tv,wr,wc,wv=[],[],[],[],[],[]
    for s in range(dim):
        wr.append(s); wc.append(s); wv.append(-mu0*nums[s])
        for j in range(L if periodic else L-1):
            k=(j+1)%L
            a=(s>>j)&1; b=(s>>k)&1
            f=s^((1<<j)|(1<<k))
            if a==b:
                tr.append(f);tc.append(s);tv.append(1.)
            else:
                wr.append(f);wc.append(s);wv.append(w0 if a else np.conj(w0))
    T=coo_matrix((tv,(tr,tc)),shape=(dim,dim),dtype=complex).tocsr()
    W=coo_matrix((wv,(wr,wc)),shape=(dim,dim),dtype=complex).tocsr()
    inv=np.zeros(dim);inv[adj>0]=1/adj[adj>0]
    R=diags(inv)
    X=R@T[:,p]
    K2=W[p][:,p]-T[p]@X
    K3=X.conj().T@T@X-X.conj().T@W[:,p]-W[p]@X
    return p,nums,adj,T,W,K2.tocsr(),K3.tocsr()


def pair_formula(L,p,w0,periodic=False):
    idx={int(s):i for i,s in enumerate(p)}
    r,c,v=[],[],[]
    for col,s0 in enumerate(p):
        s=int(s0)
        for distance,amp in ((2,-2*np.real(w0)),(3,2/3)):
            for j in range(L if periodic else L-distance):
                k=(j+distance)%L
                if ((s>>j)&1)==((s>>k)&1):
                    f=s^((1<<j)|(1<<k))
                    if f in idx:
                        r.append(idx[f]);c.append(col);v.append(amp)
    return coo_matrix((v,(r,c)),shape=(len(p),len(p)),dtype=complex).tocsr()


def low_levels(H,k=20):
    if H.shape[0]<400:
        return np.linalg.eigvalsh(H.toarray())[:k]
    vals=eigsh(H,k=k,which='SA',tol=1e-11,return_eigenvectors=False,
               v0=np.random.default_rng(17).normal(size=H.shape[0]),maxiter=200000)
    return np.sort(vals)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--max-L',type=int,default=10)
    parser.add_argument('--spectrum-L',type=int,default=10)
    parser.add_argument('--output',default='results.json');args=parser.parse_args()
    t0=time.time();out={'conventions':'spin chain H=V D+T+W0/V; Delta=1; matrix checks open and periodic, spectra periodic','matrix_checks':[],'spectra':[]}
    for periodic in (False,True):
        for L in range(7 if periodic else 4,args.max_L+1):
            for w0 in (0.5j,0.15+0.35j,0.5):
                p,num,adj,T,W,K2,K3=construct(L,w0,1.9,periodic)
                delta=K3-pair_formula(L,p,w0,periodic)
                err=float(np.max(np.abs(delta.data))) if delta.nnz else 0.
                change=num[p][:,None]-num[p][None,:]
                k3dense=K3.toarray()
                forbidden=float(np.max(np.abs(k3dense)[np.abs(change)!=2]))
                rec=dict(L=L,periodic=periodic,w0=[float(np.real(w0)),float(np.imag(w0))],P_dimension=len(p),
                         formula_max_error=err,number_change_not_two_max=forbidden)
                out['matrix_checks'].append(rec)
                assert err<1e-12 and forbidden<1e-12,rec
    L=args.spectrum_L
    for mu_eff,label in ((-2.2,'actual_candidate'),(-11.5,'C173_old_vacuum_check')):
        p,num,adj,T,W,K2,K3=construct(L,0.5j,(mu_eff+6)/2,periodic=True)
        for V in (10.,20.,40.,80.):
            H=diags(V*adj)+T+W/V
            ef=low_levels(H)
            e2=low_levels(K2/V)
            e3=low_levels(K2/V+K3/V**2)
            energy,vec=eigsh(H,k=1,which='SA',tol=1e-11,
                            v0=np.random.default_rng(3).normal(size=H.shape[0]),maxiter=200000)
            prob=np.abs(vec[:,0])**2
            rec=dict(L=L,label=label,mu_eff=mu_eff,V=V,t2_abs=1/(2*V),
                     full_ground_particles=float(prob@num),full_ground_adjacent_pairs=float(prob@adj),
                     max_low20_error_H2=float(np.max(np.abs(ef-e2))),
                     max_low20_error_H2_H3=float(np.max(np.abs(ef-e3))),
                     ground_error_H2=float(ef[0]-e2[0]),ground_error_H2_H3=float(ef[0]-e3[0]),
                     full_low20=ef.tolist(),H2_low20=e2.tolist(),H2_H3_low20=e3.tolist())
            out['spectra'].append(rec)
            print({k:v for k,v in rec.items() if 'low20' not in k},flush=True)
    out['elapsed_seconds']=time.time()-t0
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
    print('PASS',len(out['matrix_checks']),'exact-coefficient matrix checks; seconds',out['elapsed_seconds'])
if __name__=='__main__':main()
