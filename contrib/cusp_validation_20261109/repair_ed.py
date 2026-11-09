"""Small-ring diagnostic of the complete facilitated leading model.
No central-charge estimate. All particle numbers are retained, vacuum removed
only when computing the lowest nonvacuum level. Units |t2|=1.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(key,'2')
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix,diags,eye
from check_reduction import construct
from check_facilitated import facilitated
from scipy.sparse.linalg import eigsh


def model(L,r,mu,extra_real_hop=0.):
    p=[s for s in range(1,1<<L) if not s&(s>>1) and not ((s&1) and ((s>>(L-1))&1))]
    ix={s:i for i,s in enumerate(p)};rows=[];cols=[];vals=[]
    N=np.array([s.bit_count() for s in p]);v2=-2+4*r*r
    for col,s in enumerate(p):
        diag=-mu*N[col]
        for j in range(L):
            diag+=v2*((s>>j)&1)*((s>>((j+2)%L))&1)
            diag-=2/3*((s>>j)&1)*((s>>((j+3)%L))&1)
            for distance,amp in ((1,1j+extra_real_hop),(2,-1.)):
                k=(j+distance)%L;a=(s>>j)&1;b=(s>>k)&1
                if a!=b:
                    f=s^((1<<j)|(1<<k))
                    if f in ix:
                        rows.append(ix[f]);cols.append(col);vals.append(amp if a else np.conj(amp))
            f=s^(1<<j)
            if f in ix:
                a=((s>>((j+2)%L))&1)-((s>>((j-2)%L))&1)
                if a:
                    rows.append(ix[f]);cols.append(col);vals.append(2*r*a)
        rows.append(col);cols.append(col);vals.append(diag)
    H=coo_matrix((vals,(rows,cols)),shape=(len(p),len(p)),dtype=complex).tocsr()
    assert np.max(np.abs((H-H.conj().T).data),initial=0)<1e-12
    return p,N,H


def main():
    a=argparse.ArgumentParser();a.add_argument('--output',default='repair_ed_results.json')
    a.add_argument('--max-L',type=int,default=20);a.add_argument('--refined',action='store_true');a.add_argument('--only-max-L',action='store_true');a.add_argument('--q-coefficient',type=float);a.add_argument('--r-values',nargs='+',type=float);a.add_argument('--mu-values',nargs='+',type=float);args=a.parse_args();out={'scope':'All-N constrained leading model, finite rings only. Vacuum energy set to zero.','rows':[]}
    start=time.time()
    out['matrix_checks']=[]
    check_qs=(0.,args.q_coefficient if args.q_coefficient is not None else -.4746) if args.refined else (0.,)
    for L in (7,8,10):
      for r in (.2,.3):
        for q in check_qs:
            extra=q*r*r
            p,N,H=model(L,r,-2.2,extra)
            p0,nums,adj,T,W,_,_=construct(L,-r*r+extra/2+0.5j,3-2*r*r-2.2/2,periodic=True)
            F=facilitated(L);A=T+r*F;inv=np.zeros(len(adj));inv[adj>0]=1/adj[adj>0]
            K=2*(W[p0][:,p0]-A[p0]@diags(inv)@A[:,p0])+2*L*eye(len(p0))
            diff=H-K[1:,1:]
            err=float(np.max(np.abs(diff.data))) if diff.nnz else 0.
            assert p==list(p0[1:]) and err<1e-12
            out['matrix_checks'].append(dict(L=L,r=r,q_coefficient=q,max_error=err))
    mus=args.mu_values or ((-2.30,-2.27) if args.refined else (-2.27,-2.2))
    rs=args.r_values or ((.05,.1,.15) if args.refined else (0.,.2,.25,.3,.4,.5))
    qs=(args.q_coefficient,) if args.q_coefficient is not None else ((0.,-.4746) if args.refined else (0.,))
    lengths=(args.max_L,) if args.only_max_L else dict.fromkeys((12,16,18,args.max_L))
    params=[(L,mu,r,q) for L in lengths for mu in mus for r in rs for q in qs]
    for L,mu,r,q in params:
            p,N,H=model(L,r,mu,extra_real_hop=q*r*r)
            t=time.time();e,vs=eigsh(H,k=3,which='SA',tol=1e-10,maxiter=20000,
                                  v0=np.random.default_rng(42).normal(size=len(p)))
            order=np.argsort(e);e=e[order];v=vs[:,order[0]];prob=np.abs(v)**2
            nr=float(prob@N);nv=float(prob@(N*N)-nr*nr)
            rec=dict(L=L,r=r,mu=mu,q_coefficient=q,extra_real_hop=q*r*r,dim=len(p),nonvacuum_levels=e.tolist(),
                     ground_E=float(min(0,e[0])),excited_state_N=nr,excited_state_varN=nv,
                     physical_ground_N=nr if e[0]<0 else 0.,seconds=time.time()-t)
            lookup={s:i for i,s in enumerate(p)}
            translation=np.array([lookup[((s<<1)&((1<<L)-1))|(s>>(L-1))] for s in p])
            rec['levels']=[]
            for k in range(3):
                vk=vs[:,order[k]];pk=np.abs(vk)**2;nk=float(pk@N);tk=np.vdot(vk[translation],vk)
                rec['levels'].append(dict(E=float(e[k]),N=nk,varN=float(pk@(N*N)-nk*nk),
                                          translation_phase=float(np.angle(tk)),translation_abs=float(abs(tk))))
            out['rows'].append(rec);print(rec,flush=True)
    out['elapsed_seconds']=time.time()-start
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
