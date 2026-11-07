"""C135 distance brackets. All reported bounds are floating-point checks."""
from core import *
import argparse
def lower_distance(target,n,real=False,sym='none',eps=2e-8):
    start=time.time();m=LTI(n,real,sym)
    tau=m.r
    for j in range(n-2):tau=tau[::2,::2]+tau[1::2,1::2]
    target=target.real if real else target
    D=target-tau;S0=cp.Variable((4,4),symmetric=True) if real else cp.Variable((4,4),hermitian=True)
    a=S0-D>>0;b=S0+D>>0
    obj=cp.trace(S0) if real else cp.real(cp.trace(S0))
    problem=cp.Problem(cp.Minimize(obj),m.constraints+[a,b])
    problem.solve(solver='SCS',eps=eps,max_iters=16000)
    r=herm(m.r.value);x=herm(m.ec.dual_value);Y=herm(a.dual_value-b.dual_value)
    # The final eigenvalue bound is valid for any Hermitian Y and X after
    # norm repair, independent of the solver's dual scaling convention.
    scale=max(1.,la.norm(Y,2));Y/=scale;x/=scale
    K=-np.kron(Y,np.eye(2**(n-2)))+np.kron(x,I)-np.kron(I,x)
    t=float(la.eigvalsh(herm(K),subset_by_index=[0,0])[0])-2e-11
    low=float(np.trace(Y@target).real)+t
    q,mix,eq=repair(r,n)
    dist=float(np.abs(la.eigvalsh(target-marginal(q,n))).sum())
    return dict(n=n,value=float(problem.value),lower=max(0,low),primal_upper_to_LTI=dist,status=problem.status,seconds=time.time()-start,trace_residual=float(abs(np.trace(r)-1)),lti_residual=float(la.norm(left(r)-right(r))),min_rho=float(la.eigvalsh(r,subset_by_index=[0,0])[0]),dual_norm=float(la.norm(Y,2)),repair_mix=mix,scs={k:problem.solver_stats.extra_stats['info'][k] for k in ['res_pri','res_dual','gap','status']})
def hull_distance(target,dictionary):
    A=np.stack([t.reshape(-1) for t in dictionary],axis=1);w=cp.Variable(len(dictionary),nonneg=True)
    tau=cp.reshape(A@w,(4,4),order='C');D=target-tau
    S0=cp.Variable((4,4),hermitian=True)
    p=cp.Problem(cp.Minimize(cp.real(cp.trace(S0))),[cp.sum(w)==1,S0-D>>0,S0+D>>0])
    p.solve(solver='CLARABEL',tol_feas=1e-9,tol_gap_abs=1e-9,tol_gap_rel=1e-9)
    ww=np.maximum(w.value,0);ww/=ww.sum();t=(A@ww).reshape(4,4)
    return dict(value=float(p.value),upper=float(abs(la.eigvalsh(herm(target-t))).sum()),status=p.status,active=int((ww>1e-7).sum()),weights=[[int(i),float(ww[i])] for i in np.flatnonzero(ww>1e-7)])
def run(args):
    if args.large and not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Large run requires SLURM')
    files=sorted(Path('results').glob('*.json'));records=[json.loads(p.read_text()) for p in files]
    dictionary=[];labels=[]
    for rec in records:
        for L,r in rec['rings'].items():
            for p in r['petz']:
                # Reject badly normalized cutoff results, never silently renormalize.
                if p['diagnostics']['trace_residual']<1e-8:
                    dictionary.append(unpack(p['tau']));labels.append([rec['name'],int(L),p['ell']])
    # D=1 uniform MPS examples. These are exact points of Omega_inf.
    for s in [P[1],P[2],P[3]]:
        for sign in [-1,1]:
            r=(I+sign*s)/2;dictionary.append(np.kron(r,r));labels.append(['product',int(len(labels)),1])
    write_json('dictionary.json',dict(labels=labels,points=[pack(t) for t in dictionary]))
    if args.index is not None:records=[records[args.index]]
    if args.select:records=[r for r in records if r['name'] in args.select.split(',')]
    for rec in records:
        name=rec['name'];path=Path(f'results_distance/{name}.json');Path('results_distance').mkdir(exist_ok=True)
        out=json.loads(path.read_text()) if path.exists() else dict(name=name,dictionary_size=len(dictionary),rows={})
        for L in range(3,9):
            if str(L) not in rec['rings']:continue
            target=unpack(rec['rings'][str(L)]['rho2']);row=out['rows'].get(str(L),{})
            # Hull upper depends on all current dictionary points, so recompute on resume.
            row['hull']=hull_distance(target,dictionary);row['L2_upper']=L*L*row['hull']['upper']
            n=L+1;nmax=9 if args.large else (6 if abs(target.imag).max()<1e-12 else 4)
            if args.max_n:nmax=min(nmax,args.max_n)
            if n<=nmax and 'lti' not in row:
                real=abs(target.imag).max()<1e-12
                # Symmetry reduction is valid for distance only if the target shares it.
                sym=rec['sym'];q=np.diag([sum((i>>j)&1 for j in range(2)) for i in range(4)])
                if sym=='u1' and la.norm(q@target-target@q)>1e-9:sym='none'
                if sym=='parity' and la.norm(np.kron(P[3],P[3])@target-target@np.kron(P[3],P[3]))>1e-9:sym='none'
                row['lti']=lower_distance(target,n,real,sym);row['L2_lower']=L*L*row['lti']['lower']
            out['rows'][str(L)]=row;write_json(path,out)
            print(name,L,row.get('L2_lower'),row['L2_upper'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--large',action='store_true');ap.add_argument('--index',type=int);ap.add_argument('--select');ap.add_argument('--max-n',type=int);run(ap.parse_args())
