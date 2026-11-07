"""Eight-angle multistart search. Limited local searches are exploratory.
Large mode includes n=6 to n=8 LTI differences and larger dictionaries.
"""
from core import *
from scipy.optimize import minimize
import argparse
def from_angles(a):
    c=np.empty(9);s=1.
    for i in range(8):c[i]=s*np.cos(a[i]);s*=np.sin(a[i])
    c[8]=s;return c
def to_angles(c):
    c=c/la.norm(c);return np.r_[[np.arctan2(la.norm(c[i+1:]),c[i]) for i in range(7)],np.arctan2(c[8],c[7])]
def evaluate(c,kind,L,n,large=False):
    h=np.einsum('a,aij->ij',c,B);m=gauge(h)['upper'];h/=m
    if kind=='lti':
        rows=[]
        for j in [n,n+2]:rows.append(LTI(j,abs(h.imag).max()<1e-13,'none').solve(h,solver='SCS',eps=2e-7,iters=1500 if not large else 15000)[0])
        return n*n*(rows[1]['lower']-rows[0]['upper']),dict(lti=rows,coefficients=(c/m).tolist(),n=n)
    rec,v=ring(h,L);r=ti_marginal(v,L,max(2,L//2));tau,diag=petz(r)
    # Upper endpoint of the directional bracket, a true variational objective
    # in exact arithmetic, not a lower witness against (K).
    up=float(np.trace(h@tau).real);score=L*L*(up-rec['e'])
    return score,dict(L=L,ring=rec,e0_upper=up,petz_diagnostics=diag,coefficients=(c/m).tolist())
def run(args):
    if args.large and not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Large searches require SLURM')
    Path('search').mkdir(exist_ok=True);kinds=[('lti',0,4 if args.large else 2),('ring_upper',6,0),('ring_upper',8,0),('ring_upper',10,0)]
    if args.large:kinds[0]=('lti',0,6)
    if args.kind:kinds=[v for v in kinds if v[0]==args.kind and (not args.L or v[1]==args.L)]
    for kind,L,n in kinds:
        for seed in range(args.starts):
            target=Path(f'search/{kind}_L{L}_n{n}_seed{seed}.json')
            if target.exists() and json.loads(target.read_text())['evaluations']>=args.budget:continue
            rng=np.random.default_rng(135400+seed)
            c=np.r_[np.ones(3),np.zeros(6)] if seed==0 else rng.normal(size=9);c/=la.norm(c)
            base=np.ones(9)/3;v=base-c;rot=np.eye(9)-2*np.outer(v,v)/(v@v) if la.norm(v)>1e-12 else np.eye(9)
            a=to_angles(base);best=[-np.inf,None];count=0;start=time.time()
            def loss(ang):
                nonlocal count
                cc=rot@from_angles(ang);q,rec=evaluate(cc,kind,L,n,args.large);count+=1
                if q>best[0]:best[:]=[q,rec]
                return -q
            result=minimize(loss,a,method='Nelder-Mead',options=dict(maxfev=args.budget,xatol=1e-5,fatol=1e-7,initial_simplex=np.vstack([a,a+.12*np.eye(8)])))
            out=dict(kind=kind,L=L,n=n,seed=seed,score=best[0],point=best[1],evaluations=count,success=bool(result.success),message=str(result.message),seconds=time.time()-start,tag='numerical exploratory; finite-budget local search')
            write_json(target,out);print(json.dumps(out),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--large',action='store_true');ap.add_argument('--starts',type=int,default=2);ap.add_argument('--budget',type=int,default=24);ap.add_argument('--kind');ap.add_argument('--L',type=int);run(ap.parse_args())
