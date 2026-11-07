"""Resumable C135 energy scan. Large mode requires a SLURM allocation."""
from core import *
import argparse
def models():
    out=[]
    for d in [1,.5,2,-.9,-.99,-.9999,0]:
        name='Heisenberg' if d==1 else ('XX' if d==0 else f'XXZ_{d}')
        out.append((name,xxz(d),'u1'))
    out.extend([('XX_field',xxz(0,.5),'u1'),('XXZ_DM',xxz(.5,0,.7),'u1'),('chiral',xxz(1,0,0)*0+B[2]+B[5]-.25*np.eye(4),'u1'),('TFIM_1',tfim(1),'parity')])
    for n in range(2,9):
        for c in [-2,.5,5]:out.append((f'TFIM_c{c}_scale{n}',tfim(1+c/n),'parity'))
    rng=np.random.default_rng(1352026);F=np.eye(4)[[0,2,1,3]]
    for j in range(30):
        a=rng.normal(size=(4,4))
        if j>=15:a=a+1j*rng.normal(size=(4,4))
        h=herm(a)
        if j<8:h=(h+F@h@F)/2
        h=h-np.trace(h)*np.eye(4)/4
        m=gauge(h);h=h/m['upper']
        out.append((f'random_{j:02}',h,'none'))
    return out
def run(args):
    if args.large and not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Large computations require SLURM')
    Path('results').mkdir(exist_ok=True);Path('matrices').mkdir(exist_ok=True)
    mods=models();write_json('models.json',[dict(name=name,h=pack(h),sym=sym,M=gauge(h)) for name,h,sym in mods])
    if args.index is not None:mods=[mods[args.index]]
    elif args.select:mods=[v for v in mods if v[0] in args.select.split(',')]
    for name,h,sym in mods:
        path=Path(f'results/{name}.json');out=json.loads(path.read_text()) if path.exists() else dict(name=name,M=gauge(h),h=pack(h),sym=sym,lti={},rings={})
        nmax=10 if args.large else (8 if sym!='none' else (6 if np.max(abs(h.imag))<1e-13 else 4))
        if args.nmax:nmax=args.nmax
        Lmax=20 if args.large else (12 if name in ['Heisenberg','XX','XX_field','TFIM_1'] else 10)
        if args.Lmax:Lmax=args.Lmax
        for n in range(2,nmax+1):
            if str(n) in out['lti']:continue
            real=np.max(abs(h.imag))<1e-13
            rec,r,x=LTI(n,real,sym).solve(h,eps=2e-8,iters=12000)
            # Retry inaccurate small Clarabel solves with an independent solver.
            if rec['status']!='optimal' and n<=5:
                old=rec;rec,r,x=LTI(n,real,sym).solve(h,solver='SCS',eps=2e-9,iters=25000);rec['first_attempt']=old
            out['lti'][str(n)]=rec
            np.savez(f'matrices/{name}_{n}.npz',rho=r,X=x,h=h)
            write_json(path,out);print(name,'n',n,rec['value'],rec['status'],rec['seconds'],flush=True)
        for L in range(3,Lmax+1):
            if str(L) in out['rings']:continue
            rec,v=ring(h,L)
            rec['rho2']=pack(ti_marginal(v,L,2)); rec['petz']=[]
            windows=sorted(set([max(2,L//2),min(L-1,max(2,L//2+1))]))
            for ell in windows:
                rho=ti_marginal(v,L,ell);tau,diag=petz(rho)
                rec['petz'].append(dict(ell=ell,tau=pack(tau),energy=float(np.trace(h@tau).real),diagnostics=diag))
            out['rings'][str(L)]=rec;write_json(path,out)
            print(name,'L',L,rec['e'],rec['residual'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--large',action='store_true');ap.add_argument('--index',type=int);ap.add_argument('--select');ap.add_argument('--nmax',type=int);ap.add_argument('--Lmax',type=int)
    run(ap.parse_args())
