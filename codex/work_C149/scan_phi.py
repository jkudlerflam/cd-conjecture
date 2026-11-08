"""C149: size of the quantile-trial residual Phi_j = -F_j(y0)/(2 pi L) for y0 = quantiles of r_Q (m(Q)=nu),
and of the true displacement L(Z_0(x_j) - I_j/L). Floating point only. Usage: python3 scan_phi.py eps L1 L2 ..."""
import numpy as np, sys, json
from multiprocessing import Pool
from bethe_ferro import Model, residual, solve
from continuum import Kernel
from trial import continuum_trial, quantiles

def case(args):
  try:
    return case0(args)
  except Exception as e:
    return dict(eps=args[0],L=args[1],n=args[2],error=str(e))

def case0(args):
    eps, L, n = args
    eta=2*np.arcsin(np.sqrt(eps/2)); m=Model(eta)
    R,P=continuum_trial(eta,L,n,N=700)
    y=quantiles(P.inv,L,n)
    Phi=-residual(m,y,L)/(2*np.pi*L)
    d=R.Q-np.abs(y); lam=eta*L*P.rho_at(y)
    G,x,_=solve(L,n,eta)
    I=np.arange(1,n+1)-(n+1)/2
    disp=L*(P.Z(x)-I/L)
    e=d<5*eta
    out=dict(eps=eps,L=L,n=n,Q=R.Q,lam0=eta*L*R.r[0],lamQ=eta*L*R.r[-1],ZQ=R.Z[-1],
        maxPhi=float(L*np.abs(Phi).max()/eta), maxPhi_edge=float(L*np.abs(Phi[e]).max()/eta) if e.any() else 0.0,
        maxPhi_bulk=float(L*np.abs(Phi[~e]).max()/eta) if (~e).any() else 0.0,
        sgn_right_edge=float(np.sign(Phi[-1])), maxdisp=float(np.abs(disp).max()),
        maxdisp_lam2=float(np.abs(disp[lam>=2]).max()) if (lam>=2).any() else 0.0)
    return out

if __name__=='__main__':
    eps=float(sys.argv[1]); Ls=[int(a) for a in sys.argv[2:]]
    jobs=[(eps,L,n) for L in Ls for n in sorted(set([L//16,L//8,L//4,3*L//8,L//2-8,L//2-2]))]
    with Pool(8) as p:
        res=p.map(case,jobs)
    for r in res: print(json.dumps({k:(round(v,4) if isinstance(v,float) else v) for k,v in r.items()}))
    json.dump(res,open(f'scan_phi_{eps:g}.json','w'),indent=1)
