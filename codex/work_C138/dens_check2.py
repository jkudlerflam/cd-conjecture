"""(D*) and (D1*) over a wider anisotropy range, finite and infinite Fermi boundary.
 Dmax  = eta L sup_x |sigma_L - r_Q|     on [-Q_L,Q_L] (fine grid plus roots)
 D1max = eta^2 L sup_x |sigma_L' - r_Q'| on [-Q_L+eta/4, Q_L-eta/4] (r_Q' from the mesh)
Floating point only."""
import json, sys
import numpy as np
from bethe_ferro import solve, Model
from continuum import Kernel, Restricted, grading
from scipy.optimize import brentq
out=[]
import sys
EPS=[float(e) for e in sys.argv[1:]] or [0.9, 0.5, 0.1, 1e-2, 1e-4, 1e-6]
for eps in EPS:
    eta=2*np.arcsin(np.sqrt(eps/2)); mod=Model(eta); ker=Kernel(eta)
    for L in [16, 64, 256, 1024]:
        for n in sorted(set([1, 2, int(np.sqrt(L)), L//8, L//4, 3*L//8, L//2-2])):
            if n < 2: continue
            G,x,d=solve(L,n,eta); nu=n/L
            f=lambda Q: mod.p(Q)-np.sum(mod.phi(Q-x))/L
            QL=brentq(f,x[-1]+1e-12,x[-1]+60,xtol=1e-15)
            R=Restricted(ker,QL,N=1000,c=grading(ker,QL,1000))
            X=np.concatenate([np.linspace(0,QL,6000), x[x>=0]])
            sig=ker.a1(X)+ker.T(X[:,None]-x[None,:]).sum(axis=1)/L
            r=R.interp(R.r,X)
            hh=1e-6*max(eta,1e-3)
            Tp=(ker.T(X[:,None]-x[None,:]+hh)-ker.T(X[:,None]-x[None,:]-hh))/(2*hh)
            dsig=(ker.a1(X+hh)-ker.a1(X-hh))/(2*hh)+Tp.sum(axis=1)/L
            gr=np.gradient(R.r,R.y); gr[0]=0.0   # r_Q is even: r_Q'(0)=0 (one-sided difference at x=0 is a mesh artifact)
            dr=np.interp(X,R.y,gr)
            inner=X<QL-min(eta/4,QL/4)
            row=dict(eps=eps,eta=eta,L=L,n=n,etaL=eta*L,QL=QL,
                     Dmax=float(eta*L*np.max(np.abs(sig-r))),
                     D1max=float(eta**2*L*np.max(np.abs(dsig-dr)[inner])))
            print(json.dumps(row),flush=True); out.append(row)
json.dump(out,open('dens_check2%s.json'%('' if not sys.argv[1:] else '_'+'_'.join(sys.argv[1:])),'w'),indent=1)
print('MAX Dmax',max(r['Dmax'] for r in out),'MAX D1max',max(r['D1max'] for r in out))
