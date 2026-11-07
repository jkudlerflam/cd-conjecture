"""(D*) at finite Fermi boundary: eta L sup|sigma_L - r_{Q_L}| on [-Q_L,Q_L] (fine grid),
and the profile of the cumulative error G (centered) to see whether its growth is smooth."""
import json
import numpy as np
from bethe_ferro import solve, Model
from continuum import Kernel, Restricted, grading
from scipy.optimize import brentq
out=[]
for (k,L,n) in [(3,512,192),(4,1024,384),(4,2048,768),(5,1024,384),(5,2048,768),(6,2048,512),(6,2048,128)]:
    eps=10.0**(-k); eta=2*np.arcsin(np.sqrt(eps/2)); mod=Model(eta); ker=Kernel(eta)
    G,x,d=solve(L,n,eta); nu=n/L
    f=lambda Q: mod.p(Q)-np.sum(mod.phi(Q-x))/L
    QL=brentq(f,x[-1]+1e-12,x[-1]+50,xtol=1e-15)
    R=Restricted(ker,QL,N=1400,c=grading(ker,QL,700))
    X=np.concatenate([np.linspace(0,QL,20000), x[x>=0]])
    sig=ker.a1(X)+ker.T(X[:,None]-x[None,:]).sum(axis=1)/L
    r=R.interp(R.r,X)
    D=eta*L*np.abs(sig-r)
    # where is the max: distance from edge, in units of eta
    i=np.argmax(D)
    # smooth part of G: average of centered G over each cell pair, at roots j: G(x_j^-)+G(x_j^+) /2
    yy=R.y; Cr=np.concatenate([[0],np.cumsum((R.r[1:]+R.r[:-1])/2*np.diff(yy))])
    Ix=np.sign(x)*np.interp(np.abs(x),yy,Cr)
    j=np.arange(1,n+1); Gmid=L*((j-0.5)/L-nu/2-Ix)
    row=dict(k=k,L=L,n=n,etaL=eta*L,QL=QL,Dmax=float(D.max()),Dmax_at_dist_over_eta=float((QL-X[i])/eta),
             D_interior=float(D[X<QL-0.2].max()), LGmid_max=float(np.abs(Gmid).max()),
             LGmid_at_center=float(Gmid[n//2]), LGmid_at_edge=float(Gmid[-1]))
    print(json.dumps(row),flush=True); out.append(row)
json.dump(out,open('dens_check.json','w'),indent=1)
