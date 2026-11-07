"""Mesh-convergence check of the finite-boundary diagnostics (dN and L(F-L(E+h(nu-m))))."""
import sys, json
import numpy as np
from bethe_ferro import solve, Model
from continuum import Kernel, Restricted, grading
from scipy.optimize import brentq
out=[]
for (k,L,n) in [(6,2048,512),(6,2048,896),(5,2048,768),(4,2048,768)]:
    eps=10.0**(-k); eta=2*np.arcsin(np.sqrt(eps/2)); mod=Model(eta); ker=Kernel(eta)
    G,x,d=solve(L,n,eta); nu=n/L; F=mod.Delta*L/4+G
    f=lambda Q: mod.p(Q)-np.sum(mod.phi(Q-x))/L
    QL=brentq(f,x[-1]+1e-12,x[-1]+50,xtol=1e-15)
    res=[]
    for N in [700,1400,2800]:
        c=grading(ker,QL,700)
        R=Restricted(ker,QL,N=N,c=c)
        res.append((N,R.m,R.E,R.h))
    (N1,m1,E1,h1),(N2,m2,E2,h2),(N3,m3,E3,h3)=res
    mR=(4*m3-m2)/3; ER=(4*E3-E2)/3; hR=(4*h3-h2)/3
    mR2=(4*m2-m1)/3; ER2=(4*E2-E1)/3
    row=dict(k=k,L=L,n=n,QL=QL,dN_rich_2800=L*(nu-mR),dN_rich_1400=L*(nu-mR2),
             TII_over_eta_2800=L*(F-L*(ER+hR*(nu-mR)))/eta,
             TII_over_eta_1400=L*(F-L*(ER2+h2*(nu-mR2)))/eta,
             raw_m=[r[1] for r in res])
    print(json.dumps(row),flush=True); out.append(row)
json.dump(out,open('conv_check.json','w'),indent=1)
