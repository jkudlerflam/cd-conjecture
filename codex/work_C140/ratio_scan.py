"""min over R of sigma_L/r0 and Lsigma_L at roots, 0<Delta<=1. Floating point diagnostics. J=1."""
import json, sys
import numpy as np
from bethe_scan import solve
def funcs(delta):
    if delta<1:
        g=np.arccos(delta); b=np.pi/g
        a1=lambda x: np.sin(g)/(np.pi*(np.cosh(2*x)-np.cos(g)))
        K=lambda x: np.sin(2*g)/(2*np.pi*(np.sinh(x)**2+np.sin(g)**2))
        r0=lambda x: 1/(2*g*np.cosh(b*x))
        J0=(np.pi-2*g)/((np.pi-g)*2*g)
    else:
        b=np.pi
        a1=lambda x: 2/(np.pi*(1+4*x*x)); K=lambda x: 1/(np.pi*(1+x*x))
        r0=lambda x: 1/(2*np.cosh(np.pi*x)); J0=1/2
    return a1,K,r0,b,J0
out=[]
deltas=[float(v) for v in sys.argv[1].split(',')]
Ls=[int(v) for v in sys.argv[2].split(',')]
for delta in deltas:
    a1,K,r0,b,J0=funcs(delta)
    for L in Ls:
        ns=sorted(set([1,2,3,L//8,L//4,L//3,L//2-2,L//2-1,L//2]))
        for n in ns:
            if n<1 or n>L//2: continue
            E,d,x=solve(L,n,delta,return_roots=True)
            X=x.max()
            ys=np.concatenate([np.linspace(0,X+12/b,6000),x[x>=0]])
            s=np.empty_like(ys)
            for i in range(0,len(ys),500):
                yy=ys[i:i+500]; s[i:i+500]=a1(yy)-K(yy[:,None]-x[None,:]).sum(1)/L
            rat=s/r0(ys)
            row=dict(delta=delta,L=L,n=n,min_ratio=float(rat.min()),
                     at_over_xn=float(ys[np.argmin(rat)]/X) if X>0 else 0.,
                     min_Lsig_root_over_b=float((L*s[len(ys)-len(x[x>=0]):]).min()/b),
                     min_sig=float(s.min()))
            out.append(row); print(row,flush=True)
json.dump(out,open(f'ratio_scan_{sys.argv[3]}.json','w'),indent=0)
