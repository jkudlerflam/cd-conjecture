"""Location of the outermost root: L r0(x_n)/b against the proved b(1-q_J)/4 and the needed (2-ln2)/pi^2.
Also min over gaps of sigma_L/r0 against theta(gamma)=1-8(2-ln2)(pi-gamma)/pi^3. Floating point. J=1."""
import json
import numpy as np
from bethe_scan import solve
out=[]
for delta in [0.1,0.5,0.9,0.98,0.99,0.999,0.9999]:
    g=np.arccos(delta); b=np.pi/g; qJ=(np.pi-2*g)/(2*(np.pi-g))
    a1=lambda x: np.sin(g)/(np.pi*(np.cosh(2*x)-np.cos(g)))
    K=lambda x: np.sin(2*g)/(2*np.pi*(np.sinh(x)**2+np.sin(g)**2))
    r0=lambda x: 1/(2*g*np.cosh(b*x))
    theta=1-8*(2-np.log(2))*(np.pi-g)/np.pi**3
    worst=(np.inf,None); wg=(np.inf,None)
    for L in [2,3,4,5,8,9,16,17,64,65,256,257,1024]:
        for n in sorted(set([1,2,3,L//8,L//4,L//3,L//2-1,L//2])):
            if n<2 or n>L//2: continue
            E,d,x=solve(L,n,delta,return_roots=True)
            v=L*r0(x[-1])/b
            if v<worst[0]: worst=(v,(L,n))
            ys=np.linspace(0,x[-1],4001)
            s=a1(ys)-K(ys[:,None]-x[None,:]).sum(1)/L
            m=float((s/r0(ys)).min())
            if m<wg[0]: wg=(m,(L,n))
    row=dict(delta=delta,min_Lr0xn_over_b=worst[0],at=worst[1],proved_lower=(1-qJ)/4,
             needed=(2-np.log(2))/np.pi**2,min_gap_sig_over_r0=wg[0],at_gap=wg[1],theta_proved=theta)
    out.append(row); print(row,flush=True)
json.dump(out,open('location_scan.json','w'),indent=0)
