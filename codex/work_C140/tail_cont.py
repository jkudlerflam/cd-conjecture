"""Test of the sufficient continuum tail inequality (T_c):
   s_Q(x) >= 2 (int_Q^inf s_Q) phi_Q(x) for x>Q,  phi_Q = -d s_Q/dQ / s_Q(Q).
Floating point. J=1."""
import json, sys
import numpy as np
from continuum_rep import Model, Restricted
out=[]
for delta in [0.1,0.5,0.9,0.99]:
    M=Model(delta); g=M.g
    for Q in [0.05,0.2,0.5,1.0,2.0,3.0]:
        if Q>12*g: continue
        dQ=1e-4*max(g,Q)
        xs=Q+np.concatenate([np.linspace(1e-6,3*g,200),np.linspace(3*g,12,400)])
        S=Restricted(M,Q,N=10); Sp=Restricted(M,Q+dQ,N=10); Sm=Restricted(M,Q-dQ,N=10)
        sQ=S.s_full(xs)
        ds=(Sp.s_full(xs)-Sm.s_full(xs))/(2*dQ)
        sQQ=S.s_full(np.array([Q]))[0]
        phi=-ds/sQQ
        tail=(0.5-S.m)/(1-M.qJ)/2
        ratio=sQ/(2*tail*phi)
        row=dict(delta=delta,Q=Q,min_ratio=float(ratio.min()),at=float(xs[np.argmin(ratio)]-Q),
                 ratio_at_Q=float(ratio[0]),ratio_far=float(ratio[-1]))
        out.append(row); print(row,flush=True)
json.dump(out,open('tail_cont.json','w'),indent=0)
