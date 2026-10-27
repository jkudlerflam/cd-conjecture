import json
import numpy as np
from scipy.optimize import minimize_scalar
from pair_odd_screen import energy

g=2.0
s=0.5
tp=g*g*(2+s)/(4*(1+s)**2)
U=-4+2*tp-4*tp*s*s+2*g*g*s*s/(1+s)
q=2*np.arcsin(s)
res=minimize_scalar(lambda K:energy(K,U,0,g,tp,n=160),bounds=(0.001,np.pi-0.001),method='bounded',options={'xatol':1e-13})
delta=1e-4
curv=(energy(q+delta,U,0,g,tp,160)+energy(q-delta,U,0,g,tp,160)-2*energy(q,U,0,g,tp,160))/delta**2
out=dict(g=g,tp=tp,U=U,V=0,Q=q,K_numeric=res.x,E_numeric=res.fun,E_at_0=energy(0,U,0,g,tp,160),E_at_pi=energy(np.pi,U,0,g,tp,160),pair_curvature=curv,pair_mass=1/curv,atomic_mass=0.5)
print(json.dumps(out,indent=2))
with open('pair_threshold_results.json','w') as f:json.dump(out,f,indent=2)
