import json,numpy as np
from scipy.signal import find_peaks
from scipy.optimize import minimize_scalar
from pair_odd_screen import energy
out=[]
ks=np.linspace(0,np.pi,1001)
for delta in [0,.15,.3,.45,.6,.75,.9,1.05,1.2]:
 es=np.array([energy(k,-14/9+delta,0,2,10/9,n=100) for k in ks])
 ix=find_peaks(-es)[0]
 minima=[]
 for j in ix:
  if j<2 or j>len(ks)-3:continue
  rr=minimize_scalar(lambda k:energy(k,-14/9+delta,0,2,10/9,n=160),bounds=(ks[j-1],ks[j+1]),method='bounded')
  minima.append(dict(Q=float(rr.x),E=float(rr.fun),binding_to_same_K_continuum=float(-4*np.cos(rr.x/2)-rr.fun)))
 row=dict(delta=delta,minima=minima)
 out.append(row);print(json.dumps(row),flush=True)
with open('detuned_pair_bands_results.json','w') as f:json.dump(out,f,indent=2)
