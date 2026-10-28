"""Exact pairing-only factorization-line check, NumPy/SciPy only."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
import json
from checks import ham
V=4.;delta=1.;mu=V-np.sqrt(V*V-4*delta*delta);h=(V-mu)/2
out=[]
for L in (4,6,8):
    H,b=ham(L,delta/2,-delta/2,V/4,0,h)
    e=np.linalg.eigvalsh(H.toarray());n=int(np.sum(abs(e-e[0])<1e-9))
    k=4*np.pi/L
    branch=V-np.sqrt(V*V-4*delta*delta*np.sin(k/2)**2)
    r=dict(L=L,V=V,delta=delta,mu=float(mu),h=float(h),E0=float(e[0]),ground_degeneracy=n,
           ground_energy_error=float(e[0]+V*L/4),first_gap=float(e[n]-e[0]),
           predicted_defect_branch=float(branch),branch_spectrum_distance=float(np.min(abs(e-e[0]-branch))))
    out.append(r);print(r)
open('factorization.json','w').write(json.dumps(out,indent=2))
