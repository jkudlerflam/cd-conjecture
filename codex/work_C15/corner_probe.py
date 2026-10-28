"""Finite-size probe near the singular two-minimum corner. No critical tuning assumed."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json,numpy as np
from checks import full_case
out=[]
for L in (10,12,14,16):
    r=full_case('singular_corner_probe',L,dict(Jx=.49,Jy=-.51,Jz=2.,D=.01,g=float(np.sqrt(15.))))
    out.append(r); print(L,r['c_chord'],r['S_half'],r['levels'],flush=True)
json.dump(out,open('corner_probe.json','w'),indent=2)
