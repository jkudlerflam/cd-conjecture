"""Full matrix intertwining check for every momentum on an eight-site ring."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import sys,json
from pathlib import Path
import numpy as np
import spectrum as sp
HERE=Path(__file__).resolve().parent
from independent_ed import ed
p=dict(Jx=.5,Jy=-.5,Jz=.625,Dz=.1,hz=-.875);L=8
rows=[]
for P in [0,1]:
    H,basis=ed.build_sector(L,p,parity=P);inv=np.full(1<<L,-1,int);inv[basis]=np.arange(len(basis))
    reps,ids,shifts,periods=data=sp.orbits(L,P)
    for k in range(L):
        h,defect=sp.build(L,k,p,data);keep=np.flatnonzero((periods.astype(int)*k)%L==0)
        U=np.zeros((len(basis),len(keep)),complex)
        for col,i in enumerate(keep):
            for s in basis[ids[basis]==i]:
                U[inv[s],col]=np.exp(-2j*np.pi*k/L*int(shifts[s]))/np.sqrt(float(periods[i]))
        err=np.linalg.norm(H@U-U@h);norm=np.linalg.norm(U.conj().T@U-np.eye(len(keep)))
        assert err<1e-10 and norm<1e-12,(P,k,err,norm)
        rows.append(dict(parity=P,k=k,intertwiner_error=err,isometry_error=norm))
(HERE/'spectrum_validation.json').write_text(json.dumps(rows,indent=2)+'\n')
print('Largest intertwiner error',max(r['intertwiner_error'] for r in rows))
