"""Low levels within selected momentum blocks; no conformal labels assumed."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import json,time
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import eigsh
from spectrum import orbits,build
HERE=Path(__file__).resolve().parent
p=dict(Jx=.5,Jy=-.5,Jz=.625,Dz=.1,hz=-.875)
out=HERE/'low_levels.json';rows=json.loads(out.read_text())['rows'] if out.exists() else []
for L in [16,20,24]:
    for P in [0,1]:
        data=orbits(L,P)
        for k in ([0,L//2] if P==0 else [0]):
            if any(r['L']==L and r['parity']==P and r['k']==k for r in rows):continue
            start=time.time();H,defect=build(L,k,p,data)
            val,vec=eigsh(H,k=4,which='SA',tol=1e-11,maxiter=100000,
                v0=np.random.default_rng(391+L+P+k).normal(size=H.shape[0]))
            order=np.argsort(val);val=val[order];vec=vec[:,order]
            residuals=[float(np.linalg.norm(H@vec[:,i]-val[i]*vec[:,i])) for i in range(4)]
            assert max(residuals)<1e-8
            row=dict(L=L,parity=P,k=k,energies=val.tolist(),residuals=residuals,
                dimension=H.shape[0],hermiticity_defect=defect,seconds=time.time()-start)
            rows.append(row);tmp=out.with_suffix('.tmp');tmp.write_text(json.dumps(dict(couplings=p,rows=rows,complete=len(rows)==9),indent=2)+'\n');tmp.replace(out)
            print(json.dumps(row),flush=True)
