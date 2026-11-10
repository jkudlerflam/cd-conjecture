"""Complete momentum/parity scan of the lowest ring energy in each sector."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import eigsh
from spectrum import orbits,build
HERE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--sizes',type=int,nargs='+',default=[16,20,24]);a=ap.parse_args()
p=dict(Jx=.5,Jy=-.5,Jz=.625,Dz=.1,hz=-.875)
for L in a.sizes:
    out=HERE/f'dispersion_L{L}.json';rows=json.loads(out.read_text())['sectors'] if out.exists() else []
    for P in [0,1]:
        data=orbits(L,P)
        for k in range(L):
            if any(r['parity']==P and r['k']==k for r in rows):continue
            start=time.time();H,defect=build(L,k,p,data)
            val,v=eigsh(H,k=1,which='SA',tol=1e-11,maxiter=100000,v0=np.random.default_rng(100+P*L+k).normal(size=H.shape[0]))
            residual=float(np.linalg.norm(H@v[:,0]-val[0]*v[:,0]));assert residual<1e-8
            row=dict(parity=P,k=k,E=float(val[0]),residual=residual,dimension=H.shape[0],seconds=time.time()-start,hermiticity_defect=defect)
            rows.append(row);result=dict(L=L,couplings=p,sectors=rows,complete=len(rows)==2*L)
            if result['complete']:
                ground=min(rows,key=lambda r:r['E']);E=ground['E'];result.update(ground=ground,ground_sectors=[r for r in rows if r['E']<E+1e-8],
                    parity_gap=min(r['E'] for r in rows if r['parity']==1)-min(r['E'] for r in rows if r['parity']==0))
            tmp=out.with_suffix('.tmp');tmp.write_text(json.dumps(result,indent=2)+'\n');tmp.replace(out)
            print(L,P,k,row['E'],flush=True)
