"""Independent free-fermion checks and infinite-MPS Ising control summaries."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.sparse.linalg import eigsh
from independent_ed import ed
HERE=Path(__file__).resolve().parent
rows=[]
for model,h in [('ising_critical',1.),('ising_gapped',1.01)]:
    exact=-quad(lambda k:np.sqrt(1+h*h-2*h*np.cos(k)),0,np.pi,epsabs=1e-12,epsrel=1e-12)[0]/np.pi
    L=8;p=dict(Jx=-1.,Jy=0.,Jz=0.,Dz=0.,hz=-h)
    H,basis=ed.build_sector(L,p,parity=0)
    E=float(eigsh(H,k=1,which='SA',tol=1e-13,return_eigenvectors=False)[0])
    k=(2*np.arange(L)+1)*np.pi/L
    analytic=-float(np.sum(np.sqrt(1+h*h-2*h*np.cos(k))))
    assert abs(E-analytic)<1e-10
    points=[]
    for f in sorted((HERE/'results').glob(f'{model}_*.json')):
        r=json.loads(f.read_text())
        if r['status']!='done':continue
        err=r['E_physical']-exact
        assert err>-1e-8,(f.name,err)
        points.append(dict(source=str(f.relative_to(HERE)),chi=r['chi'],converged=r['converged'],
            one_site=r['config'].get('one_site',False),
            energy_error=err,S=float(np.mean(r['entropy']['1'])),xi=r['correlation_lengths'][0],
            order=float(np.sqrt(np.mean(np.array(r['X'])**2+np.array(r['Y'])**2)))))
    rows.append(dict(model=model,h=h,exact_energy_density=exact,exact_fermion_gap=2*abs(h-1),
        L8_even_energy_error=abs(E-analytic),points=sorted(points,key=lambda p:p['chi'])))
(HERE/'ising_controls.json').write_text(json.dumps(dict(controls=rows),indent=2)+'\n')
for r in rows:print(r)
