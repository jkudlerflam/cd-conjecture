"""Independent matrix and spectrum checks of exact antiunitary pairing."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import sys,json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
from independent_ed import ed
p=dict(Jx=.5,Jy=-.5,Jz=.625,Dz=.1,hz=-.875);checks=[]
L=8
for P in [0,1]:
    H,basis=ed.build_sector(L,p,parity=P);H=H.toarray();N=len(basis)
    z=np.array([1-2*((basis>>(L-1-j))&1) for j in range(L)])
    v=np.exp(1j*np.pi/4*np.sum((-1.)**np.arange(L)[:,None]*z,axis=0))
    inv={int(s):i for i,s in enumerate(basis)};T=np.zeros((N,N))
    for i,s in enumerate(basis):T[inv[(int(s)>>1)|((int(s)&1)<<(L-1))],i]=1
    V=np.diag(v);sgn=(-1)**P
    errors=dict(hamiltonian=float(np.max(np.abs(V@H.conj()-H@V))),
        unitary_square=float(np.max(np.abs(v*v-sgn))),
        antiunitary_square=float(np.max(np.abs(v*v.conj()-1))),
        translation=float(np.max(np.abs(T@V-sgn*V@T))))
    assert max(errors.values())<1e-12
    checks.append(dict(L=L,parity=P,errors=errors))
pairings=[]
for L in [16,20,24]:
    d=json.loads((HERE/f'dispersion_L{L}.json').read_text());assert d['complete']
    lookup={(r['parity'],r['k']):r['E'] for r in d['sectors']}
    err=max(abs(E-lookup[P,((L//2 if P else 0)-k)%L]) for (P,k),E in lookup.items())
    assert err<1e-8
    pairings.append(dict(L=L,max_energy_pairing_error=err))
out=dict(matrix_checks=checks,spectral_checks=pairings)
(HERE/'symmetry_validation.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
