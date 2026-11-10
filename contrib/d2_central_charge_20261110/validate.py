"""Validate periodic DMRG against independent full-Hilbert-space ED."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import sys,json,numpy as np
from pathlib import Path
from scipy.sparse.linalg import eigsh
HERE=Path(__file__).resolve().parent
from independent_ed import ed
rows=[]
for P in (0,1):
    f=HERE/'results'/f'candidate_periodic_L8_P{P}_s17_chi16.json';r=json.loads(f.read_text())
    H,b=ed.build_sector(8,{k:v for k,v in r['couplings'].items() if k!='Dy'},parity=P)
    vals,vec=eigsh(H,k=4,which='SA',tol=1e-13,v0=np.random.default_rng(42).normal(size=H.shape[0]))
    order=np.argsort(vals);vals=vals[order];vec=vec[:,order]
    full=np.zeros(256,complex);full[b]=vec[:,0]
    S=ed.entropy_profile(full,8);gap=float(vals[1]-vals[0])
    row=dict(parity=P,energy_error=float(abs(vals.min()-r['E_physical'])),gap=gap,
        entropy_comparison_valid=gap>1e-8,
        entropy_error=float(np.max(np.abs(S-r['entropy']['1']))) if gap>1e-8 else None)
    assert row['energy_error']<1e-9
    assert row['entropy_error'] is None or row['entropy_error']<1e-7
    rows.append(row)
result=dict(checks=rows,note='Odd L8 ring ground state is degenerate. An initial unconditional entropy comparison failed for this reason; individual degenerate eigenvectors need not have equal interval entropy. Its energy is checked and entropy equality is not asserted.')
(HERE/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
