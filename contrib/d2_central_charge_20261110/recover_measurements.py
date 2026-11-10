"""Recover observables, not convergence certificates, from final saved states."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import hashlib,importlib.util,json,pickle
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('_cc_runner',HERE/'run.py');R=importlib.util.module_from_spec(spec);spec.loader.exec_module(R)
out=[]
for f in (HERE/'results').glob('*.json'):
    r=json.loads(f.read_text())
    if r['status']!='postprocessing_failed':continue
    ck=Path(r['checkpoint']);sha=hashlib.sha256(ck.read_bytes()).hexdigest()
    with ck.open('rb') as h:s=pickle.load(h)
    assert s['checkpoint_kind']=='final'
    for k in ['L','geometry','parity','couplings']:assert s[k]==r[k]
    psi=s['psi'];p=r['couplings'];L=r['L'];inf=r['geometry']=='infinite'
    model=R.PauliChain(dict(L=L,bc_MPS='infinite' if inf else 'finite',bc_x='periodic',couplings=p,conserve=None if inf else 'parity'))
    E=float(R.local_energy(psi,p).real.mean()+model.energy_offset/L) if inf else float(model.H_MPO.expectation_value(psi)+model.energy_offset)
    ent={str(n):psi.entanglement_entropy(n=n).tolist() for n in [1,2,3]}
    row=dict(source=str(f.relative_to(HERE)),source_sha256=hashlib.sha256(f.read_bytes()).hexdigest(),
        checkpoint=str(ck),checkpoint_sha256=sha,E_physical=E,entropy=ent,
        status='observables_recovered_no_sweep_convergence_certificate',sweeps=s['sweeps'])
    if not inf:row['entropy_stats']={k:R.stats(v,L,6*int(k)/(int(k)+1)) for k,v in ent.items()}
    assert sha==hashlib.sha256(ck.read_bytes()).hexdigest();out.append(row)
(HERE/'recovered_measurements.json').write_text(json.dumps(out,indent=2)+'\n')
for r in out:print(r['source'],r['E_physical'],r.get('entropy_stats',{}).get('1'))
