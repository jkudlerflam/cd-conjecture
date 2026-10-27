"""Gauge sensitivity of the Fourier-density estimator, same total H."""
import json
from pathlib import Path
from d2 import detailed, CALIBRATIONS, dump
out=json.loads(Path('results/density_check.json').read_text()) if Path('results/density_check.json').exists() else {}
for name in ('Ising','XXZ_0.5'):
    p=CALIBRATIONS[name]; out.setdefault(name,[])
    for L in (12,16,20,22):
        if any(r['L']==L for r in out[name]): continue
        row=detailed(p['J'],p['field'],L,m=0 if name=='Ising' or L%4==0 else L//2,levels=16,nup=None if name=='Ising' else L//2,field_split=1.)
        out[name].append(row)
        print(name,L,row['stress'][0]['selected']['c_MV'],flush=True)
        dump('results/density_check.json',out)
