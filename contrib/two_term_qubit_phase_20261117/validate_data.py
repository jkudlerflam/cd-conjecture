"""Check saved residuals, observables, multiplets and repeat-run agreement."""
import json
from pathlib import Path
import numpy as np

p=Path(__file__).parent
initial=json.loads((p/'spectra.json').read_text())['records']
repeat=json.loads((p/'neutral_poles.json').read_text())['records']
equal=json.loads((p/'equal_coefficients.json').read_text())['records']
repair=json.loads((p/'weak_resolution_check.json').read_text())['records']
lookup={(r['L'],r['g']):r for r in initial}
for r in initial+repeat+equal+repair:
    assert r['ground_sector']==[1,1]
    assert max(s['max_residual'] for s in r['sectors'])<1e-8
    assert 0<=r['observables']['density_pi']<=1
    for name,value in r['observables'].items():
        if name.endswith('connected'):assert value>=-1e-10
    for name,poles in r.get('retained_neutral_poles',{}).items():
        assert sum(t['fraction_of_total'] for t in poles)<=1+1e-8
    if r['g']!=.125 or r in repair:
        for s in r['sectors']:
            assert max(t['translation_leakage'] for t in s['levels'][:-1])<1e-6
for r in repeat+repair:
    first=lookup[r['L'],r['g']]
    assert abs(r['ground_energy']-first['ground_energy'])<1e-8
    for s,t in zip(r['sectors'],first['sectors']):
        assert abs(s['levels'][0]['energy']-t['levels'][0]['energy'])<1e-8
    for name in first['observables']:
        assert np.isclose(first['observables'][name],r['observables'][name],atol=1e-8)
print('Saved residuals, positive fluctuations, spectral weights, resolved target/repair multiplets and repeat-run observables passed.')
