"""Audit saved C19 records without repeating ED."""
import json
from pathlib import Path
import numpy as np
from uniform_cell import exact_filter_check

for n in [1,2,4,8,12]:assert exact_filter_check(n)
for name in ['checks_results.json','completion_results.json','reflection_chain_results.json','exchange_field_results.json']:
    data=json.loads(Path(name).read_text())
    rows=data if isinstance(data,list) else data.get('ed',[])+data.get('reflection_ed',[])
    for row in rows:
        assert row['eigenpair_residual']<3e-11,(name,row)
        if 'endpoint_leakage' in row: assert row['endpoint_leakage']<1e-11
        if 'heff_real' in row:
            h=np.array(row['heff_real'])+1j*np.array(row['heff_imag'])
            assert np.linalg.norm(h-h.conj().T)<1e-12
            assert np.max(abs(np.linalg.eigvalsh(h)-row['low_eigenvalues']))<1e-12
    if isinstance(data,dict):
        for row in data.get('algebra',[])+data.get('reflection',[]):
            assert row['bond_error']<1e-11
            assert row['kernel_error']<1e-11
            assert row['min_weight']>0
            assert row['parity_error']<1e-11
    print(name,': PASS',len(rows),'ED records')
f=json.loads(Path('exchange_field_results.json').read_text())
r=[x for x in f['ed'] if x['n']==1 and x['cells']==2]
ratio=r[0]['error']/r[1]['error']
assert 7.5<ratio<8.5
print('Final exchange-only two-cell error ratio:',ratio)
print('Exact filter identities and all saved records pass.')
