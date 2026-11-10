"""Reproduce comparisons with pinned ED/OBC results and exact on-site gauge."""
import json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
r=json.loads((HERE/'results/candidate_periodic_L16_P0_s17_chi128.json').read_text())
reference=next(json.loads(x) for x in (HERE/'inputs/bridge_rings.jsonl').read_text().splitlines() if json.loads(x)['L']==16)
sector=next(x for x in reference['sectors'] if x['parity']==0)
eerr=abs(r['E_physical']-sector['energies'][0]);serr=float(np.max(np.abs(np.array(r['entropy']['1'])-sector['raw']['entropy_profile'])))
assert eerr<1e-8 and serr<1e-7
(HERE/'validation_L16.json').write_text(json.dumps(dict(energy_error=eerr,entropy_profile_error=serr,reference='inputs/bridge_rings.jsonl'),indent=2)+'\n')
a,b=[json.loads((HERE/f'results/candidate_periodic_L8_P0_s17_chi16{suffix}.json').read_text()) for suffix in ['', '_g']]
eerr=abs(a['E_physical']-b['E_physical']);serr=float(np.max(np.abs(np.array(a['entropy']['1'])-b['entropy']['1'])))
assert eerr<1e-10 and serr<1e-10
(HERE/'gauge_validation.json').write_text(json.dumps(dict(energy_error=eerr,entropy_profile_error=serr),indent=2)+'\n')
f=HERE/'results/candidate_open_L48_P0_s17_chi128_g.json'
r=json.loads(f.read_text()) if f.exists() else {'status':'pending'}
if r['status']=='done':
    old=json.loads((HERE/'inputs/size_comparisons.json').read_text())['b_r5_2_d1_10_q0_hm7_8']
    ref=next(x for x in old['ground_branches'] if x['L']==48)
    result=dict(new_source=str(f.relative_to(HERE)),previous_source=ref['source'],
        previous_record_sha256=ref['record_sha256'],new_converged=r['converged'],
        energy_error=abs(r['E_physical']-ref['E_physical']),
        fractional_entropy_error=abs(r['entropy_stats']['1']['fractional']-ref['S_fractional']),
        central3_entropy_error=abs(r['entropy_stats']['1']['central3']-ref['S_central3']))
    result['agreement']=result['energy_error']<1e-7 and max(result['fractional_entropy_error'],result['central3_entropy_error'])<1e-5
    assert result['agreement'],result
    (HERE/'obc_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
else:print('OBC check pending; ED and gauge comparisons pass')
