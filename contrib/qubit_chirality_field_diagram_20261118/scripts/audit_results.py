"""Cross-check independently obtained operators, spectra and normalizations."""
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parents[1]


def main():
    records=[]
    for name in ('coarse.json','refine.json','larger.json'):
        records+=json.loads((HERE/'results'/name).read_text())['records']
    errors=dict(energy_reconstruction=0.,structure_sum_rule=0.,
                structure_positivity_violation=0.,spectral_operator_match=0.,
                projected_energy_match=0.,eigenpair_residual=0.)
    for r in records:
        errors['energy_reconstruction']=max(errors['energy_reconstruction'],
            abs(r['e']-r['observables']['bond']['mean']+r['h']*r['magnetization']))
        for obs in r['observables'].values():
            sf=np.array(obs['structure']);c=np.array(obs['connected_re'])
            errors['structure_sum_rule']=max(errors['structure_sum_rule'],abs(sf.mean()-c[0]))
            errors['structure_positivity_violation']=max(errors['structure_positivity_violation'],-min(sf))
        errors['eigenpair_residual']=max(errors['eigenpair_residual'],max(s['residual'] for s in r['sectors']))
    count=0
    for name in ('operator_levels.json','operator_levels20.json'):
        for r in json.loads((HERE/'results'/name).read_text())['records']:
            ref=next(t for t in records if t['L']==r['L'] and t['h']==r['h'])
            errors['projected_energy_match']=max(errors['projected_energy_match'],abs(r['ground_energy']-ref['E']))
            errors['eigenpair_residual']=max(errors['eigenpair_residual'],r['ground_residual'],max(b['residual'] for b in r['blocks']))
            block=next(b for b in r['blocks'] if b['shift']==r['L']//2)
            for name in ('Z','bond','kappa','B'):
                target=ref['observables'][name]['pi_over_L']*r['L']
                errors['spectral_operator_match']=max(errors['spectral_operator_match'],abs(target-block['operators'][name]['total_weight']))
            for block in r['blocks']:
                for obs in block['operators'].values():
                    assert sum(obs['low_level_weights'])<=obs['total_weight']+1e-8
            for obs in r['transverse'].values():
                assert abs(obs['connected_re'][0]-1)<1e-8
                assert min(obs['structure'])>-1e-8
            count+=1
    assert max(errors.values())<1e-8,errors
    result=dict(passed=True,ED_points=len(records),operator_points=count,max_errors=errors,
                scope='Finite matrix and contraction checks only; no thermodynamic certification.')
    (HERE/'results'/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
