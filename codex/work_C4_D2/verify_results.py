"""Cross-check independent full-space and orbit-space runs of the same models."""
import json
from pathlib import Path
from d2 import dump
raw={}
for name in ('global_scan','boundary_scan','optimization','optimization16','shortlist','global_shortlist'):
    path=Path('results/'+name+'.json')
    if not path.exists(): continue
    for p in json.loads(path.read_text()):
        for r in p.get('rows',[]): raw[(str(p['id']),r['L'])]=r
sources=['best_opt','candidate_refined','candidate854','candidate628','candidate871',
         'extra_candidates','extra_candidates_b','extra_candidates_c','structured_refined',
         'critical_refined','optimized16_refined','max_raw_refined','ferromagnet_half','structured_dm_refined']
energy_errors=[]; entropy_errors=[]; residuals=[]; sizes=[]; targets=[]
for name in sources:
    path=Path('results/'+name+'.json')
    if not path.exists(): continue
    for key,p in json.loads(path.read_text()).items():
        for r in p['rows']:
            residuals.append(r['residual']); sizes.append(r['L'])
            for s in r['stress']:
                if 'target_eigen_residual' in s: targets.append(s['target_eigen_residual'])
            a=raw.get((str(key),r['L']))
            if a is None: continue
            energy_errors.append(abs(a['E0']-r['E0']))
            if a.get('ground_degeneracy_returned')==1 and a.get('gap',0)>1e-6:
                entropy_errors.append(max(abs(x-y) for x,y in zip(a['entropy'],r['entropy'])))
result=dict(full_vs_orbit_energy_checks=len(energy_errors),
            max_full_vs_orbit_energy_error=max(energy_errors),
            unique_ground_entropy_checks=len(entropy_errors),
            max_full_vs_orbit_entropy_error=max(entropy_errors),
            max_ground_eigen_residual=max(residuals),
            max_stored_target_eigen_residual=max(targets,default=None),largest_L=max(sizes))
ff=json.loads(Path('results/free_fermion_control.json').read_text())
ff_errors=[abs(a['E0']-b['E0']) for item in ff.values() for a,b in zip(item['exact'],item['rows'])]
result['XX_exact_energy_checks']=len(ff_errors)
result['max_XX_exact_energy_error']=max(ff_errors)
assert result['max_full_vs_orbit_energy_error']<1e-8,result
assert result['max_full_vs_orbit_entropy_error']<1e-6,result
assert result['max_ground_eigen_residual']<1e-7,result
assert result['max_XX_exact_energy_error']<1e-8,result
dump('results/verification.json',result)
print(result)
