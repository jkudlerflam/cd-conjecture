import json,math
from pathlib import Path
from analyze import load
rows=load();tests={}
expected={'scan_sp4.jsonl':288,'scan_su2.jsonl':576,'scan_orbital_dense.jsonl':546,'scan_u1.jsonl':300,'refine_sp4.jsonl':69,'refine_su2.jsonl':28,'refine_orbital.jsonl':13,'refine_u1.jsonl':20,'saturation.jsonl':48,'sp4_fine.jsonl':768}
for fn,n in expected.items():
 got=len(Path(fn).read_text().splitlines());assert got==n,(fn,got,n);tests[fn+'_rows']=got
for fn,n in [('u1_global_velocity.jsonl',300),('refined_global_velocity.jsonl',20),('deformation_global_velocity.jsonl',72)]:
 rr=[json.loads(x) for x in open(fn)];keys={(r['label'],r['L']) for r in rr}
 assert len(keys)==n,(fn,len(keys),n)
 assert min(r['global_adjacent_gap'] for r in rr)>-1e-7,fn
 tests[fn+'_distinct_audits']=len(keys)
assert len(rows)==3526,len(rows)
assert all(r['family'] not in ['spin1','tb1'] or r['L']%2==0 or r['charge']==[0] or r['source']=='corrected_u1_ground.jsonl' for r in rows if not r['source'].startswith('scan_u1'))
assert all(math.isfinite(r['E0']) for r in rows)
assert max(r['max_residual'] for r in rows)<1e-7
# Verify the local-energy bound at the exact product-state saturation control.
sat=[r for r in rows if r['source']=='saturation.jsonl' and r['coeff'][3]>=4]
tests['saturation_energy_max_error']=max(abs(r['E0']-2*r['L']) for r in sat)
assert tests['saturation_energy_max_error']<1e-7
for label in ['orb_s1_x0.75_y0.75']:
 rr=[r for r in rows if r['label']==label]
 tests['dimer_product_energy_max_error']=max(abs(r['E0']+9*r['L']/16) for r in rr)
 assert tests['dimer_product_energy_max_error']<1e-7
tests['retained_records']=len(rows);tests['passed']=True
Path('verification.json').write_text(json.dumps(tests,indent=2)+'\n');print(json.dumps(tests,indent=2))
