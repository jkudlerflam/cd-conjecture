import json,numpy as np
from pathlib import Path
from analyze import load
rows=[json.loads(x) for x in open('sp4_towers.jsonl')];res=[]
for t in [2*np.pi/3,np.pi]:
 rs=sorted([r for r in rows if abs(r['theta']-t)<1e-9],key=lambda r:r['L'])
 if len(rs)!=3:continue
 L=np.array([r['L'] for r in rs],float);xx=L**-2
 va=np.array([r['adjoint_current_velocity'] for r in rs]);vm=np.array([r['vector_mixed_current_velocity'] for r in rs])
 if any(v is None for v in va) or any(v is None for v in vm):continue
 einf=np.array([r['E0']/r['L'] for r in rs]);A=-np.polyfit(xx,einf,2)[1]
 va0=np.polyfit(xx,va,1)[1];vm0=np.polyfit(xx,vm,1)[1];v10=2*vm0-va0
 predicted=np.pi*(2*va0+vm0)/6
 row={'theta':t,'sizes':L.tolist(),'adjoint_velocities':va.tolist(),'mixed_vector_velocities':vm.tolist(),'v5_extrap':va0,'v1_extrap':v10,'A_energy':A,'A_two_velocity_prediction':predicted,'relative_difference':A/predicted-1,'naive_single_velocity_c':6*A/(np.pi*va0),'c_total_under_five_plus_one_Majorana_ansatz':3,'ansatz_is_not_a_lattice_proof':True}
 res.append(row)
Path('sp4_two_velocity_test.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
