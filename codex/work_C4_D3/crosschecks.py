import json,numpy as np
from pathlib import Path
from analyze import load
from d3 import Sector,swap,local_components,projectors
rows=load();out={}
def get(family,label,L,source=None):
 rr=[r for r in rows if r['family']==family and r['label']==label and r['L']==L and (source is None or r['source']==source)]
 return rr[-1] if rr else None
for L in [6,9,12]:
 a=get('uls3','uls3',L);b=get('spin1','bb_t012',L)
 out[f'uls3_spin1_energy_L{L}']=abs((a['E0']+L)/np.sqrt(2)-b['E0'])
 out[f'uls3_spin1_entropy_L{L}']=float(np.max(np.abs(np.array(a['entropy'])-b['entropy'])))
for L in [8,10,12]:
 a=get('uls4','uls4',L)
 b=get('sp4','sp4_t024',L) if L<12 else get('sp4','sp4_angle1.57079632680',L)
 if b:out[f'uls4_sp4_energy_L{L}']=abs((a['E0']+L)/2-b['E0'])
 b=get('orbital','orb_x0.25_y0.25',L) if L<12 else get('orbital','orb_s1_x0.25_y0.25',L)
 if b:out[f'uls4_orbital_energy_L{L}']=abs(a['E0']/4-L/16-b['E0'])
errs=[]
for a in rows:
 if a['source']=='refine_sp4.jsonl' and a['L'] in [8,10]:
  for b in rows:
   if b['source']=='scan_sp4.jsonl' and b['L']==a['L'] and np.max(np.abs(np.array(b['coeff'])-a['coeff']))<1e-9:
    errs.append(abs(a['E0']-b['E0']))
out['Sp4_Cartan_vs_Sz_energy_max']=max(errs)
for family in ['tb1','tb32','su2','sp4']:
 d,hs,ch=local_components(family);L=4;target=(L*ch[0]).tolist();sec=Sector(family,L,target)
 ev,_,_=sec.solve([1]*len(hs),0,1)
 out[family+'_polarized_energy_error']=float(abs(ev[0]-L*sum(h[0,0] for h in hs)))
out['production_max_residual']=max(r['max_residual'] for r in rows)
assert all(v<1e-7 for k,v in out.items() if k!='production_max_residual'),out
assert out['production_max_residual']<1e-7,out
Path('crosschecks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
