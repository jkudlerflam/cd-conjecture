#!/usr/bin/env python3
import json, numpy as np
from pathlib import Path
FILES=['calibration.jsonl','calibration4.jsonl','calibration_tb1.jsonl','calibration_tb32.jsonl','calibration_more.jsonl','calibration_tb1_14.jsonl','scan_sp4.jsonl','scan_su2.jsonl','scan_spin1.jsonl','scan_orbital.jsonl','scan_orbital_dense.jsonl','scan_u1.jsonl','corrected_spin1_L9.jsonl','refine_sp4.jsonl','corrected_u1_ground.jsonl','refine_su2.jsonl','refine_orbital.jsonl','refine_spin1.jsonl','refine_u1.jsonl','scan_su2near.jsonl','unequal_velocity_control.jsonl','spin_diagnostics.jsonl','sp4_L4.jsonl','orbital_L4.jsonl','saturation.jsonl','sp4_fine.jsonl']
def load():
 out=[]
 for fn in FILES:
  if Path(fn).exists():
   for line in Path(fn).read_text().splitlines():
    try:r=json.loads(line)
    except json.JSONDecodeError:continue
    r['source']=fn;out.append(r)
 corrected={(r['label'],r['L']) for r in out if r['source']=='corrected_u1_ground.jsonl'}
 out=[r for r in out if (r['label'],r['L']) not in corrected or r['source']=='corrected_u1_ground.jsonl']
 return out

def cf(r,cut=2):
 f=r.get('entropy_fits',[None]*3)[cut-1]
 return f['c'] if f else np.nan

def fit_energy(rs):
 rs=sorted(rs,key=lambda r:r['L']);L=np.array([r['L'] for r in rs],float);e=np.array([r['E0']/r['L'] for r in rs]);v=np.array([r['v_current_assumed'] for r in rs])
 x=L**-2
 a=np.polyfit(x,e,1);vf=np.polyfit(x,v,1)[1]
 out={'sizes':L.tolist(),'A':float(-a[0]),'v_extrap':float(vf),'cE_linear':float(-6*a[0]/(np.pi*vf)),'cE_last_v':float(-6*a[0]/(np.pi*v[-1]))}
 if len(rs)>2:
  a2=np.polyfit(x,e,2);out['cE_quadratic']=float(-6*a2[1]/(np.pi*vf))
 return out
if __name__=='__main__':
 rows=load()
 for family in ['uls3','uls4','tb1','tb32']:
  rs=[r for r in rows if r['family']==family and r['source'].startswith('calibration')]
  valid=[r for r in rs if r['L']>=6 and (family!='uls3' or r['L']%3==0) and (family!='uls4' or r['L']%4==0)]
  if family=='uls4':valid=[r for r in rs if r['L']%4==0]
  print('CAL',family,[(r['L'],round(cf(r),6),round(r['v_current_assumed'],6)) for r in rs])
  for n in range(2,len(valid)+1):print(fit_energy(valid[-n:]))
 for family in ['sp4','su2','spin1','orbital']:
  rs=[r for r in rows if r['family']==family];last=max(r['L'] for r in rs)
  top=sorted((r for r in rs if r['L']==last),key=lambda r:cf(r),reverse=True)[:20]
  print('\nTOP',family,'count',len(rs),'L',last)
  for r in top:
   others=sorted([s for s in rs if s['label']==r['label']],key=lambda s:s['L'])
   print(r['label'],[(s['L'],round(cf(s),4),round(s['v_current_assumed'] or 0,5),s['ground_k']) for s in others])
