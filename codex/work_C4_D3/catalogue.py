#!/usr/bin/env python3
import json,csv,math,numpy as np
from pathlib import Path
from collections import Counter
from analyze import load,cf,fit_energy

def classify(r):
 f=r['family'];c=np.array(r['coeff']);eps=2e-9
 if f in ('uls3','uls4','tb1','tb32'):return 'supplied_CFT_calibration'
 if f=='sp4':
  a,b=c
  if b<=eps and b<=a+eps:return 'proved_FM_or_degenerate_nonCFT'
  if abs(a)<eps and b>eps:return 'ULS_SU4_anchor'
 if f=='su2':
  if c[3]<=min(c)+eps:return 'proved_FM_or_degenerate_nonCFT'
  if abs(c[0]-c[2])<eps and abs(c[1]-c[3])<eps and c[1]>c[0]+eps:return 'ULS_SU4_anchor'
 if f=='spin1' and abs(c[0]-1)<eps and abs(c[1]-1)<eps and c[3]>=4-eps and c[3]>=1+abs(c[2]+1)-eps:
  return 'proved_vacuum_or_quadratic_saturation'
 if f=='spin1' and abs(c[2])+abs(c[3])<eps:
  a,b=c[:2];e=np.array([-2*a+4*b,-a+b,a+b])
  if e[2]<=min(e)+eps:return 'proved_FM_or_degenerate_nonCFT'
  if abs(a-b)<eps and a>eps:return 'ULS_SU3_anchor'
 if f=='orbital':
  g,a,b=c
  if g>eps and abs(a-.75*g)<eps and abs(b-.75*g)<eps:return 'proved_dimer_ground_states_gap_unproved'
  energies=[9*g/16-3*(a+b)/4,-3*g/16-3*a/4+b/4,-3*g/16+a/4-3*b/4,g/16+(a+b)/4]
  if energies[3]<=min(energies)+eps:return 'proved_polarized_or_degenerate_nonCFT'
  if g>eps and abs(a-g/4)<eps and abs(b-g/4)<eps:return 'ULS_SU4_anchor'
  if abs(g)<eps and a>eps and b>eps:return 'decoupled_two_chain_control'
 if r.get('adjacent_gap') is not None and r['adjacent_gap']<1e-8:return 'degenerate_velocity_unusable'
 if not np.isfinite(cf(r)):return 'no_entropy_estimate'
 if cf(r)>r.get('d',3 if f=='spin1' else 4)-1:return 'entropy_excess_unresolved'
 return 'no_entropy_excess_at_this_size'

def fit_envelope(rs):
 a=fit_energy(sorted(rs,key=lambda r:r['L'])[-2:])['cE_linear'];b=fit_energy(sorted(rs,key=lambda r:r['L'])[-3:])['cE_quadratic']
 return {'midpoint':(a+b)/2,'half_width':abs(a-b)/2,'endpoints':[min(a,b),max(a,b)],'definition':'two-largest linear versus three-largest quadratic E/L in L^-2; velocity linear in L^-2 on each window'}

def main():
 rows=load();vel={}
 for fn in ['u1_global_velocity.jsonl','refined_global_velocity.jsonl','deformation_global_velocity.jsonl']:
  p=Path(fn)
  if p.exists():
   for line in p.read_text().splitlines():
    try:r=json.loads(line)
    except json.JSONDecodeError:continue
    vel[(r['label'],r['L'])]=r
 fields=['family','label','L','coefficients','charge','source','classification','E0','ground_k','dimension','c_ent2','regression_se','cut_systematic','cut1_c','cut3_c','velocity_ground_sector','velocity_global_charge_audit','adjacent_gap','ground_charge_set','max_residual']
 with open('catalogue.csv','w') as f:
  w=csv.DictWriter(f,fields);w.writeheader()
  for r in rows:
   fits=r.get('entropy_fits',[None]*3);center=cf(r);vals=[x['c'] for x in fits if x]
   v=vel.get((r['label'],r['L']),{})
   if not v and 'global_momentum_minima' in r:
    def fold(k):k%=r['L'];return min(k,r['L']-k)
    gap=min(r['global_momentum_minima'][str(fold(r['ground_k']+q))]['E'] for q in [-1,1])-r['E0']
    v={'global_v_current_assumed':gap*r['L']/(2*np.pi),'ground_charge_set_nonnegative':r.get('ground_charge_set_nonnegative',[])}
   row={'family':r['family'],'label':r['label'],'L':r['L'],'coefficients':json.dumps(r['coeff']),'charge':json.dumps(r['charge']),'source':r['source'],'classification':classify(r),'E0':r['E0'],'ground_k':r['ground_k'],'dimension':r['sector_dimension'],'c_ent2':center,'regression_se':fits[1]['regression_se'] if fits[1] else None,'cut_systematic':max(abs(x-center) for x in vals) if np.isfinite(center) else None,'cut1_c':fits[0]['c'] if fits[0] else None,'cut3_c':fits[2]['c'] if fits[2] else None,'velocity_ground_sector':r['v_current_assumed'],'velocity_global_charge_audit':v.get('global_v_current_assumed'),'adjacent_gap':r['adjacent_gap'],'ground_charge_set':json.dumps(v.get('ground_charge_set_nonnegative',[])),'max_residual':r['max_residual']}
   w.writerow(row)
 cal={}
 for family,expected in [('uls3',2),('uls4',3),('tb1',1.5),('tb32',1.8)]:
  rs=[r for r in rows if r['family']==family]
  rs=[r for r in rs if (not family.startswith('uls')) or r['L']%int(family[-1])==0]
  rs=sorted(rs,key=lambda r:r['L']);last=rs[-1]
  q=last['L']//(int(family[-1]) if family.startswith('uls') else 2);k=(last['ground_k']+q)%last['L'];k=min(k,last['L']-k)
  cal[family]={'expected_c':expected,'energy_fit_envelope':fit_envelope(rs),'largest_L':last['L'],'entropy_fits':last['entropy_fits'],'v_L':last['v_current_assumed'],'observed_soft_x':last['tower_x_assumed'][str(k)][0],'expected_soft_x':(int(family[-1])-1)/int(family[-1]) if family.startswith('uls') else (3/8 if family=='tb1' else 3/10)}
 cal['unequal_velocity_control']=fit_envelope([r for r in rows if r['label']=='orb_decoupled_ratio3'])
 Path('calibrations_summary.json').write_text(json.dumps(cal,indent=2)+'\n')
 counts={'records':len(rows),'by_source':dict(Counter(r['source'] for r in rows)),'by_classification':dict(Counter(classify(r) for r in rows)),'max_residual':max(r['max_residual'] for r in rows)}
 Path('campaign_counts.json').write_text(json.dumps(counts,indent=2)+'\n');print(json.dumps(counts,indent=2))

if __name__=='__main__':main()
