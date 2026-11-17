#!/usr/bin/env python3
"""Short-field tolerance, independent start and the old V=4 chi gap."""
import json
import numpy as np
from common import HERE,ROOT,Ana


def main():
    a=json.loads((HERE/'analysis.json').read_text());rows=a['rows'];out={}
    groups={}
    for r in rows:
        if r['tag'].startswith('tolerance_'):
            groups.setdefault((r['a'],r['kappa']),[]).append(r)
    toler=[]
    for (aa,k),rs in sorted(groups.items()):
        if len(rs)!=2:continue
        r=min(rs,key=lambda r:r['E'])
        toler.append(dict(a=aa,kappa=k,delta_a=aa-.96,delta_kappa=k-1,delta_hz=(aa-.96)/8,delta_D=(k-1)/32,
            ground_P=r['P'],L=r['L'],chi=r['chi'],n_mid=r['n_mid'],Sbar=r['Sbar'],
            density_range=r['central_density_range'],valley=r['valley'],valley_half_boundary=r['valley_half_boundary'],
            both_gates=all(x['all_gates'] for x in rs),schedule_limit_hits=[x['schedule_limit_hit'] for x in rs],
            source_tags=[x['tag'] for x in rs]))
    out['tolerance_pilots']=toler
    out['tolerance_scope']='L=32 finite-state response only; no verified thermodynamic CFT window width.'
    rr=next((r for r in rows if r['tag']=='restart_V8_a24o25_L64_P0_seed7'),None)
    base=next((r for r in rows if r['tag']=='extended_V8_a24o25_L64_P0'),None)
    out['independent_start']=dict(status='pending') if rr is None else dict(status='done',
        tags=[base['tag'],rr['tag']],dE=rr['E']-base['E'],
        dSbar={str(w):rr['Sbar'][str(w)]-base['Sbar'][str(w)] for w in (1,6,12)},
        dn=rr['n_mid']-base['n_mid'],both_gates=rr['all_gates'] and base['all_gates'],
        schedule_limit_hits=[base['schedule_limit_hit'],rr['schedule_limit_hit']])
    legacy=[]
    for P in (0,1):
        new=next((r for r in rows if r['tag']==f'legacy_V4_a1_L96_P{P}'),None)
        if new is None:continue
        oldpath=ROOT/f'contrib/d2_order_one_20261110/dmrg/bridge_r4_u1_k4_L96_P{P}_s17.json'
        old=json.loads(oldpath.read_text());raw=json.loads((HERE/new['file']).read_text());stage192=[r for r in raw['hist'] if r['chi']==192][-1]
        legacy.append(dict(P=P,old_source=str(oldpath.relative_to(ROOT)),new_source=new['file'],
            old=dict(chi=old['chi'],E=old['E_physical'],Sbar={str(w):Ana.sbar(old['Sprof'],96,w) for w in (1,6,12)},density_range=old['density_bulk_range']),
            new192=dict(E=stage192['E_phys'],Sbar6_sweep=stage192['Sbar6']),
            new256=dict(E=new['E'],Sbar=new['Sbar'],density=new['n_mid'],density_range=new['central_density_range'],
                gates=new['gates'],schedule_limit_hit=new['schedule_limit_hit'],dS_chi=new['dS_chi'],variance=raw.get('variance'),variance_per_site=raw['variance']/96 if 'variance' in raw else None,
                physical_energy_check=raw.get('variance_energy'))))
    out['legacy_V4']=legacy
    out['legacy_scope']='Same microscopic V=4,kappa=1,a=1; new seed and converged-stage chi comparison. No new CFT fit from this alone.'
    (HERE/'followup_analysis.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
