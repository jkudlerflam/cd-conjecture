#!/usr/bin/env python3
"""Recompute screen summaries from raw records, including all null results."""
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def dump(name, data):
    (HERE/name).write_text(json.dumps(data,indent=2)+'\n')


def converged_stage(h):
    return (len(h)>=3 and abs(h[-1]['E']-h[-2]['E'])<=1e-8
            and abs(h[-1]['S_fractional']-h[-2]['S_fractional'])<=1e-3
            and h[-1]['res_max']<=1e-8)


def best_chi_check(runs, parity):
    checks=[]
    for r in runs:
        if r['parity']!=parity:continue
        stages={}
        for h in r['history']:
            if h['alpha']==0:stages.setdefault(h['chi'],[]).append(h)
        chis=sorted(stages)
        for a,b in zip(chis,chis[1:]):
            if converged_stage(stages[a]) and converged_stage(stages[b]):
                checks.append(dict(seed=r['seed'],chi_from=a,chi_to=b,
                    absolute_entropy_step=abs(stages[b][-1]['S_fractional']-stages[a][-1]['S_fractional'])))
    return max(checks,key=lambda x:(x['chi_to'],x['chi_from'])) if checks else None


def ring_summary():
    rows=[json.loads(s) for s in (HERE/'screen.jsonl').read_text().splitlines()]
    groups={}
    for r in rows:groups.setdefault(r['name'],{})[r['L']]=r
    out=[]
    for name,g in groups.items():
        a,b=g[10],g[12]
        selected=[next(s for s in r['sectors'] if s['parity']==r['selected_parity']) for r in (a,b)]
        # Window sensitivity is computed from preserved profiles, not new eigensolves.
        fits=[]
        for r,s in zip((a,b),selected):
            x=np.arange(1,r['L']);S=np.asarray(s['entropy_profile'])
            masks=[(x>=f*r['L'])&(x<=(1-f)*r['L']) for f in (.25,.375)]
            fits.append([float(3*np.polyfit(np.log(r['L']/np.pi*np.sin(np.pi*x[m]/r['L'])),S[m],1)[0]) for m in masks])
        out.append(dict(name=name,couplings=b['couplings'],metrics=b['metrics'],
            c_size=3*(b['entropy_half']-a['entropy_half'])/math.log(1.2),
            c_profiles=fits,densities=[a['density_mean'],b['density_mean']],
            parities=[a['selected_parity'],b['selected_parity']],
            translation_moduli=[s['translation_modulus'] for s in selected],
            density_ranges=[s['density_profile_range'] for s in selected],
            raw_lanczos_gaps=[a['ground_gap'],b['ground_gap']],
            gap_caveat='Single-vector Lanczos may omit exact degenerate partners; not a certified gap'))
    dump('stage_a_summary.json',out)
    return out


def dmrg_summary():
    rows=[json.loads(p.read_text()) for p in (HERE/'dmrg').glob('*.json')]
    groups={}
    for r in rows:
        if r['status']=='done':groups.setdefault(r['name'],{}).setdefault(r['L'],[]).append(r)
    out=[]
    for name,sz in sorted(groups.items()):
        states=[]
        for L,runs in sorted(sz.items()):
            best=min(runs,key=lambda r:r['E_physical'])
            # Multiple seeds can be present; compare lowest state within each parity.
            parity={P:min([r for r in runs if r['parity']==P],key=lambda r:r['E_physical']) for P in set(r['parity'] for r in runs)}
            parity_complete=set(parity)==({0} if name.startswith('control') else {0,1})
            p_gap=parity[1]['E_physical']-parity[0]['E_physical'] if len(parity)==2 else None
            stage_ends={}
            for h in best['history']:
                if h['alpha']==0:stage_ends[h['chi']]=h
            chis=sorted(stage_ends)
            chi_step=abs(stage_ends[chis[-1]]['S_fractional']-stage_ends[chis[-2]]['S_fractional']) if len(chis)>1 else None
            previous=[h for h in best['history'] if len(chis)>1 and h['chi']==chis[-2] and h['alpha']==0]
            previous_converged=(len(previous)>1 and abs(previous[-1]['E']-previous[-2]['E'])<=1e-8
                and abs(previous[-1]['S_fractional']-previous[-2]['S_fractional'])<=1e-3
                and previous[-1]['res_max']<=1e-8) if len(chis)>1 else None
            corr=best['measurements'].get('long_connected_correlations',{})
            states.append(dict(L=L,parity=best['parity'],seed=best['seed'],chi=best['chi'],
                parity_complete=parity_complete,ground_branch_status='selected' if parity_complete else 'provisional_missing_parity',
                energy=best['E_physical'],parity_gap=p_gap,S_central3=best['S_central3'],
                S_fractional=best['S_fractional'],density=best['density_bulk'],
                density_range=best['density_bulk_range'],profiles=best['profiles'],
                convergence=best['convergence'],chi_entropy_step=chi_step,
                preceding_chi_sweep_converged=previous_converged,
                best_available_same_seed_chi_check=best_chi_check(runs,best['parity']),
                correlations_at_L_over_2={k:v[-1][1] for k,v in corr.items() if k!='origin'},
                sector_states=[dict(parity=p,energy=r['E_physical'],S_fractional=r['S_fractional'],
                                    S_central3=r['S_central3'],sweep_converged=r['convergence'].get('sweep_converged')) for p,r in parity.items()]))
        pairs=[]
        for a,b in zip(states,states[1:]):
            if not a['parity_complete'] or not b['parity_complete']:
                continue
            factor=6/math.log(b['L']/a['L'])
            sector_slopes={}
            for pa in a['sector_states']:
                pb=next((x for x in b['sector_states'] if x['parity']==pa['parity']),None)
                if pb:sector_slopes[pa['parity']]={key:factor*(pb[key]-pa[key]) for key in ['S_central3','S_fractional']}
            pairs.append(dict(sizes=[a['L'],b['L']],c_central3=factor*(b['S_central3']-a['S_central3']),
                c_fractional=factor*(b['S_fractional']-a['S_fractional']),same_parity_slopes=sector_slopes,
                sweep_converged=all(s['convergence'].get('sweep_converged',False) for s in [a,b])))
        out.append(dict(name=name,states=states,pairs=pairs))
    dump('dmrg_summary.json',out)
    for r in out:
        print(r['name'])
        for p in r['pairs']:print(' ',p['sizes'],'c3/frac',round(p['c_central3'],5),round(p['c_fractional'],5),'sweeps',p['sweep_converged'])
        for s in r['states']:print(' ',s['L'],'P',s['parity'],'n',round(s['density'],5),'chi',s['chi'],'dSchi',s['chi_entropy_step'],'conv',s['convergence'])
    return out


if __name__=='__main__':
    ring_summary();dmrg_summary()
