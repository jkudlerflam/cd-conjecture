"""Convergence-gated periodic size and infinite entanglement summaries."""
import hashlib,json,math
from collections import defaultdict
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
rows=[]
for f in sorted((HERE/'results').glob('*.json')):
    r=json.loads(f.read_text());r['source']=str(f.relative_to(HERE))
    r['source_sha256']=hashlib.sha256(f.read_bytes()).hexdigest()
    if r['status']=='done' and 'final_sweep_norm_small' not in r['convergence_checks']:
        r['convergence_checks']['final_sweep_norm_small']=r['sweep_stats']['norm_err'][-1]<1e-5
        r['converged']=all(r['convergence_checks'].values())
    rows.append(r)
inventory=[{k:r.get(k) for k in ['source','status','geometry','model','L','parity','seed','chi','gauge','converged','E_physical','convergence_checks']} for r in rows]
periodic=defaultdict(list);inf=defaultdict(list)
for r in rows:
    if r['status']!='done':continue
    if r['geometry']=='periodic':periodic[(r['model'],r['L'],r['parity'])].append(r)
    elif r['geometry']=='infinite':inf[(r['model'],r['L'],r['seed'],r.get('gauge',False),r.get('config',{}).get('one_site',False))].append(r)

def brief(r):
    out={k:r[k] for k in ['source','source_sha256','model','L','parity','seed','chi','E_physical','converged']}
    out['density_mean']=r['density_mean'];out['density_range']=r['density_range']
    out['gauge']=r.get('gauge',False);out['one_site']=r['config'].get('one_site',False)
    if r['geometry']!='infinite':out['entropy_stats']=r['entropy_stats']
    else:
        out['S']=float(np.mean(r['entropy']['1']));out['xi']=r['correlation_lengths'][0]
        out['S_renyi']={alpha:float(np.mean(r['entropy'][alpha])) for alpha in ['1','2','3']}
        out['correlation_lengths']=r['correlation_lengths']
        out['X']=r['X'];out['Y']=r['Y']
    return out

# Energy comparisons do not require a stationary choice of vector within a
# degenerate eigenspace. All other original gates remain mandatory. These
# rows cannot enter an entropy estimator unless the entropy gate also passes.
energy_keys=['recent_energy_stable','local_residual_small','norm_before_small',
             'final_sweep_norm_small','energy_agreement']
energy_sectors=[]
for key,rr in sorted(periodic.items()):
    good=[r for r in rr if all(r['convergence_checks'][k] for k in energy_keys)]
    if good:energy_sectors.append(brief(min(good,key=lambda r:r['E_physical'])))
energy_byL=defaultdict(dict)
for r in energy_sectors:energy_byL[(r['model'],r['L'])][r['parity']]=r
parity_gaps=[]
for key,sectors in sorted(energy_byL.items()):
    if set(sectors)=={0,1}:
        gap=sectors[1]['E_physical']-sectors[0]['E_physical']
        parity_gaps.append(dict(model=key[0],L=key[1],gap=gap,L_gap=key[1]*gap,
            sources=[sectors[p]['source'] for p in [0,1]],
            chis=[sectors[p]['chi'] for p in [0,1]],
            both_entropy_gates_pass=all(sectors[p]['converged'] for p in [0,1])))
bond_checks=[];gap_by_chi=defaultdict(dict)
for key,rr in sorted(periodic.items()):
    accepted=[]
    for chi in sorted(set(r['chi'] for r in rr)):
        energy_good=[r for r in rr if r['chi']==chi and all(r['convergence_checks'][k] for k in energy_keys)]
        if energy_good:
            winner=min(energy_good,key=lambda r:r['E_physical'])
            gap_by_chi[(key[0],key[1],chi)][key[2]]=brief(winner)
        full_good=[r for r in rr if r['chi']==chi and r['converged']]
        if full_good:accepted.append(brief(min(full_good,key=lambda r:r['E_physical'])))
    for a,b in zip(accepted,accepted[1:]):
        bond_checks.append(dict(model=key[0],L=key[1],parity=key[2],chis=[a['chi'],b['chi']],
            sources=[a['source'],b['source']],delta_E=b['E_physical']-a['E_physical'],
            delta_S={alpha:{w:b['entropy_stats'][alpha][w]-a['entropy_stats'][alpha][w]
                for w in ['central3','fractional']} for alpha in ['1','2','3']}))
matched_gaps=[]
for key,sectors in sorted(gap_by_chi.items()):
    if set(sectors)=={0,1}:
        matched_gaps.append(dict(model=key[0],L=key[1],chi=key[2],
            gap=sectors[1]['E_physical']-sectors[0]['E_physical'],sources=[sectors[p]['source'] for p in [0,1]]))

sector=[]
for key,rr in sorted(periodic.items()):
    good=[r for r in rr if r['converged']]
    if good:sector.append(brief(min(good,key=lambda r:r['E_physical'])))
byL=defaultdict(list)
for r in sector:byL[(r['model'],r['L'])].append(r)
selected=[]
selection_exclusions=[]
for key,rr in sorted(byL.items()):
    chosen=min(rr,key=lambda r:r['E_physical']).copy()
    chosen['both_parities_converged']=set(r['parity'] for r in rr)=={0,1}
    energy_candidates=energy_byL.get(key,{})
    chosen['both_parity_energy_gates_pass']=set(energy_candidates)=={0,1}
    if energy_candidates:
        energy_winner=min(energy_candidates.values(),key=lambda r:r['E_physical'])
        if energy_winner['parity']!=chosen['parity']:
            selection_exclusions.append(dict(model=key[0],L=key[1],reason='Lower energy-gated parity has no accepted entropy state',source=energy_winner['source']))
            continue
    selected.append(chosen)
pairs=[]
for model in sorted(set(r['model'] for r in selected)):
    rr=sorted([r for r in selected if r['model']==model],key=lambda r:r['L'])
    for a,b in zip(rr,rr[1:]):
        pairs.append(dict(model=model,L=[a['L'],b['L']],parities=[a['parity'],b['parity']],
            both_parities_converged=a['both_parities_converged'] and b['both_parities_converged'],
            both_parity_energy_gates_pass=a['both_parity_energy_gates_pass'] and b['both_parity_energy_gates_pass'],
            c={alpha:{w:6*float(alpha)/(float(alpha)+1)*(b['entropy_stats'][alpha][w]-a['entropy_stats'][alpha][w])/math.log(b['L']/a['L'])
                for w in ['central3','fractional']} for alpha in ['1','2','3']}))
fes=[]
for key,rr in sorted(inf.items()):
    good=sorted([brief(r) for r in rr if r['converged']],key=lambda r:r['chi'])
    entry=dict(model=key[0],cell=key[1],seed=key[2],gauge=key[3],one_site=key[4],points=good,pairs=[],fits=[])
    for a,b in zip(good,good[1:]):
        entry['pairs'].append(dict(chis=[a['chi'],b['chi']],c=6*(b['S']-a['S'])/math.log(b['xi']/a['xi']),
            c_renyi={alpha:12*float(alpha)/(float(alpha)+1)*(b['S_renyi'][alpha]-a['S_renyi'][alpha])/math.log(b['xi']/a['xi']) for alpha in ['1','2','3']},
            xi_chi_power=math.log(b['xi']/a['xi'])/math.log(b['chi']/a['chi'])))
    for n in range(3,len(good)+1):
        tail=good[-n:];x=np.log([r['xi'] for r in tail]);y=np.array([r['S'] for r in tail]);s,b=np.polyfit(x,y,1)
        entry['fits'].append(dict(chis=[r['chi'] for r in tail],c=6*s,rms=float(np.sqrt(np.mean((s*x+b-y)**2))),xi_ratio=math.exp(max(x)-min(x)),
            c_renyi={alpha:12*float(alpha)/(float(alpha)+1)*float(np.polyfit(x,[r['S_renyi'][alpha] for r in tail],1)[0]) for alpha in ['1','2','3']}))
    fes.append(entry)
result=dict(inventory=inventory,periodic_sectors=sector,periodic_selected=selected,periodic_pairs=pairs,infinite=fes,
    periodic_energy_gate=energy_keys,periodic_energy_sectors=energy_sectors,periodic_parity_gaps=parity_gaps,
    periodic_selection_exclusions=selection_exclusions,
    periodic_bond_checks=bond_checks,periodic_matched_chi_gaps=matched_gaps,
    interpretation='Conditional entropy estimators. No thermodynamic central charge is assumed; failed convergence rows excluded.')
(HERE/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
for r in selected:print('PBC',r['L'],r['parity'],r['chi'],r['entropy_stats']['1'], 'both parity',r['both_parities_converged'])
for r in pairs:print('SIZE',r['L'],r['c']['1'])
for r in fes:print('FES',r['model'],r['cell'],r['seed'],'gauge',r['gauge'],'one-site',r['one_site'],'points',[(x['chi'],round(x['S'],6),round(x['xi'],3)) for x in r['points']], 'fits',r['fits'])
for r in parity_gaps:print('PARITY GAP',r)
