#!/usr/bin/env python3
"""Reproducible completed-record inventory and fixed-window size comparisons."""
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'d2_order_one_20261110'
sys.path.insert(0,str(OLD))
import dmrg_screen as D
import ed_screen as E

MAIN='b_r5_2_d1_10_q0_hm7_8'
WANTED={MAIN,'bridge_r4_u1_k4','pref_ref','pref_rat','pref_zrat','attr_ref','simple_attr_7_m1_5_m7'}
ALIASES={'reference_peotta_exploratory':'attr_ref'}
FOLDERS=['validation','bridge_scan/dmrg','bridge_scan/validation_snapshots',
         'attractive_scan/results','attractive_scan/snapshots','simple_attractive/results',
         'independent_seed','independent_sizes']

def main():
    source_dirs=[OLD/'dmrg']+[HERE/f for f in FOLDERS]
    cases={};inventory=[];all_rows=[]
    for directory in source_dirs:
        for path in sorted(directory.glob('*.json')):
            raw=path.read_bytes();r=json.loads(raw)
            if 'history' not in r:continue
            name=ALIASES.get(r['name'],r['name'])
            if directory!=OLD/'dmrg':
                inventory.append(dict(name=name,L=r['L'],parity=r['parity'],seed=r['seed'],
                    status=r.get('status'),source=str(path.relative_to(HERE)),
                    snapshot='snapshots' in str(directory),
                    last_chi=r['history'][-1]['chi'] if r['history'] else None))
            if name not in WANTED or r.get('status')!='done':continue
            c=r['couplings']
            if name in cases:
                assert all(abs(c[k]-cases[name][k])<1e-14 for k in c),(name,c,cases[name])
            cases[name]=c
            row={k:r[k] for k in ('L','parity','seed','chi','E_physical','S_central3','S_fractional','profiles','convergence')}
            row.update(name=name,source=os.path.relpath(path,HERE),record_sha256=hashlib.sha256(raw).hexdigest())
            all_rows.append(row)
    selected={}
    for r in all_rows:
        key=(r['name'],r['L'],r['parity'],r['seed'])
        rank=(r['convergence']['sweep_converged'],r['chi'],-r['E_physical'])
        old=selected.get(key)
        if old is None or rank>(old['convergence']['sweep_converged'],old['chi'],-old['E_physical']):selected[key]=r
    summaries={}
    for name in sorted(cases):
        rr=[r for r in selected.values() if r['name']==name and r['convergence']['sweep_converged']]
        sizes=sorted({r['L'] for r in rr});gs=[]
        for L in sizes:
            atL=[r for r in rr if r['L']==L]
            g=min(atL,key=lambda r:r['E_physical']).copy()
            sectors=sorted({r['parity'] for r in atL})
            g['completed_parities']=sectors;g['both_parities_complete']=sectors==[0,1]
            g['converged_parities']=sectors
            if g['both_parities_complete']:
                ep={P:min(r['E_physical'] for r in atL if r['parity']==P) for P in sectors}
                g['E_odd_minus_even']=ep[1]-ep[0]
            gs.append(g)
        pairs=[]
        for i,a in enumerate(gs):
            for b in gs[i+1:]:
                if b['L']!=2*a['L'] and b!=gs[i+1]:continue
                pairs.append(dict(sizes=[a['L'],b['L']],parities=[a['parity'],b['parity']],
                    both_parities_complete=a['both_parities_complete'] and b['both_parities_complete'],
                    selected_sweeps_converged=True,
                    c_central3=6*(b['S_central3']-a['S_central3'])/math.log(b['L']/a['L']),
                    c_fractional=6*(b['S_fractional']-a['S_fractional'])/math.log(b['L']/a['L'])))
        summaries[name]=dict(couplings=cases[name],metrics=E.coupling_metrics(cases[name]),ground_branches=gs,pairs=pairs)
    D.atomic_json(HERE/'run_inventory.json',dict(generated_at=datetime.now(timezone.utc).isoformat(),
        status_counts=dict(Counter(r['status'] for r in inventory if not r['snapshot'])),records=inventory))
    D.atomic_json(HERE/'size_comparisons.json',summaries)
    lines=['# Completed-record size comparisons','',
        'Generated from exact source paths and hashes in `size_comparisons.json`.',
        'Only completed, sweep-converged records enter this table. Higher-chi checkpoints that are',
        'still running do not displace preserved completed snapshots. Lowest',
        'variational energy selects the parity; this is not a rigorous ground-state',
        'certificate. A missing second parity is explicitly marked provisional.','',
        '| Point | R | Sizes | Central-three slope | Fractional slope | Both parities converged |',
        '|---|---:|---:|---:|---:|---|']
    for name,s in summaries.items():
        for p in s['pairs']:
            a,b=p['sizes'];complete='yes' if p['both_parities_complete'] else 'provisional'
            lines.append(f"| {name} | {s['metrics']['raw_coefficient_spread']:.4g} | {a}→{b} | {p['c_central3']:.6f} | {p['c_fractional']:.6f} | {complete} |")
    lines += ['', 'These are finite-size entropy estimators, not established central charges.',
        'Interpret alongside the order, parity decomposition, Renyi, periodic-ring,',
        'bond-dimension and energy-variance checks. Deferred grids are not results.', '']
    (HERE/'SIZE_COMPARISONS.md').write_text('\n'.join(lines))
    print(json.dumps(dict(status_counts=dict(Counter(r['status'] for r in inventory if not r['snapshot'])),main=summaries.get(MAIN,{})),indent=2))

if __name__=='__main__':main()
