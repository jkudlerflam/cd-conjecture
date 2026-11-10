#!/usr/bin/env python3
"""Read-only parity selection and size/chi diagnostics of finished records."""
import json,math,hashlib
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent

def small(r):
    c=r['measurements']['long_connected_correlations']
    tails={a:dict(last=c[a][-1],rms_last_quarter=float(np.sqrt(np.mean([v*v for d,v in c[a] if d>=3*r['L']/8]))),max_abs_last_quarter=float(max(abs(v) for d,v in c[a] if d>=3*r['L']/8))) for a in ('X','Y','Z')}
    hist=r['history']; stages=[]
    for chi in sorted(set(h['chi'] for h in hist)):
        q=[h for h in hist if h['chi']==chi and h['alpha']==0][-1]
        stages.append({k:q.get(k) for k in ('chi','E_physical','S_central3','S_fractional','res_max','trunc')})
    return {k:r[k] for k in ('name','L','parity','chi','E_physical','S_central3','S_fractional','profiles','density_bulk','density_bulk_range','parity_expectation','convergence','source_record','source_record_sha256')}|dict(correlation_tail=tails,stage_terminal=stages)

def main():
    raw=[]
    for p in sorted((HERE/'results').glob('*.json')):
        r=json.loads(p.read_text())
        if r['status']=='done':
            r['source_record']=str(p);r['source_record_sha256']=hashlib.sha256(p.read_bytes()).hexdigest();raw.append(r)
    new_completed_count=len(raw)
    # Reference smaller sizes remain read-only in their original contribution.
    for L in (48,96):
        for P in (0,1):
            old=HERE.parents[1]/'d2_order_one_20261110'/'dmrg'/f'reference_peotta_exploratory_L{L}_P{P}_s17.json'
            r=json.loads(old.read_text());r['name']='attr_ref';r['source_record']=str(old);r['source_record_sha256']=hashlib.sha256(old.read_bytes()).hexdigest();raw.append(r)
    rows=[small(r) for r in raw]; selected={}; pairs=[]
    for n in sorted(set(r['name'] for r in raw)):
        selected[n]={}
        for L in sorted(set(r['L'] for r in raw if r['name']==n)):
            rr=[r for r in raw if r['name']==n and r['L']==L]
            if set(r['parity'] for r in rr)!={0,1}:continue
            rr=sorted(rr,key=lambda r:r['E_physical'])
            win=small(rr[0]);win['parity_gap']=rr[1]['E_physical']-rr[0]['E_physical'];win['other_parity_delta_S_fractional']=rr[1]['S_fractional']-rr[0]['S_fractional'];selected[n][L]=win
        ll=sorted(selected[n])
        for a,b in zip(ll,ll[1:]):
            aa,bb=selected[n][a],selected[n][b]
            pairs.append(dict(name=n,sizes=[a,b],selected_parities=[aa['parity'],bb['parity']],c_central3=6*(bb['S_central3']-aa['S_central3'])/math.log(b/a),c_fractional=6*(bb['S_fractional']-aa['S_fractional'])/math.log(b/a),both_sweep_converged=aa['convergence']['sweep_converged'] and bb['convergence']['sweep_converged']))
    data=dict(completed_new=new_completed_count,completed_baseline=len(rows)-new_completed_count,rows=rows,selected=selected,pairs=pairs)
    (HERE/'summary.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(completed_new=new_completed_count,completed_baseline=len(rows)-new_completed_count,pairs=pairs,selected={n:{L:{k:r[k] for k in ('parity','E_physical','S_fractional','parity_gap','other_parity_delta_S_fractional','profiles')} for L,r in dd.items()} for n,dd in selected.items()}),indent=2))
if __name__=='__main__':main()
