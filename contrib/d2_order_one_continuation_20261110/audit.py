#!/usr/bin/env python3
"""Audit completed numerical records, cached selections and live checkpoint metadata."""
from collections import Counter
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'d2_order_one_20261110'))
import dmrg_screen as D
import numpy as np
sys.path.insert(0,str(HERE))
from analyze import FOLDERS

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--final',action='store_true',help='Also require no active runs or checkpoint races.')
    args=parser.parse_args()
    rows=[];errors=[]
    for folder in FOLDERS:
        for path in sorted((HERE/folder).glob('*.json')):
            raw=path.read_bytes();r=json.loads(raw)
            if 'history' not in r:continue
            row=dict(source=str(path.relative_to(HERE)),status=r.get('status'),checks=[])
            if r.get('status')!='done':
                rows.append(row);continue
            try:
                solver=HERE.parents[1]/r['solver']
                assert hashlib.sha256(solver.read_bytes()).hexdigest()==r['solver_sha256']
                if not r['name'].startswith('control_'):assert r['distance']==1
                row['checks'].append('solver source hash and target interaction range')
                assert len(r['Sprof'])==r['L']-1
                S=np.asarray(r['Sprof']);assert np.all(np.isfinite(S)) and np.min(S)>-1e-10
                canonical=D.entropy_stats(S,r['L'])
                for k in ('S_central3','S_fractional'):assert abs(canonical[k]-r[k])<1e-9,(path,k)
                assert r['history'][-1]['chi']==r['chi']
                assert abs(r['history'][-1]['E_physical']-r['E_physical'])<1e-7
                density=np.asarray(r['measurements']['dens'])
                assert len(density)==r['L'] and np.all(np.isfinite(density))
                assert np.min(density)>-1e-9 and np.max(density)<1+1e-9
                row['checks'].append('shape, finite values, canonical entropy, final energy, density bounds')
                row['sweep_converged']=r['convergence']['sweep_converged']
                if 'snapshots' in folder:
                    row['checkpoint_check']='Historical JSON may name an overwritten live checkpoint; consult preservation manifests.'
                else:
                    cp=Path(r['checkpoint'])
                    with np.load(cp,allow_pickle=False) as z:meta=json.loads(str(z['meta']))
                    if path.read_bytes()!=raw:
                        row['checkpoint_check']='record changed during audit; retry after completion'
                    elif len(meta['history'])!=len(r['history']):
                        row['checkpoint_check']='checkpoint update may be active; retry after completion'
                    else:
                        assert meta['couplings']==r['couplings']
                        assert meta['L']==r['L'] and meta['parity']==r['parity'] and meta['seed']==r['seed']
                        row['checks'].append('checkpoint identity and history length')
            except Exception as exc:
                row['error']=str(exc);errors.append(dict(source=row['source'],error=repr(exc)))
            rows.append(row)
    selections=json.loads((HERE/'size_comparisons.json').read_text())
    for name,case in selections.items():
        for r in case['ground_branches']:
            path=HERE/r['source']
            if hashlib.sha256(path.read_bytes()).hexdigest()!=r['record_sha256']:
                errors.append(dict(source=r['source'],error='selected source hash changed; regenerate analyze.py'))
    result=dict(status_counts=dict(Counter(r['status'] for r in rows)),
                completed_count=sum(r['status']=='done' for r in rows),
                failed_sweep_gates=[r['source'] for r in rows if r.get('sweep_converged') is False],
                pending_checkpoint_checks=[r['source'] for r in rows if 'retry after completion' in r.get('checkpoint_check','')],
                errors=errors,records=rows)
    D.atomic_json(HERE/'audit_report.json',result)
    print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))
    if errors or result['failed_sweep_gates']:raise SystemExit(1)
    if args.final and (result['pending_checkpoint_checks'] or result['status_counts'].get('running',0)):
        raise SystemExit('Final audit requires quiescent snapshots and no active numerical runs.')

if __name__=='__main__':main()
