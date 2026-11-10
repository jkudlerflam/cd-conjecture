#!/usr/bin/env python3
"""Audit completed attractive-route artifacts, without altering checkpoints."""
import hashlib,json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def main():
    rows=[];points=json.loads((HERE/'points.json').read_text())
    expected={(n,L,P) for n,ll in [('attr_ref',[192]),('pref_ref',[48,96,192]),('pref_rat',[48,96]),('pref_zrat',[48,96])] for L in ll for P in (0,1)}
    solver_hash=hashlib.sha256((ROOT/'codex/work_C190/z2dmrg.py').read_bytes()).hexdigest()
    for f in sorted((HERE/'results').glob('*.json')):
        r=json.loads(f.read_text());assert r['status']=='done',(f,'not done')
        assert r['couplings']==points[r['name']],(f,'coupling mismatch')
        assert r['solver_sha256']==solver_hash,(f,'solver mismatch')
        ck=Path(r['checkpoint']);assert ck.exists(),ck
        with np.load(ck,allow_pickle=False) as z:meta=json.loads(str(z['meta']))
        assert meta['status']=='done' and len(meta['history'])==len(r['history']),(f,'checkpoint mismatch')
        assert meta['couplings']==r['couplings'] and meta['parity']==r['parity'] and meta['L']==r['L']
        assert abs(r['parity_expectation']-(-1)**r['parity'])<1e-9,(f,'parity mismatch')
        assert r['convergence']['sweep_converged'],(f,'sweep gate failed')
        assert abs(r['convergence']['canonical_vs_sweep_S'])<1e-6,(f,'canonical entropy mismatch')
        rows.append(dict(file=f.name,sha256=hashlib.sha256(f.read_bytes()).hexdigest(),name=r['name'],L=r['L'],parity=r['parity'],seed=r['seed'],chi=r['chi'],sweeps=len(r['history']),sweep_seconds=sum(h['t'] for h in r['history']),checkpoint=str(ck),checkpoint_sha256=hashlib.sha256(ck.read_bytes()).hexdigest(),checkpoint_bytes=ck.stat().st_size,convergence=r['convergence']))
    assert {(r['name'],r['L'],r['parity']) for r in rows}==expected,'Incomplete or unexpected registered output set'
    assert not list((HERE/'results').glob('pref_rat_L192*')),'Deferred job has raw output'
    out=dict(run_count=len(rows),configuration_count=len(expected)//2,all_sweep_gates_pass=True,all_parity_pairs_complete=True,all_checkpoint_metadata_consistent=True,solver_sha256=solver_hash,summed_sweep_seconds=sum(r['sweep_seconds'] for r in rows),rows=rows,deferred=[dict(name='attr_rat',scope='all registered48/96 calculations',reason='preferred-ratio comparison took priority before any output'),dict(name='pref_rat',scope='192 both parities, chi192',reason='new reference192 drift and order/cat evidence; interrupted at wait gate before any sweep'),dict(name='pref_ref/pref_rat',scope='96 chi128-to192 refinements',reason='same parent-directed deferral'),dict(name='pref_zrat',scope='any192 promotion',reason='not registered; strong z-order/cat warning')],source_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in HERE.glob('*.py')})
    (HERE/'run_inventory.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('rows','source_sha256','deferred')},indent=2))
if __name__=='__main__':main()
