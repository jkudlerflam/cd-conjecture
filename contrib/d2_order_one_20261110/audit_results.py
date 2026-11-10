#!/usr/bin/env python3
"""Read-only audit of completed runs and artifacts; no new optimization."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def main():
    rows=[]
    for p in sorted((HERE/'dmrg').glob('*.json')):
        r=json.loads(p.read_text());c=r.get('convergence',{})
        ck=Path(r['checkpoint'])
        rows.append(dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
            name=r['name'],L=r['L'],parity=r['parity'],seed=r['seed'],status=r['status'],
            chi=r.get('chi'),sweeps=len(r['history']),
            elapsed_sweeps_seconds=sum(h['t'] for h in r['history']),
            checkpoint=str(ck),checkpoint_exists=ck.exists(),
            checkpoint_bytes=ck.stat().st_size if ck.exists() else None,
            convergence=c))
    assert all(r['status']=='done' for r in rows),'A run remains active or incomplete'
    assert all(r['checkpoint_exists'] for r in rows),'Missing local checkpoint'
    # Every physical target configuration has both parity sectors available.
    keys={(r['name'],r['L']) for r in rows if not r['name'].startswith('control')}
    assert all({r['parity'] for r in rows if (r['name'],r['L'])==k}=={0,1} for k in keys)
    data=dict(run_count=len(rows),configuration_count=len(keys),
        all_sweep_converged=all(r['convergence'].get('sweep_converged',False) for r in rows),
        failed_sweep_gates=[r['file'] for r in rows if not r['convergence'].get('sweep_converged',False)],
        summed_sweep_seconds=sum(r['elapsed_sweeps_seconds'] for r in rows),
        checkpoint_bytes=sum(r['checkpoint_bytes'] for r in rows),runs=rows,
        source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.glob('*.py')})
    (HERE/'run_inventory.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({k:v for k,v in data.items() if k not in ['runs','source_hashes']},indent=2))


if __name__=='__main__':main()
