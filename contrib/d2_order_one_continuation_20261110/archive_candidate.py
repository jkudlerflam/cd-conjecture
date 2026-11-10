#!/usr/bin/env python3
"""Copy the winning candidate's completed MPS states to durable local storage.

Binary checkpoints are intentionally Git-ignored. The small manifest is tracked.
No source checkpoint is modified, and every copy is verified by SHA256.
"""
from pathlib import Path
import hashlib
import json
import shutil
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'d2_order_one_20261110'))
import dmrg_screen as D
import numpy as np
NAME='b_r5_2_d1_10_q0_hm7_8'

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def main():
    sources=[]
    for folder in ('bridge_scan/dmrg','independent_seed','independent_sizes'):
        for path in sorted((HERE/folder).glob(NAME+'*.json')):
            r=json.loads(path.read_text())
            assert r['status']=='done' and r['convergence']['sweep_converged'],path
            sources.append((path,Path(r['checkpoint']),r))
    target=next(r for _,_,r in sources if r['L']==192 and r['parity']==0)
    assert target['chi']>=256,'Archive only after the largest-size chi refinement completes.'
    preservation=json.loads((HERE/'bridge_scan/largest_size_preservation.json').read_text())
    path=Path(preservation['raw_record']);cp=Path(preservation['preserved_checkpoint'])
    assert sha(path)==preservation['raw_record_sha256'] and sha(cp)==preservation['checkpoint_sha256']
    sources.append((path,cp,json.loads(path.read_text())))
    out=HERE/'checkpoints_local';out.mkdir(exist_ok=True)
    rows=[]
    for path,cp,r in sources:
        before=(sha(path),sha(cp))
        with np.load(cp,allow_pickle=False) as z:meta=json.loads(str(z['meta']))
        assert meta['couplings']==r['couplings'] and len(meta['history'])==len(r['history'])
        assert all(meta[k]==r[k] for k in ('L','parity','seed'))
        name=f'{NAME}_L{r["L"]}_P{r["parity"]}_s{r["seed"]}_chi{r["chi"]}_mps.npz'
        destination=out/name
        if not destination.exists():shutil.copy2(cp,destination)
        assert sha(destination)==before[1],destination
        assert (sha(path),sha(cp))==before,'Source changed during archival'
        rows.append(dict(L=r['L'],parity=r['parity'],seed=r['seed'],chi=r['chi'],
            source_record=str(path.relative_to(HERE)),record_sha256=before[0],source_checkpoint=str(cp),
            archived_checkpoint=str(destination.relative_to(HERE)),checkpoint_sha256=before[1],bytes=destination.stat().st_size))
    D.atomic_json(HERE/'local_checkpoint_manifest.json',dict(name=NAME,
        note='Durable workspace copies; NPZ files are Git-ignored and are not included in remote publication.',
        total_bytes=sum(r['bytes'] for r in rows),states=rows))
    print(json.dumps(dict(states=len(rows),total_bytes=sum(r['bytes'] for r in rows)),indent=2))

if __name__=='__main__':main()
