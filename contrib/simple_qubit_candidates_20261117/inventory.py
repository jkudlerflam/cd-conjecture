#!/usr/bin/env python3
"""Audit completed configs, source hashes, resource use and reproduction order."""
import hashlib,json
from pathlib import Path
from common import HERE,ROOT,Z


def main():
    data={};running=[]
    for f in sorted((HERE/'outputs').glob('*.json')):
        d=json.loads(f.read_text())
        if d['status']=='done':data[d['cfg']['tag']]=(f,d)
        else:running.append(d['cfg']['tag'])
    ordered=[];seen=set()
    def visit(tag):
        if tag in seen:return
        seen.add(tag);d=data[tag][1]
        if d['cfg'].get('from'):visit(d['cfg']['from'])
        ordered.append('configs/'+tag+'.json')
    for tag in data:visit(tag)
    configs={p.stem for p in (HERE/'configs').glob('*.json')}
    rows=[]
    for tag,(f,d) in data.items():
        ck=HERE/'runtime'/(tag+'.npz')
        rows.append(dict(tag=tag,file=str(f.relative_to(HERE)),sha256=hashlib.sha256(f.read_bytes()).hexdigest(),
            L=d['cfg']['L'],P=d['cfg']['P'],chi=d['hist'][-1]['chi'],sweeps=len(d['hist']),
            sweep_seconds=sum(r['t'] for r in d['hist']),runtime_seconds=d['runtime_seconds'],peak_RSS_GB=d['peak_RSS_GB'],
            checkpoint_exists=ck.exists(),checkpoint_bytes=ck.stat().st_size if ck.exists() else None))
    sources=[Path(Z.__file__),ROOT/'codex/work_C190/models.py',ROOT/'codex/work_C190/analyze190.py']+list(HERE.glob('*.py'))
    out=dict(completed=len(data),running=running,unrun_configs=sorted(configs-set(data)-set(running)),
        source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        summed_sweep_seconds=sum(r['sweep_seconds'] for r in rows),
        maximum_peak_RSS_GB=max((r['peak_RSS_GB'] for r in rows),default=0),
        reproduction_order=ordered,runs=rows)
    (HERE/'inventory.json').write_text(json.dumps(out,indent=2)+'\n')
    (HERE/'completed_configs.txt').write_text('\n'.join(ordered)+'\n')
    print('completed',len(data),'running',running,'peakGB',out['maximum_peak_RSS_GB'],'sweepseconds',out['summed_sweep_seconds'])


if __name__=='__main__':main()
