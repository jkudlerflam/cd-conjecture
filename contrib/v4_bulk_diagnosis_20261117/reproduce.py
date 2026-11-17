#!/usr/bin/env python3
"""Recompute completed production cases and dependencies in a separate folder."""
import argparse,json,os,subprocess,sys
from pathlib import Path
from phase_common import HERE,write_json


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-root',type=Path,required=True)
    ap.add_argument('--tag',action='append');args=ap.parse_args();root=args.run_root.resolve()
    if root==HERE:ap.error('Use a separate directory to preserve the archived outputs.')
    tags=args.tag or [Path(s).stem for s in (HERE/'completed_configs.txt').read_text().splitlines()]
    order=[];seen=set()
    def visit(tag):
        if tag in seen:return
        assert Path(tag).name==tag;seen.add(tag)
        cfg=json.loads((HERE/'configs'/(tag+'.json')).read_text())
        if cfg.get('parent_output'):visit(Path(cfg['parent_output']).stem)
        order.append(tag)
    for tag in tags:visit(tag)
    (root/'logs').mkdir(parents=True,exist_ok=True)
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    checks=[]
    for tag in order:
        print('reproducing',tag,flush=True)
        with (root/'logs'/(tag+'.txt')).open('a') as log:
            subprocess.run([sys.executable,'-B',str(HERE/'run.py'),str(HERE/'configs'/(tag+'.json')),
                '--run-root',str(root)],stdout=log,stderr=subprocess.STDOUT,env=env,check=True)
        d=json.loads((root/'outputs'/(tag+'.json')).read_text());ref=json.loads((HERE/'outputs'/(tag+'.json')).read_text())
        checks.append(dict(tag=tag,status=d['status'],archived_status=ref['status'],
            delta_E=d['hist'][-1]['E_phys']-ref['hist'][-1]['E_phys'] if ref['hist'] else None,
            comparison_scope='Completed archive' if ref['status']=='done' else 'Interrupted archive; not a completed convergence reference'))
        write_json(root/'reproduction_comparison.json',checks)


if __name__=='__main__':main()
