#!/usr/bin/env python3
"""Reproduce selected archived runs, including parents, in a separate directory."""
import argparse,json,os,subprocess,sys
from pathlib import Path
from common import HERE


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--run-root',type=Path,required=True)
    ap.add_argument('--tag',action='append',help='Repeat for selected runs; default is every completed run.')
    args=ap.parse_args();root=args.run_root.resolve()
    if root==HERE:
        ap.error('Choose a separate directory so the archived records are preserved.')
    root.mkdir(parents=True,exist_ok=True);(root/'logs').mkdir(exist_ok=True)
    tags=args.tag or [Path(x).stem for x in (HERE/'completed_configs.txt').read_text().splitlines()]
    ordered=[];seen=set()
    def visit(tag):
        if tag in seen:return
        if Path(tag).name!=tag:raise ValueError('Use a config tag without path separators.')
        seen.add(tag);cfg=json.loads((HERE/'configs'/(tag+'.json')).read_text())
        if cfg.get('from'):visit(cfg['from'])
        ordered.append(tag)
    for tag in tags:visit(tag)
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    checks=[]
    for tag in ordered:
        print('Reproducing',tag,flush=True)
        with (root/'logs'/(tag+'.txt')).open('a') as log:
            subprocess.run([sys.executable,'-B',str(HERE/'run.py'),str(HERE/'configs'/(tag+'.json')),
                            '--run-root',str(root)],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        actual=json.loads((root/'outputs'/(tag+'.json')).read_text())
        archived=HERE/'outputs'/(tag+'.json')
        ref=json.loads(archived.read_text()) if archived.exists() else {}
        if ref.get('status')=='done':
            checks.append(dict(tag=tag,status=actual['status'],
                delta_E=actual['hist'][-1]['E_phys']-ref['hist'][-1]['E_phys'],
                delta_Sbar={w:actual['Sbar'][w]-ref['Sbar'][w] for w in ('1','6','12')},
                delta_n=actual['hist'][-1]['n_mid']-ref['hist'][-1]['n_mid']))
        else:
            checks.append(dict(tag=tag,status=actual['status'],comparison='No completed archived reference.'))
        (root/'comparison.json').write_text(json.dumps(checks,indent=2)+'\n')
    print('Comparison saved to',root/'comparison.json')


if __name__=='__main__':main()
