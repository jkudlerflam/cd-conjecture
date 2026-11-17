#!/usr/bin/env python3
"""Serial exact controls/ED/channel measurements: consumes ONE worker slot."""
import argparse,json,os,subprocess,sys
from threadpoolctl import threadpool_limits
from phase_common import HERE,OLD
from measure_channel import run as measure_channel


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--channels-only',action='store_true')
    ap.add_argument('--benchmark-320',action='store_true');args=ap.parse_args()
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    if not args.channels_only:
        for script in ('validate','controls','ring_ed'):
            print('start',script,flush=True)
            with (HERE/(('validation' if script=='validate' else script)+'.txt')).open('w') as f:
                subprocess.run([sys.executable,'-B',str(HERE/(script+'.py'))],env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
            print('done',script,flush=True)
    with threadpool_limits(limits=1):
        for p in sorted((HERE/'observations').glob('*.json')):
            if (HERE/'channels'/p.name).exists():continue
            d=json.loads(p.read_text())
            cp=(OLD/'runtime'/f"legacy_V4_a1_L96_P{d['P']}.npz") if d['tag'].startswith('reused_') else HERE/'runtime'/f"{d['tag']}_chi{d['chi']}.npz"
            assert cp.exists();measure_channel(cp,p.stem)
    if args.benchmark_320:
        with (HERE/'benchmark_144_320.txt').open('w') as f:
            subprocess.run([sys.executable,'-B',str(HERE/'benchmark.py'),'--L','144','--chi','320','--numpy',
                '--memory-budget-gb','8','--out',str(HERE/'benchmark_144_320.json')],env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
        print('done benchmark320',flush=True)


if __name__=='__main__':main()
