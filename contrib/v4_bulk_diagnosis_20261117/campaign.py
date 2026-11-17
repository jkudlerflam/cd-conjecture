#!/usr/bin/env python3
"""Run an explicit dependency-ready list with TWO single-thread workers."""
import argparse,concurrent.futures,subprocess,os,sys
from pathlib import Path
from phase_common import HERE


def main():
    ap=argparse.ArgumentParser();ap.add_argument('configs',nargs='+');args=ap.parse_args()
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    def one(p):
        p=Path(p);print('start',p.stem,flush=True)
        with (HERE/'logs'/(p.stem+'.txt')).open('a') as f:
            subprocess.run([sys.executable,'-B',str(HERE/'run.py'),str(p)],env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
        print('done',p.stem,flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(one,args.configs))


if __name__=='__main__':main()
