#!/usr/bin/env python3
"""Run dependency-ready configs with at most two single-thread processes.

This does not schedule dependencies; use reproduce.py for a fresh full rerun.
"""
import argparse,concurrent.futures,os,subprocess,sys
from pathlib import Path


def main():
    p=argparse.ArgumentParser();p.add_argument('configs',nargs='+');args=p.parse_args()
    here=Path(__file__).resolve().parent;logs=here/'logs';logs.mkdir(exist_ok=True)
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    def one(cp):
        path=Path(cp);print('start',path.stem,flush=True)
        with (logs/(path.stem+'.txt')).open('a') as f:
            subprocess.run([sys.executable,'-B',str(here/'run.py'),str(path)],stdout=f,stderr=subprocess.STDOUT,env=env,check=True)
        print('done',path.stem,flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(one,args.configs))


if __name__=='__main__':main()
