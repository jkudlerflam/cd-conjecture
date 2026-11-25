"""One-sweep timing from the identical retained MPS at 1/4/8 threads."""
import argparse
import json
import os
from pathlib import Path
import resource
import time
import numpy as np
import h5py
from threadpoolctl import threadpool_info,threadpool_limits
from common import write_json
import ims,ofd
from tenpy.algorithms import dmrg
from tenpy.tools import hdf5_io


def main(state,output):
    rows=[]
    for repeat in (0,1):
        for threads in ((1,4,8) if repeat==0 else (8,4,1)):
            with h5py.File(state,'r') as f:psi=hdf5_io.load_from_hdf5(f)['psi']
            P=ims.params({'model':'target','d':'0.1875','h':'1.783'});M=ofd.fmodel(P,psi.L)
            chi=int(max(psi.chi))
            with threadpool_limits(limits=threads):
                params={'trunc_params':{'chi_max':chi,'svd_min':1e-12},'mixer':False,
                        'min_sweeps':0,'max_sweeps':0,'N_sweeps_check':1,'max_trunc_err':1.,
                        'lanczos_params':{'N_min':2,'N_max':20,'E_tol':1e-13,'P_tol':1e-13}}
                t=time.time();eng=dmrg.TwoSiteDMRGEngine(psi,M,params);E,_=eng.run();elapsed=time.time()-t
                row={'repeat':repeat,'threads':threads,'L':psi.L,'chi':chi,'sweeps':eng.sweeps,'seconds':elapsed,
                     'E':float(E),'max_rss_mb':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
                     'threadpools':threadpool_info()}
                rows.append(row);write_json(output,{'input_state':state,'rows':rows});print(json.dumps(row),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('state');p.add_argument('--output',required=True);a=p.parse_args()
    main(a.state,a.output)
