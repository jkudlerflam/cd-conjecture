"""Force a deadline stop at L=8; reload and converge to the literal ED energy."""
import argparse
import json
import os
from pathlib import Path
import tempfile
import time
import numpy as np
from common import write_json
import ims,ofd
from run_stage import LoggedEngine,save,load
from tenpy.algorithms import dmrg
from tenpy.algorithms.exact_diag import ExactDiag


def check():
    P=ims.params({'model':'target','d':'0.1875','h':'1.783'});M=ofd.fmodel(P,8)
    psi=ofd.init_state(M,P,8,0,1)
    opts={'trunc_params':{'chi_max':16,'svd_min':1e-14},'mixer':True,'min_sweeps':4,
          'max_sweeps':30,'N_sweeps_check':1,'max_E_err':1e-12,'max_S_err':1e-10,
          'lanczos_params':{'N_max':30,'E_tol':1e-14,'P_tol':1e-14},'max_trunc_err':1.}
    with tempfile.TemporaryDirectory() as d:
        class StopAfterFirstSweep(LoggedEngine):
            def run_iteration(self):
                result=super().run_iteration()
                self.deadline=time.time()
                return result
        e=StopAfterFirstSweep(psi,M,dict(opts));e.log_path=Path(d)/'sweeps.jsonl';e.deadline=time.time()+3600
        Epartial,_=e.run();assert e.shelve
        save(Path(d)/'state.h5',psi,{'L':8,'parity':0,'h':1.783})
        p2=load(Path(d)/'state.h5')['psi']
        overlap=float(abs(psi.overlap(p2)))
        opts['mixer']=False
        e2=dmrg.TwoSiteDMRGEngine(p2,M,opts);E,_=e2.run()
    ed=ExactDiag(M,charge_sector=[0]);ed.build_full_H_from_mpo();ed.full_diagonalization()
    Eed=float(min(ed.E));error=float(abs(E-Eed))
    assert error<1e-10 and e2.is_converged() and abs(overlap-1)<1e-12
    return {'passed':True,'L':8,'parity':0,'forced_shelve':bool(e.shelve),'partial_sweeps':e.sweeps,
            'partial_E':float(Epartial),'resumed_E':float(E),'sector_ED_E':Eed,'energy_error':error,
            'roundtrip_overlap_abs':overlap,'resumed_converged':bool(e2.is_converged()),
            'semantics':'canonical MPS reload; rebuilt optimizer; no bitwise trajectory claim'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
    r=check();write_json(a.output,r);print(r)
