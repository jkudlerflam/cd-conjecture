"""One finite-DMRG stage, based on ofd.py, with retained restart and sweep log.

Restarts reuse the variational state; the new engine rebuilds environments
and restarts its sweep counter/mixer. This is not a bitwise engine resume.
No prior stage is deleted. The deadline returns a canonical checkpoint
with status=partial rather than falsely marking convergence.
"""
import argparse
import gzip
import json
import os
from pathlib import Path
import resource
import sys
import time
import h5py
import numpy as np
from common import source_hashes, write_json
import ims
import ofd
import meas_quartet
import measure_matrix
from tenpy.algorithms import dmrg
from tenpy.tools import hdf5_io
from tenpy.tools.misc import setup_logging


def rss_mb():
    x=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return x/(1024 if sys.platform!='darwin' else 1024**2)


def save(path,psi,meta):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.tmp')
    with h5py.File(temp,'w') as f:
        hdf5_io.save_to_hdf5(f,{'psi':psi,'metadata_json':json.dumps(meta)})
        f.flush()
    temp.replace(path)


def load(path):
    with h5py.File(path,'r') as f: return hdf5_io.load_from_hdf5(f)


class LoggedEngine(dmrg.TwoSiteDMRGEngine):
    def run_iteration(self):
        t=time.time(); result=super().run_iteration()
        st=self.sweep_stats
        row={k:float(st[k][-1]) for k in ('E','Delta_E','S','Delta_S','max_trunc_err','max_E_trunc','norm_err')}
        row.update(sweep=int(self.sweeps),seconds=time.time()-t,max_rss_mb=rss_mb(),chi=int(max(self.psi.chi)))
        # First Delta_E can be NaN. JSON null preserves that it was not measured.
        row={k:(v if not isinstance(v,float) or np.isfinite(v) else None) for k,v in row.items()}
        with self.log_path.open('a') as f: f.write(json.dumps(row)+'\n')
        print('SWEEP '+json.dumps(row),flush=True)
        self.last_sweep_seconds=row['seconds']
        return result

    def stopping_criterion(self,iteration_start_time):
        # Budget a further slow sweep plus cleanup/checkpoint, before starting it.
        if time.time()+max(2*getattr(self,'last_sweep_seconds',0),90) >= self.deadline:
            self.shelve=True
            return True
        return super().stopping_criterion(iteration_start_time)


def run(c):
    t0=time.time(); out=Path(c['results']); out.mkdir(parents=True,exist_ok=True)
    stem=f"L{c['L']}_p{c['parity']}_chi{c['chi']}_{c.get('tag','target')}"
    P=ims.params({'model':'target','d':'0.1875','h':str(c.get('h',1.783))})
    M=ofd.fmodel(P,c['L'])
    if c.get('input_state'):
        old=load(c['input_state']); psi=old['psi']
        assert psi.L==c['L']
        if 'metadata_json' in old:
            meta=json.loads(old['metadata_json'])
            for k in ('L','parity','h'):
                assert meta[k]==(c.get('h',1.783) if k=='h' else c[k]), (k,meta[k],c.get(k))
    else:
        psi=ofd.init_state(M,P,c['L'],c['parity'],c.get('seed',1))
    parity=float(np.real(psi.expectation_value_term([('Sigmaz',i) for i in range(c['L'])])))
    assert abs(parity-(-1)**c['parity'])<1e-9,parity
    first=not c.get('input_state')
    deadline=float(os.environ.get('JOB_END_EPOCH',time.time()+3600))-c.get('reserve_seconds',300)
    assert time.time()+120<deadline,'Insufficient time before checkpoint deadline'
    tolE=c.get('tolE',1e-9); tolS=c.get('tolS',1e-5)
    pars={'trunc_params':{'chi_max':c['chi'],'svd_min':1e-12},
          'mixer':True,'mixer_params':{'amplitude':1e-4 if first else 1e-6,'decay':2.,'disable_after':12 if first else 4},
          'min_sweeps':14 if first else 4,'max_sweeps':c.get('max_sweeps',40),'N_sweeps_check':1,
          'max_E_err':tolE,'max_S_err':tolS,'max_trunc_err':1.,
          'max_hours':max(0.,(deadline-time.time())/3600),
          'lanczos_params':{'N_min':2,'N_max':20,'E_tol':1e-13,'P_tol':1e-13}}
    eng=LoggedEngine(psi,M,pars); eng.log_path=out/(stem+'_sweeps.jsonl'); eng.deadline=deadline
    E,_=eng.run()
    S=psi.entanglement_entropy(); Z=psi.expectation_value('Sigmaz').real
    Epost=float(np.real(M.H_MPO.expectation_value(psi)))
    converged=bool(eng.is_converged() and eng.mixer is None)
    st=eng.sweep_stats
    meta={'L':c['L'],'parity':c['parity'],'h':c.get('h',1.783),'delta':3/16,'chi':c['chi'],
          'chi_actual':int(max(psi.chi)),'E':Epost,'E_sweep':float(E),'parity_expectation':parity,
          'Smid':float(S[c['L']//2-1]),'Sbar':float(np.mean(S[c['L']//4-1:3*c['L']//4])),
          'converged':converged,'shelved':bool(eng.shelve),'sweeps':int(eng.sweeps),
          'dE_last':float(st['Delta_E'][-1]),'dS_last':float(st['Delta_S'][-1]),
          'truncation':float(st['max_trunc_err'][-1]),'dmrg_seconds':time.time()-t0,
          'input_state':c.get('input_state'),'source_sha256':source_hashes(),'config':c,
          'restart_semantics':'MPS warm restart with rebuilt environments; counters and mixer reset'}
    # Save first, so measurement interruption cannot lose optimization work.
    save(c['output_state'],psi,meta)
    loaded=load(c['output_state'])['psi']
    meta['roundtrip_energy_abs']=float(abs(M.H_MPO.expectation_value(loaded)-Epost))
    meta['roundtrip_Z_max_abs']=float(np.max(abs(loaded.expectation_value('Sigmaz').real-Z)))
    meta['roundtrip_S_max_abs']=float(np.max(abs(loaded.entanglement_entropy()-S)))
    meta['roundtrip_norm_overlap_abs']=float(abs(abs(psi.overlap(loaded))-1))
    assert max(meta[k] for k in meta if k.startswith('roundtrip_'))<1e-8
    del loaded
    np.savez_compressed(out/(stem+'_profiles.npz'),S=S,Z=Z)
    write_json(out/(stem+'_summary.json'),meta)
    print('STAGE '+json.dumps({k:v for k,v in meta.items() if k not in ('source_sha256','config')}),flush=True)
    if c.get('measure') and time.time()+c.get('measurement_min_seconds',180)<float(os.environ['JOB_END_EPOCH'])-90:
        refs=c.get('refs',[c['L']//2-17,c['L']//2-8,c['L']//2])
        m=meas_quartet.measure(psi,c['L'],refs,meas_quartet.DS_DEFAULT)
        m.update(chi=int(max(psi.chi)),parity=c['parity'],h=c.get('h',1.783),E=Epost)
        with gzip.open(out/(stem+'_quartet.json.gz'),'wt') as f: json.dump(m,f)
        if c.get('matrix'):
            m=measure_matrix.measure(psi,c.get('matrix_refs',[c['L']//2-17]),rmax=c.get('matrix_rmax',96))
            m.update(parity=c['parity'],h=c.get('h',1.783),E=Epost)
            with gzip.open(out/(stem+'_matrix.json.gz'),'wt') as f: json.dump(m,f)
        meta['measured']=True
    else: meta['measured']=False
    meta.update(total_seconds=time.time()-t0,max_rss_mb=rss_mb(),checkpoint_bytes=Path(c['output_state']).stat().st_size)
    write_json(out/(stem+'_summary.json'),meta)
    print('FINISHED '+json.dumps({k:meta[k] for k in ('L','chi','E','converged','measured','total_seconds','max_rss_mb')}),flush=True)
    return meta


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');a=p.parse_args()
    setup_logging(to_stdout='WARNING');run(json.loads(Path(a.config).read_text()))
