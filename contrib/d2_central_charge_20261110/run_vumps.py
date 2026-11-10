"""Separate, checkpointed VUMPS stationarity checks at fixed bond dimension."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import argparse,hashlib,json,pickle,time,traceback
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import tenpy
from tenpy.algorithms.vumps import SingleSiteVUMPSEngine,VUMPSEngine
from tenpy.networks.uniform_mps import UniformMPS
from tenpy.networks.mpo import MPOTransferMatrix
import run as R
HERE=Path(__file__).resolve().parent

class Progress(SingleSiteVUMPSEngine):
    def run_iteration(self):
        out=super().run_iteration();s=self.sweep_stats
        print(json.dumps(R.plain({k:s[k][-1] for k in ['sweep','E','Delta_E','S','Delta_S','norm_err','max_split_err','max_N_lanczos']})),flush=True)
        return out

def check(source):
    start=time.time();source=Path(source);deadline=time.monotonic()+7200
    while time.monotonic()<deadline:
        if source.exists():
            r=json.loads(source.read_text())
            if r['status']=='done':break
            if r['status']!='running':raise RuntimeError('Source is not usable: '+r['status'])
        time.sleep(10)
    else:raise RuntimeError('Source did not finish in two hours')
    label=source.stem+'_vumps';target=HERE/'vumps_results'/f'{label}.json'
    assert not target.exists(),target
    ck=Path(r['checkpoint']);sha=hashlib.sha256(ck.read_bytes()).hexdigest();assert sha==r['checkpoint_sha256']
    hashes={}
    for f in [Path(__file__),HERE/'run.py',R.BASE_SOURCE]:
        raw=f.read_bytes();h=hashlib.sha256(raw).hexdigest();hashes[str(f.relative_to(HERE))]=h
        snap=HERE/'source_snapshots'/f'{h}_{f.name}'
        if not snap.exists():snap.write_bytes(raw)
    record=dict(status='running',source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        source_checkpoint_sha256=sha,source_converged=r['converged'],source_hashes=hashes,
        model=r['model'],L=r['L'],chi=r['chi'],gauge=r.get('gauge',False),couplings=r['couplings'],tenpy_version=tenpy.__version__)
    R.atomic(target,record)
    try:
        with ck.open('rb') as f:s=pickle.load(f)
        psi=s['psi'];p=r['couplings'];L=r['L'];gauge=r.get('gauge',False)
        model=(R.GaugeChain if gauge else R.PauliChain)(dict(L=L,bc_MPS='infinite',bc_x='periodic',couplings=p,conserve=None))
        opts=dict(mixer=False,min_sweeps=6,max_sweeps=40,N_sweeps_check=1,max_E_err=1e-11,
            max_S_err=1e-7,max_split_err=1e-9,norm_tol=1e-8,check_overlap=True,
            trunc_params=dict(chi_max=r['chi'],svd_min=1e-13),
            lanczos_params=dict(N_min=12,N_max=500,N_cache=500,reortho=True,P_tol=1e-20,min_gap=1.,E_tol=1e-12))
        record['engine_options']=json.loads(json.dumps(opts))
        u=UniformMPS.from_MPS(psi)
        env,_,_=MPOTransferMatrix.find_init_LP_RP(model.H_MPO,u,calc_E=True,tol=1e-12,maxiter=100000)
        eng=Progress(u,model,opts,resume_data={'init_env_data':env})
        E,psi=eng.run();norm=float(np.linalg.norm(psi.norm_test()));psi.canonical_form()
        u=UniformMPS.from_MPS(psi)
        env,Es,_=MPOTransferMatrix.find_init_LP_RP(model.H_MPO,u,calc_E=True,tol=1e-12,maxiter=100000)
        left,right=VUMPSEngine.tangent_projector_test(SimpleNamespace(psi=u,model=model),env)
        tangent=float(max(left+right));direct=float(R.local_energy(psi,p,gauge).real.mean())
        st=eng.sweep_stats;checks=dict(recent_energy_stable=len(st['E'][-4:])==4 and np.ptp(st['E'][-4:])<1e-8,
            recent_entropy_stable=len(st['S'][-4:])==4 and np.ptp(st['S'][-4:])<1e-5,
            norm_small=norm<1e-7,final_sweep_norm_small=st['norm_err'][-1]<1e-7,
            split_small=st['max_split_err'][-1]<1e-7,tangent_small=tangent<1e-7,
            energy_agreement=max(abs(float(np.real(E))-direct),max(abs(np.asarray(Es)-direct)))<1e-8)
        outck=HERE/'checkpoints_local'/f'{label}.pkl'
        R.save_checkpoint(outck,psi,dict(L=L,geometry='infinite',parity=0,couplings=p,gauge=gauge,checkpoint_kind='final',sweeps=eng.sweeps,algorithm='VUMPS'))
        ent={str(n):psi.entanglement_entropy(n=n).tolist() for n in [1,2,3]}
        record.update(status='done',converged=all(checks.values()),convergence_checks=checks,sweep_stats=st,
            E_physical=direct+model.energy_offset/L,E_direct_shifted=direct,E_returned=float(np.real(E)),
            norm=norm,tangent_left=left,tangent_right=right,max_tangent=tangent,
            environment_energies_shifted=np.real(Es).tolist(),entropy=ent,
            correlation_lengths=np.atleast_1d(psi.correlation_length(target=4,charge_sector=None)).tolist(),
            checkpoint=str(outck),checkpoint_sha256=hashlib.sha256(outck.read_bytes()).hexdigest())
    except Exception as e:
        record.update(status='failed',error=repr(e),traceback=traceback.format_exc())
    record['elapsed_seconds']=time.time()-start
    assert sha==hashlib.sha256(ck.read_bytes()).hexdigest()
    R.atomic(target,record)
    print(json.dumps(R.plain({k:v for k,v in record.items() if k not in ['sweep_stats','traceback','entropy','engine_options']})),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--sources',nargs='+',required=True);a=ap.parse_args()
    for source in a.sources:check(source)
