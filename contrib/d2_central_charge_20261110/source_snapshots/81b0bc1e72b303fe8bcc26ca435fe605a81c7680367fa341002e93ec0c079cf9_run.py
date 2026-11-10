#!/usr/bin/env python3
"""Independent full-spin periodic and infinite DMRG, fixed Pauli candidate."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[k]='1'
import argparse, hashlib, importlib.util, json, logging, pickle, sys, time
from pathlib import Path
import numpy as np
import tenpy
from tenpy.networks.mps import MPS

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BASE_SOURCE=HERE/'source_snapshots/162f416856b4a578abf3f2966b2d601c4d1538a6aee17fb153d3a75dd2927350_run.py'
spec=importlib.util.spec_from_file_location('_cc_frozen_tenpy',BASE_SOURCE)
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
PauliChain,ProgressEngine,plain,save_checkpoint=base.PauliChain,base.ProgressEngine,base.plain,base.save_checkpoint

def local_energy(psi,p):
    terms=[]
    pairs=[('Sigmax','Sigmax',p['Jx']),('Sigmay','Sigmay',p['Jy']),('Sigmaz','Sigmaz',p['Jz']),
        ('Sigmax','Sigmay',p['Dz']),('Sigmay','Sigmax',-p['Dz'])]
    for i in range(psi.L):
        e=p['hz']*psi.expectation_value_term([('Sigmaz',i)])
        for a,b,j in pairs:
            if j:e+=j*psi.expectation_value_term([(a,i),(b,i+1)])
        terms.append(e-p['Jz']-p['hz'])
    return np.asarray(terms)

PARAMS=dict(Jx=.5,Jy=-.5,Jz=.625,Dz=.1,Dy=0.,hz=-.875)
XX=dict(Jx=.5,Jy=.5,Jz=0.,Dz=0.,Dy=0.,hz=0.)

def atomic(path,obj):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(plain(obj),indent=2,allow_nan=False)+'\n');tmp.replace(path)

def stats(S,L,factor):
    S=np.asarray(S); x=np.arange(1,L)
    lo,hi=3*L/8,5*L/8; xx=np.r_[lo,x[(x>lo)&(x<hi)],hi]
    out=dict(central3=float(S[(x>=L//2-1)&(x<=L//2+1)].mean()),
        fractional=float(np.trapezoid(np.interp(xx,x,S),xx)/(hi-lo)))
    out['profiles']={}
    for f in (.25,.375):
        mask=(x>=f*L)&(x<=(1-f)*L)
        chord=np.log((L/np.pi)*np.sin(np.pi*x[mask]/L))
        a,b=np.polyfit(chord,S[mask],1)
        out['profiles'][str(f)]=dict(c=factor*float(a),rms=float(np.sqrt(np.mean((a*chord+b-S[mask])**2))))
    return out

def run(a):
    start=time.time(); p=PARAMS.copy() if a.model=='candidate' else XX.copy()
    source_dir=HERE/'source_snapshots';source_dir.mkdir(exist_ok=True)
    for source in [Path(__file__),BASE_SOURCE]:
        raw=source.read_bytes();sha=hashlib.sha256(raw).hexdigest()
        target=source_dir/f'{sha}_{source.name}'
        if not target.exists():target.write_bytes(raw)
    infinite=a.geometry=='infinite'; L=a.L
    label=f'{a.model}_{a.geometry}_L{L}_P{a.parity}_s{a.seed}_chi{a.chi}'
    if a.attempt:label+=f'_a{a.attempt}'
    out=HERE/'results'/f'{label}.json'; ck=HERE/'checkpoints_local'/f'{label}.pkl'
    if out.exists():
        old=json.loads(out.read_text())
        if old.get('status')=='done': print('EXISTS',out);return old
        raise RuntimeError(f'Partial record exists: {out}; use a new seed or explicit recovery')
    model=PauliChain(dict(L=L,bc_MPS='infinite' if infinite else 'finite',
        bc_x='periodic',couplings=p,conserve=None if infinite else 'parity'))
    rng=np.random.default_rng(a.seed)
    if a.load:
        with open(a.load,'rb') as f: saved=pickle.load(f)
        assert saved['couplings']==p and saved['L']==L and saved['geometry']==a.geometry and saved['parity']==a.parity
        psi=saved['psi'];psi.canonical_form()
    else:
        density=.266 if a.model=='candidate' else .5
        if infinite:
            phases=rng.uniform(0,2*np.pi,L) if a.seed!=17 else np.full(L,-np.pi/2)
            state=[np.array([np.sqrt(1-density),np.sqrt(density)*np.exp(1j*q)]) for q in phases]
        else:
            n=round(L*density);n+=(a.parity-n)%2
            occupied=set(rng.choice(L,n,replace=False))
            state=['down' if i in occupied else 'up' for i in range(L)]
        psi=MPS.from_product_state(model.lat.mps_sites(),state,bc=model.lat.bc_MPS,dtype=complex,unit_cell_width=L)
    meta=dict(L=L,geometry=a.geometry,parity=a.parity,couplings=p)
    record=dict(**meta,model=a.model,seed=a.seed,chi=a.chi,status='running',checkpoint=str(ck),
        config=vars(a),tenpy_version=tenpy.__version__,source_hashes={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()
            for f in [Path(__file__),BASE_SOURCE]})
    atomic(out,record)
    opts=dict(mixer=not bool(a.load),max_sweeps=a.sweeps,min_sweeps=min(12,a.sweeps),
        N_sweeps_check=2 if infinite else 1,max_E_err=1e-11,max_S_err=1e-7,
        P_tol_to_trunc=None,diag_method='default',
        trunc_params=dict(chi_max=a.chi,svd_min=1e-13),
        lanczos_params=dict(N_min=12,N_max=250,N_cache=250,reortho=True,P_tol=1e-20,min_gap=1.,E_tol=1e-12),
        mixer_params=dict(amplitude=1e-5,decay=2.,disable_after=8))
    if not a.load: opts['chi_list']={0:min(16,a.chi),2:min(32,a.chi),4:a.chi}
    if infinite: opts['update_env']=20
    record['engine_options']=json.loads(json.dumps(opts))
    kwargs=dict(resume_data={'init_env_data':{'start_env_sites':L}}) if infinite else {}
    eng=ProgressEngine(psi,model,opts,checkpoint_path=str(ck),checkpoint_every=2,
        checkpoint_metadata=meta,local_residual_tol=1e-8,**kwargs)
    E,psi=eng.run()
    norm_before=float(np.linalg.norm(psi.norm_test()))
    psi.canonical_form()
    save_checkpoint(ck,psi,dict(meta,checkpoint_kind='final',sweeps=eng.sweeps))
    st=eng.sweep_stats
    n=(1-psi.expectation_value('Sigmaz').real)/2
    ent={str(alpha):psi.entanglement_entropy(n=alpha).tolist() for alpha in (1,2,3)}
    direct=float(local_energy(psi,p).real.mean()) if infinite else float(model.H_MPO.expectation_value(psi))
    local=eng.local_iteration_stats[-1]
    recentE=np.asarray(st['E'][-4:]);recentS=np.asarray(st['S'][-4:])
    checks=dict(recent_energy_stable=len(recentE)==4 and np.ptp(recentE)<(1e-8 if infinite else 1e-7),
        recent_entropy_stable=len(recentS)==4 and np.ptp(recentS)<1e-5,
        local_residual_small=local['max_absolute_residual']<1e-8,
        norm_before_small=norm_before<1e-5,
        final_sweep_norm_small=float(st['norm_err'][-1])<1e-5,
        energy_agreement=abs(E-direct)<(1e-7 if infinite else 1e-8))
    record.update(status='done',convergence_checks=checks,converged=all(checks.values()),
        E_returned=float(E),E_direct_shifted=direct,E_physical=direct+model.energy_offset/(L if infinite else 1),
        norm_before=norm_before,norm_after=float(np.linalg.norm(psi.norm_test())),
        entropy=ent,density=n.tolist(),density_mean=float(n.mean()),density_range=float(np.ptp(n)),
        X=psi.expectation_value('Sigmax').real,Y=psi.expectation_value('Sigmay').real,
        sweep_stats=st,local_solver=eng.local_iteration_stats,elapsed_seconds=time.time()-start,
        checkpoint_sha256=hashlib.sha256(ck.read_bytes()).hexdigest())
    if infinite:
        record['correlation_lengths']=np.atleast_1d(psi.correlation_length(target=4,charge_sector=None)).tolist()
    else:
        record['entropy_stats']={key:stats(value,L,3*2*float(key)/(float(key)+1)) for key,value in ent.items()}
    atomic(out,record)
    print(json.dumps(plain({k:v for k,v in record.items() if k not in ('sweep_stats','entropy','local_solver')})),flush=True)
    return record

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--geometry',choices=['periodic','infinite'],required=True)
    ap.add_argument('--model',choices=['candidate','xx'],default='candidate')
    ap.add_argument('--L',type=int,required=True);ap.add_argument('--parity',type=int,choices=[0,1],default=0)
    ap.add_argument('--seed',type=int,default=17);ap.add_argument('--chi',type=int,default=64)
    ap.add_argument('--sweeps',type=int,default=40);ap.add_argument('--load')
    ap.add_argument('--attempt',type=int,default=0)
    a=ap.parse_args();logging.basicConfig(level=logging.WARNING);run(a)

if __name__=='__main__':main()
