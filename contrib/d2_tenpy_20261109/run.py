"""Independent full-spin TeNPy DMRG; no number projection or effective truncation.

Optional exact parity sectors are available only when transverse DM is absent.

Pauli convention, uniform fields, open boundaries. A scalar is subtracted so
the all-up product has energy zero; E_unshifted restores the literal chain.
"""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(key, '1')
import argparse
import json
import logging
import pickle
import time
from pathlib import Path
import numpy as np
import tenpy
from tenpy.models.model import CouplingMPOModel
from tenpy.networks.site import SpinHalfSite
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg
from tenpy.linalg import np_conserved as npc


def couplings(name, hz=None):
    if name == 'c183':
        p = dict(Jx=.4979983, Jy=-.5020017, Jz=4., Dz=.0152122,
                 Dy=.13693064, hz=-7.9463)
    elif name == 'dyadic':
        p = dict(Jx=255/512, Jy=-257/512, Jz=4., Dz=1/64,
                 Dy=9/64, hz=-8137/1024)
    elif name == 'original':
        p = dict(Jx=.49966754391072715, Jy=-.5003324560892728,
                 Jz=10., Dz=.00625, Dy=.075, hz=-19.9774375)
    elif name == 'xx':
        p = dict(Jx=.5, Jy=.5, Jz=0., Dz=0., Dy=0., hz=0.)
    elif name == 'c181':
        p = dict(Jx=.500098723, Jy=-.499901277, Jz=10., Dz=.006249220,
                 Dy=0., hz=-19.97589125)
    elif name == 'c181153':
        w=np.exp(1.53j)/80
        p=dict(Jx=(1+w.real)/2,Jy=(-1+w.real)/2,Jz=10.,Dz=w.imag/2,Dy=0.,hz=-19.97589125)
    elif name == 'two_xx':
        # Two independent even/odd sublattices. This is an estimator control,
        # NOT a nearest-neighbour d=2 counterexample.
        p = dict(Jx=.5, Jy=.5, Jz=0., Dz=0., Dy=0., hz=0., distance=2)
    else:
        raise ValueError(name)
    if hz is not None:
        p['hz'] = hz
    return p


class PauliChain(CouplingMPOModel):
    def init_sites(self, model_params):
        return SpinHalfSite(conserve=model_params.get('conserve',None))

    def init_terms(self, model_params):
        p = model_params['couplings']
        L = self.lat.N_sites
        distance = p.get('distance', 1)
        pairs = [('Sigmax','Sigmax',p['Jx']), ('Sigmay','Sigmay',p['Jy']),
                 ('Sigmaz','Sigmaz',p['Jz']), ('Sigmax','Sigmay',p['Dz']),
                 ('Sigmay','Sigmax',-p['Dz']), ('Sigmax','Sigmaz',p['Dy']),
                 ('Sigmaz','Sigmax',-p['Dy'])]
        for a,b,j in pairs:
            if j:
                self.add_coupling(j,0,a,0,b,distance)
        self.add_onsite(p['hz'],0,'Sigmaz')
        bonds = L if self.lat.bc_MPS == 'infinite' else L-distance
        self.energy_offset = bonds*p['Jz']+L*p['hz']
        self.add_onsite(-self.energy_offset/L,0,'Id')


def save_checkpoint(path, psi, metadata):
    """Atomically save a trusted local state without altering a running engine."""
    target=Path(path);target.parent.mkdir(parents=True,exist_ok=True)
    temporary=target.with_name(target.name+'.tmp')
    with temporary.open('wb') as handle:
        pickle.dump(dict(metadata,psi=psi),handle,pickle.HIGHEST_PROTOCOL)
    os.replace(temporary,target)


def local_summary(rows, tolerance):
    residuals=[r['absolute_residual'] for r in rows]
    cap_count=sum(r['iteration_cap_reached'] for r in rows)
    finite=all(np.isfinite(r) for r in residuals)
    maximum=max(residuals,default=None)
    status=('not_measured' if not rows else 'nonfinite_residual' if not finite else
            'iteration_cap_reached' if cap_count else
            'residual_threshold_unmet' if maximum>tolerance else
            'residual_threshold_pass_no_global_certificate')
    return dict(status=status,updates=len(rows),max_absolute_residual=maximum,
                residual_threshold=tolerance,residual_threshold_pass=bool(rows) and finite and maximum<=tolerance,
                iteration_cap_count=cap_count,max_N_lanczos=max((r['N_lanczos'] for r in rows),default=None))


class ProgressEngine(dmrg.TwoSiteDMRGEngine):
    def __init__(self,*args,local_residual_tol=1e-8,checkpoint_path=None,
                 checkpoint_every=0,checkpoint_metadata=None,**kwargs):
        self.local_residual_tol=local_residual_tol
        self.local_solver_checks=[];self.local_iteration_stats=[]
        self.checkpoint_path=checkpoint_path;self.checkpoint_every=checkpoint_every
        self.checkpoint_metadata=checkpoint_metadata or {};self.checkpoint_sweeps=[]
        super().__init__(*args,**kwargs)

    def diag(self,theta_guess):
        E,theta,N,overlap_change=super().diag(theta_guess)
        # One additional Hamiltonian application; no additional eigensolve.
        # This checks the untruncated local solution, not the global MPS residual.
        residual=float(npc.norm(self.eff_H.matvec(theta)-E*theta)/npc.norm(theta))
        limit=int(self.lanczos_params.get('N_max',20))
        cap=bool(N>0 and N>=limit)
        self.local_solver_checks.append(dict(sweep_before_update=int(self.sweeps),bond=int(self.i0),
            N_lanczos=int(N),N_max=limit,iteration_cap_reached=cap,
            absolute_residual=residual,residual_threshold_pass=bool(np.isfinite(residual) and residual<=self.local_residual_tol)))
        return E,theta,N,overlap_change

    def run_iteration(self):
        first=len(self.local_solver_checks)
        result=super().run_iteration()
        st=self.sweep_stats
        summary=local_summary(self.local_solver_checks[first:],self.local_residual_tol)
        summary['sweep']=int(st['sweep'][-1]);self.local_iteration_stats.append(summary)
        if self.checkpoint_every and summary['sweep']%self.checkpoint_every==0:
            save_checkpoint(self.checkpoint_path,self.psi,dict(self.checkpoint_metadata,
                checkpoint_kind='sweep_boundary',checkpoint_sweep=summary['sweep']))
            self.checkpoint_sweeps.append(summary['sweep'])
        progress={k:plain(st[k][-1]) for k in ('sweep','E','Delta_E','Delta_S','max_trunc_err','max_chi')}
        progress['local_solver']=summary
        print(json.dumps(plain(progress),allow_nan=False),flush=True)
        return result


def entropy_fits(S):
    S=np.asarray(S); L=len(S)+1; cuts=np.arange(1,L)
    x=np.log(2*L/np.pi*np.sin(np.pi*cuts/L))
    out={}
    for frac in (.125,.20,.25,.33):
        mask=(cuts>=frac*L)&(cuts<=(1-frac)*L)
        A=np.column_stack([x[mask]/6,np.ones(mask.sum())])
        coeff=np.linalg.lstsq(A,S[mask],rcond=None)[0]
        out[str(frac)]=dict(c=float(coeff[0]),rms=float(np.sqrt(np.mean((A@coeff-S[mask])**2))))
    return out


def plain(x):
    if isinstance(x,dict): return {str(k):plain(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [plain(v) for v in x]
    if isinstance(x,np.ndarray): return plain(x.tolist())
    if isinstance(x,np.generic): return plain(x.item())
    if isinstance(x,complex): return [x.real,x.imag]
    if isinstance(x,float) and not np.isfinite(x): return None
    if hasattr(x,'eps') and hasattr(x,'ov'):
        return dict(eps=plain(x.eps),ov=plain(x.ov))
    return x


def run(cfg):
    start=time.time(); L=cfg['L']; p=couplings(cfg['name'],cfg.get('hz'))
    checkpoint_every=cfg.get('checkpoint_every',0)
    if checkpoint_every<0 or (checkpoint_every and not cfg.get('save_state')):
        raise ValueError('checkpoint_every must be nonnegative and needs save_state when enabled.')
    residual_tol=cfg.get('local_residual_tol',1e-8)
    if not np.isfinite(residual_tol) or residual_tol<=0:
        raise ValueError('local_residual_tol must be finite and positive.')
    parity=cfg.get('parity')
    if parity and p['Dy']:
        raise ValueError('Transverse DM breaks spin parity; projection forbidden.')
    model=PauliChain(dict(L=L,bc_MPS='finite',bc_x='open',couplings=p,conserve='parity' if parity else None))
    rng=np.random.default_rng(cfg.get('seed',1))
    density=cfg.get('init_density',.12)
    # Random complex product state provides a nonzero overlap with all sectors.
    state=[np.array([np.sqrt(1-density),np.sqrt(density)*np.exp(2j*np.pi*rng.random())]) for _ in range(L)]
    if parity:
        count=max(2,round(density*L))
        desired=0 if parity=='even' else 1
        count += (desired-count)%2
        occupied=set(np.linspace(1,L-2,count,dtype=int).tolist())
        state=['down' if j in occupied else 'up' for j in range(L)]
    psi=MPS.from_product_state(model.lat.mps_sites(),state,bc='finite',dtype=complex)
    if cfg.get('resume_state'):
        # Only load checkpoints created by this script and trusted by the caller.
        with open(cfg['resume_state'],'rb') as handle:
            saved=pickle.load(handle)
        if saved['L']!=L or saved['couplings']!=p or saved['parity']!=parity:
            raise ValueError('Checkpoint Hamiltonian, length or parity does not match.')
        psi=saved['psi']
        # Periodic snapshots may contain the running engine's noncanonical gauge.
        psi.canonical_form()
    chi=cfg.get('chi',96); sweeps=cfg.get('sweeps',16)
    options=dict(mixer=True,max_sweeps=sweeps,min_sweeps=min(8,sweeps),
                 N_sweeps_check=1,max_E_err=1e-11,max_S_err=1e-7,
                 trunc_params=dict(chi_max=chi,svd_min=1e-12,trunc_cut=1e-12),
                 chi_list={0:min(16,chi),2:min(32,chi),4:min(64,chi),6:chi},
                 lanczos_params=dict(N_max=80,P_tol=1e-13),
                 mixer_params=dict(amplitude=1e-5,decay=2.,disable_after=8))
    if cfg.get('resume_state'):
        options['chi_list']={0:min(max(psi.chi),chi),2:chi}
        options['mixer']=False
    if cfg.get('strict_local'):
        options['P_tol_to_trunc']=None
        options['lanczos_params'].update(N_max=160,N_cache=160,reortho=True)
    checkpoint_metadata=dict(L=L,couplings=p,parity=parity)
    eng=ProgressEngine(psi,model,options,local_residual_tol=residual_tol,
        checkpoint_path=cfg.get('save_state'),checkpoint_every=checkpoint_every,
        checkpoint_metadata=checkpoint_metadata)
    E,psi=eng.run()
    psi.canonical_form()
    if cfg.get('save_state'):
        save_checkpoint(cfg['save_state'],psi,dict(checkpoint_metadata,checkpoint_kind='final',
            checkpoint_sweep=int(eng.sweep_stats['sweep'][-1])))
    S=psi.entanglement_entropy()
    n=(1-np.real(psi.expectation_value('Sigmaz')))/2
    X=np.real(psi.expectation_value('Sigmax')); Y=np.real(psi.expectation_value('Sigmay'))
    center=L//2
    row=psi.correlation_function('Sm','Sp',sites1=[L//4],sites2=list(range(L//4,L-L//4)))[0]
    t2_scale=abs(p['Jx']-p['Jy'])**2/(8*abs(p['Jz'])) if p['Jz'] else None
    final_residual=eng.local_iteration_stats[-1]['max_absolute_residual']
    record=dict(**cfg,couplings=p,tenpy_version=tenpy.__version__,E_shifted=float(E),
                E_unshifted=float(E+model.energy_offset),energy_expectation=float(model.H_MPO.expectation_value(psi)),
                entropy=S,density=n,X=X,Y=Y,bulk_density=float(n[L//4:3*L//4].mean()),
                S_center=float(S[center-1]),S_middle3=float(S[max(0,center-2):center+1].mean()),
                S_middle_third=float(S[L//3:2*L//3].mean()),profile_fits=entropy_fits(S),
                norm_error=float(np.linalg.norm(psi.norm_test())),chi_used=max(psi.chi),
                sweep_stats=eng.sweep_stats,update_stats=eng.update_stats,
                local_solver_diagnostics=dict(scope='Untruncated local solutions; not a global ground-state certificate.',
                    residual_energy_units='literal Pauli Hamiltonian units',
                    effective_t2_energy_scale=t2_scale,
                    scale_definition='Leading blockade scale (Jx-Jy)^2/(8*abs(Jz)); absent for Jz=0.',
                    final_max_residual_over_t2=final_residual/t2_scale if t2_scale else None,
                    residual_threshold_over_t2=residual_tol/t2_scale if t2_scale else None,
                    diag_method=eng.diag_method,strict_local_requested=bool(cfg.get('strict_local')),
                    final_P_tol=float(eng.lanczos_params.get('P_tol',1e-13)),
                    all_updates=local_summary(eng.local_solver_checks,eng.local_residual_tol),
                    final_iteration=eng.local_iteration_stats[-1],iterations=eng.local_iteration_stats,
                    checks=eng.local_solver_checks),
                checkpoint_sweeps=eng.checkpoint_sweeps,correlation_origin=L//4,one_body_row=row,
                runtime_seconds=time.time()-start)
    return plain(record),psi,model


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--name',default='c183');ap.add_argument('--L',type=int,default=32)
    ap.add_argument('--hz',type=float);ap.add_argument('--chi',type=int,default=96)
    ap.add_argument('--sweeps',type=int,default=16);ap.add_argument('--seed',type=int,default=1)
    ap.add_argument('--out',required=True);ap.add_argument('--init-density',type=float,default=.12)
    ap.add_argument('--parity',choices=['even','odd'])
    ap.add_argument('--save-state');ap.add_argument('--resume-state');ap.add_argument('--strict-local',action='store_true')
    ap.add_argument('--checkpoint-every',type=int,default=0)
    ap.add_argument('--local-residual-tol',type=float,default=1e-8,
                    help='Diagnostic absolute local residual threshold; does not set Lanczos P_tol.')
    args=vars(ap.parse_args());out=Path(args.pop('out'))
    logging.basicConfig(level=logging.WARNING)
    rec,_,_=run(args)
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(rec,indent=2)+'\n')
    print(json.dumps({k:rec[k] for k in ('name','L','hz','chi','E_shifted','bulk_density','S_middle3','profile_fits','runtime_seconds')}),flush=True)
