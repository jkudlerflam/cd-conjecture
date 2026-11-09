"""Independent full-spin TeNPy DMRG; no number/parity projection or effective truncation.

Pauli convention, uniform fields, open boundaries. A scalar is subtracted so
the all-up product has energy zero; E_unshifted restores the literal chain.
"""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(key, '1')
import argparse
import json
import logging
import time
from pathlib import Path
import numpy as np
import tenpy
from tenpy.models.model import CouplingMPOModel
from tenpy.networks.site import SpinHalfSite
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg


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
        return SpinHalfSite(conserve=None)

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
    return x


def run(cfg):
    start=time.time(); L=cfg['L']; p=couplings(cfg['name'],cfg.get('hz'))
    model=PauliChain(dict(L=L,bc_MPS='finite',bc_x='open',couplings=p))
    rng=np.random.default_rng(cfg.get('seed',1))
    density=cfg.get('init_density',.12)
    # Random complex product state provides a nonzero overlap with all sectors.
    state=[np.array([np.sqrt(1-density),np.sqrt(density)*np.exp(2j*np.pi*rng.random())]) for _ in range(L)]
    psi=MPS.from_product_state(model.lat.mps_sites(),state,bc='finite',dtype=complex)
    chi=cfg.get('chi',96); sweeps=cfg.get('sweeps',16)
    options=dict(mixer=True,max_sweeps=sweeps,min_sweeps=min(8,sweeps),
                 N_sweeps_check=1,max_E_err=1e-11,max_S_err=1e-7,
                 trunc_params=dict(chi_max=chi,svd_min=1e-12,trunc_cut=1e-12),
                 chi_list={0:min(16,chi),2:min(32,chi),4:min(64,chi),6:chi},
                 lanczos_params=dict(N_max=80,P_tol=1e-13),
                 mixer_params=dict(amplitude=1e-5,decay=2.,disable_after=8))
    eng=dmrg.TwoSiteDMRGEngine(psi,model,options)
    E,psi=eng.run()
    psi.canonical_form()
    S=psi.entanglement_entropy()
    n=(1-np.real(psi.expectation_value('Sigmaz')))/2
    X=np.real(psi.expectation_value('Sigmax')); Y=np.real(psi.expectation_value('Sigmay'))
    center=L//2
    row=psi.correlation_function('Sm','Sp',sites1=[L//4],sites2=list(range(L//4,L-L//4)))[0]
    record=dict(**cfg,couplings=p,tenpy_version=tenpy.__version__,E_shifted=float(E),
                E_unshifted=float(E+model.energy_offset),energy_expectation=float(model.H_MPO.expectation_value(psi)),
                entropy=S,density=n,X=X,Y=Y,bulk_density=float(n[L//4:3*L//4].mean()),
                S_center=float(S[center-1]),S_middle3=float(S[max(0,center-2):center+1].mean()),
                S_middle_third=float(S[L//3:2*L//3].mean()),profile_fits=entropy_fits(S),
                norm_error=float(np.linalg.norm(psi.norm_test())),chi_used=max(psi.chi),
                sweep_stats=eng.sweep_stats,correlation_origin=L//4,one_body_row=row,
                runtime_seconds=time.time()-start)
    return plain(record),psi,model


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--name',default='c183');ap.add_argument('--L',type=int,default=32)
    ap.add_argument('--hz',type=float);ap.add_argument('--chi',type=int,default=96)
    ap.add_argument('--sweeps',type=int,default=16);ap.add_argument('--seed',type=int,default=1)
    ap.add_argument('--out',required=True);ap.add_argument('--init-density',type=float,default=.12)
    args=vars(ap.parse_args());out=Path(args.pop('out'))
    logging.basicConfig(level=logging.WARNING)
    rec,_,_=run(args)
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(rec,indent=2)+'\n')
    print(json.dumps({k:rec[k] for k in ('name','L','hz','chi','E_shifted','bulk_density','S_middle3','profile_fits','runtime_seconds')}),flush=True)
