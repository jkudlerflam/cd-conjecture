"""Independent iDMRG diagnostic: ARPACK, explicit local energy, and Hermiticity.

This is a numerical-method check, not a central-charge measurement.
"""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import argparse
import json
import logging
import pickle
import time
from pathlib import Path
import numpy as np
from run import PauliChain, couplings, plain
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg
from tenpy.linalg import np_conserved as npc


class AuditedEngine(dmrg.TwoSiteDMRGEngine):
    def __init__(self, *args, **kwargs):
        self.checkpoint_path=kwargs.pop('checkpoint_path',None)
        self.hermiticity_checks = []
        self.local_solver_checks = {}
        super().__init__(*args, **kwargs)

    def diag(self, theta_guess):
        # Auditing small actual effective matrices also tests infinite environments.
        size = int(np.prod(theta_guess.shape))
        if size <= 1024 and len(self.hermiticity_checks) < 20:
            h = self.eff_H.to_matrix().to_ndarray()
            skew=h-h.conj().T
            scalar=np.trace(skew)/len(skew)
            self.hermiticity_checks.append(dict(sweep=self.sweeps, size=size,
                relative_defect=float(np.linalg.norm(skew)/max(1.,np.linalg.norm(h))),
                skew_scalar=plain(scalar),
                traceless_skew_norm=float(np.linalg.norm(skew-scalar*np.eye(len(skew))))))
            if len(self.hermiticity_checks)==20:
                print('Initial effective-H maximum relative Hermiticity defect:',
                      max(r['relative_defect'] for r in self.hermiticity_checks),flush=True)
        result=super().diag(theta_guess)
        energy,theta,iterations,_=result
        htheta=self.eff_H.matvec(theta)
        residual=float(npc.norm(htheta-energy*theta))
        rayleigh=npc.inner(theta,htheta,'labels',do_conj=True)/npc.norm(theta)**2
        orthogonal_residual=float(npc.norm(htheta-rayleigh*theta))
        record=self.local_solver_checks.setdefault(str(self.sweeps),dict(max_residual=0.,
            max_orthogonal_residual=0.,max_abs_imag_rayleigh=0.,max_iterations=-1,updates=0))
        record['max_residual']=max(record['max_residual'],residual)
        record['max_orthogonal_residual']=max(record['max_orthogonal_residual'],orthogonal_residual)
        record['max_abs_imag_rayleigh']=max(record['max_abs_imag_rayleigh'],float(abs(np.imag(rayleigh))))
        record['max_iterations']=max(record['max_iterations'],iterations)
        record['updates']+=1
        return result

    def run_iteration(self):
        result=super().run_iteration()
        print(json.dumps({k:plain(self.sweep_stats[k][-1]) for k in ('sweep','E','norm_err','max_trunc_err')}),flush=True)
        print('Last local solves:',self.local_solver_checks.get(str(self.sweeps-1)),flush=True)
        if self.checkpoint_path and self.sweeps%10==0:
            temporary=self.checkpoint_path+'.tmp'
            with open(temporary,'wb') as handle: pickle.dump(self.psi,handle)
            os.replace(temporary,self.checkpoint_path)
        return result


def local_energy(psi, p):
    pairs = [('Sigmax','Sigmax',p['Jx']),('Sigmay','Sigmay',p['Jy']),
             ('Sigmaz','Sigmaz',p['Jz']),('Sigmax','Sigmay',p['Dz']),
             ('Sigmay','Sigmax',-p['Dz']),('Sigmax','Sigmaz',p['Dy']),
             ('Sigmaz','Sigmax',-p['Dy'])]
    terms=[]
    for i in range(psi.L):
        e=p['hz']*psi.expectation_value_term([('Sigmaz',i)])
        for op1,op2,j in pairs:
            if j: e+=j*psi.expectation_value_term([(op1,i),(op2,i+p.get('distance',1))])
        terms.append(e-p['Jz']-p['hz'])
    return np.asarray(terms)


def hardcore_coherent_state(sites, density, phase=-np.pi/2):
    """Uniform exact-exclusion MPS; coefficients are z^N times the specified phase.

    A1 @ A1 = 0, and the classical transfer matrix gives the requested density.
    This is an initial variational state, not an eigenstate of the full spin model.
    """
    if not 0<density<.5:
        raise ValueError('Hardcore density must lie strictly between zero and one half.')
    z=np.sqrt(((1-2*density)**-2-1)/4)
    tensor=np.array([[[1,0],[1,0]],[[0,z*np.exp(1j*phase)],[0,0]]],dtype=complex)
    psi=MPS.from_Bflat(sites,[tensor.copy() for _ in sites],bc='infinite',dtype=complex,
                      form=None,unit_cell_width=len(sites))
    psi.canonical_form()
    return psi


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--name',default='xx')
    ap.add_argument('--chi',type=int,default=16)
    ap.add_argument('--cell',type=int,default=4)
    ap.add_argument('--sweeps',type=int,default=30)
    ap.add_argument('--diag',default='arpack')
    ap.add_argument('--env-init',choices=['iter','contract','TM'],default='contract')
    ap.add_argument('--update-env',type=int,default=1)
    ap.add_argument('--lanczos-max',type=int,default=80)
    ap.add_argument('--mixer',action='store_true')
    ap.add_argument('--seed',type=int,default=921)
    ap.add_argument('--density',type=float)
    ap.add_argument('--phase',type=float,help='Uniform product-state phase; default is random site phases.')
    ap.add_argument('--hardcore',action='store_true',help='Exact blockade coherent-gas MPS instead of product state.')
    ap.add_argument('--warm-L',type=int,default=0)
    ap.add_argument('--warm-sweeps',type=int,default=12)
    ap.add_argument('--load-state',help='Trusted local pickle made by this diagnostic.')
    ap.add_argument('--save-state',help='Temporary MPS checkpoint; do not track in Git.')
    ap.add_argument('--out',required=True)
    a=ap.parse_args(); logging.basicConfig(level=logging.WARNING); start=time.time()
    p=couplings(a.name)
    model=PauliChain(dict(L=a.cell,bc_MPS='infinite',bc_x='periodic',couplings=p))
    rng=np.random.default_rng(a.seed)
    density=a.density if a.density is not None else (.5 if a.name=='xx' else .15)
    warm_record=None
    if a.load_state:
        with open(a.load_state,'rb') as handle: psi=pickle.load(handle)
        if psi.L!=a.cell: raise ValueError('Loaded state unit cell differs from requested cell.')
        psi.canonical_form()
    elif a.warm_L:
        from run import run
        warm_record,finite_psi,_=run(dict(name=a.name,L=a.warm_L,chi=a.chi,sweeps=a.warm_sweeps,
                                        seed=a.seed,init_density=density))
        first=(a.warm_L-a.cell)//2
        tensors=[finite_psi.get_B(i,form='B').transpose(['p','vL','vR']).to_ndarray()
                 for i in range(first,first+a.cell)]
        if tensors[0].shape[1]!=tensors[-1].shape[2]:
            raise ValueError('Central finite tensors have incompatible closing bond dimensions.')
        # The end virtual bases need not be the same: this is a variational initial
        # state, not an exact infinite continuation of the finite-chain solution.
        psi=MPS.from_Bflat(model.lat.mps_sites(),tensors,bc='infinite',dtype=complex,
                          form=None,unit_cell_width=a.cell)
        psi.canonical_form()
    elif a.hardcore:
        phase=-np.pi/2 if a.phase is None else a.phase
        psi=hardcore_coherent_state(model.lat.mps_sites(),density,phase)
    else:
        phases=[a.phase if a.phase is not None else 2*np.pi*rng.random() for _ in range(a.cell)]
        initial=[np.array([np.sqrt(1-density),np.sqrt(density)*np.exp(1j*phase)]) for phase in phases]
        psi=MPS.from_product_state(model.lat.mps_sites(),initial,bc='infinite',dtype=complex,unit_cell_width=a.cell)
    initial_record=dict(local_energy_density=float(np.real(local_energy(psi,p).mean())),
                        density=(1-np.real(psi.expectation_value('Sigmaz')))/2,
                        norm_error=float(np.linalg.norm(psi.norm_test())))
    print('Initial state:',json.dumps(plain(initial_record)),flush=True)
    options=dict(diag_method=a.diag,mixer=a.mixer,max_sweeps=a.sweeps,min_sweeps=10,
        N_sweeps_check=2,max_E_err=1e-10,max_S_err=1e-6,update_env=a.update_env,P_tol_to_trunc=None,
        trunc_params=dict(chi_max=a.chi,svd_min=1e-12),
        lanczos_params=dict(P_tol=1e-11,N_min=30,N_max=a.lanczos_max,N_cache=a.lanczos_max,reortho=True),start_env=4,
        mixer_params=dict(amplitude=1e-5,decay=2.,disable_after=10))
    init_data = dict(start_env_sites=a.cell) if a.env_init=='contract' else dict(force_init_method=a.env_init)
    engine=AuditedEngine(psi,model,options,resume_data=dict(init_env_data=init_data),checkpoint_path=a.save_state)
    returned_E,psi=engine.run()
    norm_before=float(np.linalg.norm(psi.norm_test()))
    psi.canonical_form()
    local=local_energy(psi,p)
    direct=float(np.real(local.mean()))
    result=dict(config=vars(a),couplings=p,returned_energy_density=float(returned_E),
        initial_state=initial_record,warm_finite_record=warm_record,
        local_energy_density=direct,local_terms=local,entropy=psi.entanglement_entropy(),
        norm_before_canonicalization=norm_before,norm_after=float(np.linalg.norm(psi.norm_test())),
        correlation_length=float(psi.correlation_length()),
        density=(1-np.real(psi.expectation_value('Sigmaz')))/2,
        hermiticity_checks=engine.hermiticity_checks,sweep_stats=engine.sweep_stats,
        local_solver_checks=engine.local_solver_checks,
        runtime_seconds=time.time()-start,
        status='diagnostic_only_no_central_charge_claim')
    if a.name=='xx': result['exact_energy_density']=-2/np.pi
    result['passes_all_up_bound']=direct<=1e-9
    result['returned_minus_local']=float(returned_E-direct)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(plain(result),indent=2)+'\n')
    if a.save_state:
        with open(a.save_state,'wb') as handle: pickle.dump(psi,handle)
    print(json.dumps({k:plain(v) for k,v in result.items() if k not in ('sweep_stats','hermiticity_checks')},indent=2))


if __name__=='__main__': main()
