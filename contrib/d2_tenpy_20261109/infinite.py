"""Experimental infinite-chain optimization; numerical stability is not a phase result."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(key,'1')
import argparse,json,time,logging
from pathlib import Path
import numpy as np
from run import PauliChain,couplings,plain
from tenpy.networks.mps import MPS
from independent_infinite_check import local_energy,hardcore_coherent_state,AuditedEngine

ap=argparse.ArgumentParser();ap.add_argument('--name',default='c183');ap.add_argument('--cell',type=int,default=4)
ap.add_argument('--chi',type=int,default=64);ap.add_argument('--hz',type=float);ap.add_argument('--seed',type=int,default=1)
ap.add_argument('--init-density',type=float,default=.15);ap.add_argument('--update-env',type=int,default=20)
ap.add_argument('--hardcore',action='store_true',help='Initialize a coherent gas with exact nearest-neighbor exclusion.')
ap.add_argument('--lanczos-max',type=int,default=80);ap.add_argument('--diag',default='default')
ap.add_argument('--sweeps',type=int,default=80);ap.add_argument('--out',required=True)
a=ap.parse_args();logging.basicConfig(level=logging.WARNING);t=time.time()
p=couplings(a.name,a.hz);m=PauliChain(dict(L=a.cell,bc_MPS='infinite',bc_x='periodic',couplings=p))
rng=np.random.default_rng(a.seed)
st=[np.array([np.sqrt(1-a.init_density),np.sqrt(a.init_density)*np.exp(2j*np.pi*rng.random())]) for _ in range(a.cell)]
psi=(hardcore_coherent_state(m.lat.mps_sites(),a.init_density) if a.hardcore else
     MPS.from_product_state(m.lat.mps_sites(),st,bc='infinite',dtype=complex,unit_cell_width=a.cell))
opts=dict(mixer=True,max_sweeps=a.sweeps,min_sweeps=min(24,a.sweeps),N_sweeps_check=2,
          update_env=a.update_env,P_tol_to_trunc=None,diag_method=a.diag,
          max_E_err=1e-11,max_S_err=1e-7,
          chi_list={0:min(16,a.chi),4:min(32,a.chi),8:a.chi},
          trunc_params=dict(chi_max=a.chi,svd_min=1e-12,trunc_cut=1e-11),
          lanczos_params=dict(N_max=a.lanczos_max,N_cache=a.lanczos_max,reortho=True,P_tol=1e-13),
          mixer_params=dict(amplitude=1e-5,decay=1.5,disable_after=16))
# Contract finite-range environments rather than invoking the iterative/GMRES
# initializer whose imaginary scalar contamination was found in this installation.
engine=AuditedEngine(psi,m,opts,
                            resume_data={'init_env_data':{'start_env_sites':a.cell}})
E,psi=engine.run()
norm_before=float(np.linalg.norm(psi.norm_test()))
psi.canonical_form()
local_terms=local_energy(psi,p)
direct=float(np.real(local_terms.mean()))
stats=engine.sweep_stats
recent_energy=np.asarray(stats['E'][-4:])
recent_entropy=np.asarray(stats['S'][-4:])
last_solve=engine.local_solver_checks[max(engine.local_solver_checks,key=int)]
checks=dict(all_up_variational_bound=direct<=1e-9,
            returned_local_energy_agree=abs(float(E)-direct)<1e-7,
            four_recent_energy_span_small=len(recent_energy)==4 and np.ptp(recent_energy)<1e-8,
            four_recent_entropy_span_small=len(recent_entropy)==4 and np.ptp(recent_entropy)<1e-6,
            final_sweep_norm_small=float(stats['norm_err'][-1])<1e-5,
            final_local_solver_residual_small=last_solve['max_residual']<1e-7)
numerically_stable=all(checks.values())
result=dict(name=a.name,cell=a.cell,chi=a.chi,seed=a.seed,couplings=p,
            initial_state='hardcore_coherent' if a.hardcore else 'random_phase_product',
            initial_density=a.init_density,local_solver_checks=engine.local_solver_checks,
            status='experimental_numerically_stable' if numerically_stable else 'experimental_unconverged',
            ground_state_or_phase_validated=False,convergence_checks=checks,
            initialization='contracted_finite_range_environments',update_env=a.update_env,
            E_per_site_shifted=float(E),E_per_site_unshifted=float(E+m.energy_offset/a.cell),
            local_energy_per_site_shifted=direct,local_energy_per_site_unshifted=direct+m.energy_offset/a.cell,
            local_energy_terms=local_terms,returned_minus_local=float(E)-direct,
            entropy=psi.entanglement_entropy(),density=(1-np.real(psi.expectation_value('Sigmaz')))/2,
            X=np.real(psi.expectation_value('Sigmax')),Y=np.real(psi.expectation_value('Sigmay')),
            correlation_length=float(psi.correlation_length()),sweep_stats=engine.sweep_stats,
            norm_before_final_canonicalization=norm_before,
            norm_error=float(np.linalg.norm(psi.norm_test())),runtime_seconds=time.time()-t)
out=Path(a.out);out.parent.mkdir(exist_ok=True,parents=True);out.write_text(json.dumps(plain(result),indent=2)+'\n')
print(json.dumps({k:plain(v) for k,v in result.items() if k!='sweep_stats'}),flush=True)
