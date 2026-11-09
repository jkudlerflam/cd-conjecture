"""Infinite-chain cross-check with no symmetry constraints; reports finite-chi data."""
import argparse,json,time,logging
from pathlib import Path
import numpy as np
from run import PauliChain,couplings,plain
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg

ap=argparse.ArgumentParser();ap.add_argument('--name',default='c183');ap.add_argument('--cell',type=int,default=4)
ap.add_argument('--chi',type=int,default=64);ap.add_argument('--hz',type=float);ap.add_argument('--seed',type=int,default=1)
ap.add_argument('--sweeps',type=int,default=80);ap.add_argument('--out',required=True)
a=ap.parse_args();logging.basicConfig(level=logging.WARNING);t=time.time()
p=couplings(a.name,a.hz);m=PauliChain(dict(L=a.cell,bc_MPS='infinite',bc_x='periodic',couplings=p))
rng=np.random.default_rng(a.seed)
st=[np.array([np.sqrt(.85),np.sqrt(.15)*np.exp(2j*np.pi*rng.random())]) for _ in range(a.cell)]
psi=MPS.from_product_state(m.lat.mps_sites(),st,bc='infinite',dtype=complex,unit_cell_width=a.cell)
opts=dict(mixer=True,max_sweeps=a.sweeps,min_sweeps=24,N_sweeps_check=2,
          max_E_err=1e-11,max_S_err=1e-7,
          chi_list={0:16,4:min(32,a.chi),8:a.chi},
          trunc_params=dict(chi_max=a.chi,svd_min=1e-12,trunc_cut=1e-11),
          lanczos_params=dict(N_max=80,P_tol=1e-13),
          mixer_params=dict(amplitude=1e-5,decay=1.5,disable_after=16))
engine=dmrg.TwoSiteDMRGEngine(psi,m,opts);E,psi=engine.run();psi.canonical_form()
result=dict(name=a.name,cell=a.cell,chi=a.chi,seed=a.seed,couplings=p,
            E_per_site_shifted=float(E),E_per_site_unshifted=float(E+m.energy_offset/a.cell),
            entropy=psi.entanglement_entropy(),density=(1-np.real(psi.expectation_value('Sigmaz')))/2,
            X=np.real(psi.expectation_value('Sigmax')),Y=np.real(psi.expectation_value('Sigmay')),
            correlation_length=float(psi.correlation_length()),sweep_stats=engine.sweep_stats,
            norm_error=float(np.linalg.norm(psi.norm_test())),runtime_seconds=time.time()-t)
out=Path(a.out);out.parent.mkdir(exist_ok=True,parents=True);out.write_text(json.dumps(plain(result),indent=2)+'\n')
print(json.dumps({k:plain(v) for k,v in result.items() if k!='sweep_stats'}),flush=True)
