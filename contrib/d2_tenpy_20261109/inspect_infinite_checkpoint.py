"""Inspect a trusted local MPS checkpoint without changing the running calculation.

The checkpoint can be an unconverged variational state. No ground-state claim is
made by canonicalizing it and evaluating its observables.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import argparse
import json
import pickle
from pathlib import Path
import numpy as np
from independent_infinite_check import local_energy
from run import couplings,plain

ap=argparse.ArgumentParser()
ap.add_argument('--checkpoint',required=True)
ap.add_argument('--name',default='c183')
ap.add_argument('--hz',type=float)
ap.add_argument('--sweep',type=int,required=True)
ap.add_argument('--out',required=True)
a=ap.parse_args()
with open(a.checkpoint,'rb') as handle: psi=pickle.load(handle)
norm_before=float(np.linalg.norm(psi.norm_test()))
psi.canonical_form()
p=couplings(a.name,a.hz)
local=local_energy(psi,p)
S=psi.entanglement_entropy()
n=(1-np.real(psi.expectation_value('Sigmaz')))/2
result=dict(status='provisional_variational_checkpoint_not_converged_ground_state',
    checkpoint=a.checkpoint,checkpoint_sweep=a.sweep,couplings=p,cell=psi.L,chi=max(psi.chi),
    local_energy_density=float(np.real(local.mean())),local_energy_terms=local,
    density=n,mean_density=float(np.mean(n)),density_span=float(np.ptp(n)),
    entropy=S,mean_entropy=float(np.mean(S)),entropy_span=float(np.ptp(S)),
    X=np.real(psi.expectation_value('Sigmax')),Y=np.real(psi.expectation_value('Sigmay')),
    correlation_length=float(psi.correlation_length()),
    norm_before_canonicalization=norm_before,norm_after=float(np.linalg.norm(psi.norm_test())))
Path(a.out).write_text(json.dumps(plain(result),indent=2)+'\n')
print(json.dumps(plain(result),indent=2))
