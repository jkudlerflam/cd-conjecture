"""Independent small checks and Heisenberg/XX baselines, no optimization."""
import json, importlib.util
from pathlib import Path
import numpy as np
from optimize_q import B,G,S,I,ring,energy,evaluate,gauge_M,gauge_bracket
spec=importlib.util.spec_from_file_location('survey','../../contrib/rigorous_toolbox_map/d2_sanity.py')
survey=importlib.util.module_from_spec(spec); spec.loader.exec_module(survey)
rng=np.random.default_rng(131)
errors=[]
for N in (2,3,4,5,6):
 c=rng.normal(size=9); h=np.einsum('a,aij->ij',c,B)
 # The reference orders site 0 as the most significant bit.
 perm=np.array([int(format(s,f'0{N}b')[::-1],2) for s in range(1<<N)])
 H=ring(c,N).toarray(); reference=survey.ring_H(h,N).toarray()[perm][:,perm]
 errors.append(float(np.max(np.abs(H-reference))))
 assert errors[-1]<1e-12
merrors=[]
for _ in range(3):
 c=rng.normal(size=9); h=np.einsum('a,aij->ij',c,B)
 ours,_=gauge_M(h,True); theirs=survey.gauge_M(h)
 merrors.append(abs(ours-theirs)); assert merrors[-1]<1e-7
 shift=rng.normal(size=3); alpha=.718
 changed=h+np.einsum('a,aij->ij',shift,G)+alpha*np.eye(4)
 assert abs(gauge_M(changed,True)[0]-ours)<1e-7
heis=np.r_[np.ones(3),np.zeros(6)]
profiles=[]
for L in (2,3,4,5,6,7,8):
 r=evaluate(heis,L,True); profiles.append(r); print('Heisenberg',L,r['Q'],flush=True)
assert abs(profiles[0]['Q']-4)<1e-10
# C107 ferromagnetic XX: hopping t=1, so Jx=Jy=-2 and M=1.
c=np.r_[[-2.,-2.,0.],np.zeros(6)]; xxerrors=[]
for N in range(2,11):
 exact=-2/np.sin(np.pi/N) if N%2==0 else -1/np.sin(np.pi/(2*N))
 xxerrors.append(abs(energy(c,N)-exact)); assert xxerrors[-1]<1e-8
out=dict(max_ring_matrix_error=max(errors),max_gauge_comparison_error=max(merrors),max_XX_energy_error=max(xxerrors),heisenberg=profiles)
Path('checks_output.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS', {k:v for k,v in out.items() if k!='heisenberg'})
