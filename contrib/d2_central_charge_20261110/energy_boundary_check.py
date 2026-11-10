"""Post-hoc energy consistency across boundaries; no CFT or gap inference."""
import hashlib,json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
old=HERE/'inputs/energy_consistency.json'
source=HERE/'results/candidate_infinite_L2_P0_s17_chi256_g.json'
prior=next(r for r in json.loads(old.read_text())['cases'] if r['name']=='b_r5_2_d1_10_q0_hm7_8')
bulk=json.loads(source.read_text());assert bulk['converged']
energies={int(k):v for k,v in prior['energies'].items()};sizes=sorted(energies)
pairs=[];fits=[]
for a,b in zip(sizes,sizes[1:]):
    e=(energies[b]-energies[a])/(b-a)
    pairs.append(dict(sizes=[a,b],energy_difference_per_added_site=e,
                      difference_from_infinite=e-bulk['E_physical']))
for i in range(len(sizes)-2):
    window=sizes[i:i+3];x=np.array([[L,1.,1./L] for L in window])
    e,boundary,inv=np.linalg.solve(x,[energies[L] for L in window])
    fits.append(dict(sizes=window,assumed_form='E(L)=e*L+b+a/L',e_bulk=e,
                     e_boundary=boundary,a=inv,difference_from_infinite=e-bulk['E_physical']))
out=dict(registration='Post-hoc read-only consistency check after chi256 completed; not a new c estimator.',
         sources={str(f.relative_to(HERE)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [old,source]},
         infinite_energy_density=bulk['E_physical'],open_size_differences=pairs,conditional_fits=fits,
         interpretation='The finite open-chain differences approach the infinite-chain variational energy. The assumed inverse-size fits do not certify a global minimum, a velocity, a gap, or a central charge.')
(HERE/'energy_boundary_check.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
