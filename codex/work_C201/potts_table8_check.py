"""Sharper check of Astra's table (8): Potts spectrum (both signs) is contained in the union of twisted-XXZ
sectors restricted to the PREDICTED finite set: m=0 with phi in {+-2gamma, pi}; m=j>0 with phi=2 pi k/j.
Also check the ground-state equality (10) E_P^{(+)}(L) = min Spec(2L; cos gamma, 2 gamma, 0) to larger L.
Claude Fable 5.2 (checking Astra gpt-6.2-astra), work_C201."""
import numpy as np
from potts_xxz_inclusion import potts_H, xxz_sector
from potts_casimir import ring_E
import scipy.sparse as sp, scipy.sparse.linalg as spl

for Q, Lmax in ((2, 5), (3, 5), (4, 4)):
    gam = np.arccos(np.sqrt(Q)/2)
    for s in (+1, -1):
        for L in range(2, Lmax+1):
            Ep = np.unique(np.round(np.linalg.eigvalsh(potts_H(Q, L, s)), 9))
            N = 2*L; allx = []
            for m in range(0, L+1):
                phis = [2*gam, -2*gam, np.pi] if m == 0 else [2*np.pi*k/m for k in range(m)]
                for phi in phis:
                    allx.append(np.linalg.eigvalsh(xxz_sector(N, s*np.cos(gam), phi, L-m)))
            allx = np.concatenate(allx)
            worst = max(np.min(np.abs(allx-e)) for e in Ep)
            print(f"Q={Q} s={s:+d} L={L}: {len(Ep)} levels, worst mismatch vs predicted finite twist set = {worst:.1e}")
print("ground-state equality (10):")
for Q, Ls in ((3, (6, 7)), (4, (5, 6)), (2, (7,))):
    gam = np.arccos(np.sqrt(Q)/2)
    for L in Ls:
        EP = ring_E(Q, L, +1)
        EX = np.linalg.eigvalsh(xxz_sector(2*L, np.cos(gam), 2*gam, L))[0]
        print(f"Q={Q} L={L}: E_P={EP:.10f}  minSpec(2L;cos g,2g,0)={EX:.10f}  diff={EP-EX:.1e}")
