"""Validate blockk.py: (i) the block decomposition reproduces a direct real-site Hamiltonian on Lb = 4 blocks
(L = 12, open, penalty U across blocks, same terms), (ii) DMRG reproduces its ground energy."""
import numpy as np, itertools
from blockk import block_model, op_matrix, LOC
from dmrg_d import dmrg
w, kappa, V2, V3, mu, U = 1j, 1.0, -1.0, -2/3, -2.1, 40.0
Lb = 4; L = 3*Lb
hb, hs, hl, hr = block_model(w, kappa, V2=V2, V3=V3, mu=mu, U=U)
ons = [hl] + [hs]*(Lb-2) + [hr]
d = len(LOC)
H = np.zeros((d**Lb,)*2, complex)
for b in range(Lb):
    H += np.kron(np.kron(np.eye(d**b), ons[b]), np.eye(d**(Lb-b-1)))
for b in range(Lb-1):
    H += np.kron(np.kron(np.eye(d**b), hb), np.eye(d**(Lb-b-2)))
states = [sum(c, ()) for c in itertools.product(LOC, repeat=Lb)]
Hd = op_matrix(states, L, w, -1.0, V2, V3, mu, U, kappa, across=None, mu_sites=tuple(range(L)), left_end=True, right_end=True)
print("block vs direct max diff:", np.abs(H - Hd).max())
E = np.linalg.eigvalsh(Hd)[0]
out = dmrg(hb, Lb, chis=(16, 32, 64, 64), cutoff=1e-12, onsite=ons)
print(f"ED {E:.12f} DMRG {out['E']:.12f} diff {out['E']-E:.2e}")
