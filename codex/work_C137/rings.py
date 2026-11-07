"""Ring energies E_L/L for comparison with window values. Usage: python rings.py NAME Lmin Lmax"""
import sys, json, warnings
warnings.filterwarnings('ignore')
from common import *
import models
def ring_op(h, L):
    H = sum(bond_op(h, j, L) for j in range(L - 1))
    # wrap bond (L-1, 0): permute basis so that it acts on sites L-1 and 0
    d = 2 ** L
    idx = np.arange(d)
    # cyclic shift by one site: new site k = old site k+1
    bits = (idx[:, None] >> np.arange(L - 1, -1, -1)[None, :]) & 1
    rolled = np.roll(bits, -1, axis=1)
    perm = (rolled * (1 << np.arange(L - 1, -1, -1))[None, :]).sum(1)
    P = sp.csr_matrix((np.ones(d), (perm, idx)), shape=(d, d))
    B = bond_op(h, L - 2, L)   # acts on sites L-2, L-1 ; conjugate by shift to act on L-1, 0
    H = H + P.T @ B @ P
    return H.tocsr()
name = sys.argv[1]; Lmin, Lmax = int(sys.argv[2]), int(sys.argv[3])
path = HERE / 'results' / 'rings' ; path.mkdir(exist_ok=True)
out = {}
for L in range(Lmin, Lmax + 1):
    m = models.get(name, L)
    h = m['h'] if np.max(np.abs(m['h'].imag)) > 1e-14 else m['h'].real
    H = ring_op(h, L)
    out[L] = float(lowest(H, 2)[0][0]) / L
    print(name, L, out[L], flush=True)
import os
write_json(path / (f'{name}_L{Lmin}.json' if os.environ.get('C137_SPLIT') else f'{name}.json'), out)
