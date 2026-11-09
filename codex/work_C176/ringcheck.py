"""C176: consistency of the open-chain DMRG bulk energy with ring ED of (3.2) (C173's ring_ed.py, copied).
Ring ground-state energy per site at n = 1/6 for L = 18, 24, 30, 36 (both fermion boundary conditions), against the
open-chain bulk energy density e = [E(L2,N2) - E(L1,N1)]/(L2 - L1) from DMRG outputs given on the command line."""
import numpy as np, sys, json
import scipy.sparse.linalg as spla
from ring_ed import build_ring
w = 1j
for (L, N) in ((18, 3), (24, 4), (30, 5), (36, 6)):
    es = []
    for bc in (1, -1):
        b, H = build_ring(L, N, w, bc=bc)
        e = spla.eigsh(H, k=1, which='SA', tol=1e-10)[0][0]
        es.append(e / L)
    print(f"ring L={L} N={N} dim={H.shape[0]} E/L(bc=+1)={es[0]:.8f} E/L(bc=-1)={es[1]:.8f} min={min(es):.8f}", flush=True)
if len(sys.argv) > 2:
    r = [json.loads(open(f).readlines()[-1]) for f in sys.argv[1:3]]
    (L1, N1, E1), (L2, N2, E2) = [(x['cfg']['L'], x['cfg']['N'], x['E']) for x in r]
    print(f"open-chain bulk e from L={L1},{L2}: {(E2 - E1) / (L2 - L1):.8f}")
