"""Check C176's fibdmrg at modified V2 against brute-force ED of the open-chain effective model (L = 14, N = 3)."""
import numpy as np, itertools
from fibdmrg import Chain, DMRG
L, N, wa, V2, V3 = 14, 3, 1.0, -1.0, -2/3
w = 1j*wa; t2 = -1.0
states = [c for c in itertools.combinations(range(L), N) if all(c[i+1]-c[i] >= 2 for i in range(N-1))]
idx = {s: i for i, s in enumerate(states)}
H = np.zeros((len(states),)*2, complex)
for s, i in idx.items():
    occ = set(s)
    H[i, i] = sum(V2*((x+2) in occ) + V3*((x+3) in occ) for x in s)
    for x in s:
        for d, amp in ((1, w), (-1, np.conj(w)), (2, t2), (-2, t2)):
            y = x + d
            if not (0 <= y < L) or y in occ: continue
            new = tuple(sorted((occ - {x}) | {y}))
            if any(new[k+1]-new[k] < 2 for k in range(N-1)): continue
            H[idx[new], i] += amp   # c+_{x+d} c_x with amplitude a_d (hard-core, no passing)
E_ed = np.linalg.eigvalsh(H)[0]
ch = Chain(L, N, w, V2=V2, V3=V3, gauge=1)
dm = DMRG(ch, seed=0)
for chi in (16, 16, 32, 32, 64, 64):
    E, tr = dm.sweep(chi)
print(f"ED {E_ed:.12f}  DMRG {E:.12f}  diff {E-E_ed:.2e}")
