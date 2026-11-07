"""C147: flux pi per bond is the sector minimizer (twist Phi on one bond, Phi=0 is flux pi per bond). L<=10."""
import json
import numpy as np
from itertools import combinations
def ed_twist(L, n, Delta, Phi):
    states = [sum(1 << i for i in c) for c in combinations(range(L), n)]
    idx = {s: i for i, s in enumerate(states)}; H = np.zeros((len(states),)*2, complex)
    for a, s in enumerate(states):
        for i in range(L):
            j = (i+1) % L; bi, bj = (s >> i) & 1, (s >> j) & 1
            H[a, a] += Delta*(0.25 if bi == bj else -0.25)
            if bi != bj:
                t = s ^ ((1 << i) | (1 << j)); ph = np.exp(1j*Phi*(1 if bi else -1)) if i == L-1 else 1.0
                H[idx[t], a] += -0.5*ph
    return np.linalg.eigvalsh(H)[0]
out = {}
for Delta in [0.5, 0.99, 0.999, 1.0, -0.5, -0.99]:
    gap = np.inf
    for L in range(3, 11):
        for n in range(1, L//2+1):
            e0 = ed_twist(L, n, Delta, 0.0)
            for Phi in np.linspace(0.2, 2*np.pi-0.2, 9):
                gap = min(gap, ed_twist(L, n, Delta, Phi)-e0)
    out[str(Delta)] = gap; print(Delta, gap, flush=True)
json.dump(out, open('check_flux.json', 'w'), indent=1)
