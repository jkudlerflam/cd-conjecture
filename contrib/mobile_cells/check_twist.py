#!/usr/bin/env python3
"""Exact squeezed-space (Ogata-Shiba) check at u = t, J = Jz = 0: the ring spectrum in sector
(W, S^z) equals the union over the squeezed-chain momenta Q = 2 pi n / W of free-fermion
spectra with twisted quantization  k_a = (2 pi I_a + Q + pi (W - 1)) / L,  each with the
multiplicity of C-eigenvalue e^{iQ} on (spins of S^z sector) x (orientation register), where C
shifts the squeezed chain cyclically and flips the orientation.  Writes data/check_twist.json."""
import itertools, json
import numpy as np
from ed import Ring

def predicted(L, W, nup, t=1.0, mu=0.0):
    # multiplicities of C on spins x orientation
    sb = [s for s in range(1 << W) if bin(s).count('1') == nup]
    idx = {s: i for i, s in enumerate(sb)}
    n = len(sb)
    C = np.zeros((2 * n, 2 * n))
    for i, s in enumerate(sb):
        rs = ((s << 1) & ((1 << W) - 1)) | (s >> (W - 1))      # last -> first
        for o in (0, 1):
            C[(1 - o) * n + idx[rs], o * n + i] = 1
    ev = np.linalg.eigvals(C)
    nQ = np.round(np.angle(ev) / (2 * np.pi) * W).astype(int) % W
    mult = np.bincount(nQ, minlength=W)
    E = []
    for q in range(W):
        if mult[q] == 0:
            continue
        Q = 2 * np.pi * q / W
        ks = (2 * np.pi * np.arange(L) + Q + np.pi * (W - 1)) / L
        eps = -2 * t * np.cos(ks)
        for occ in itertools.combinations(range(L), W):
            E += [eps[list(occ)].sum() + mu * W / 2] * mult[q]
    return np.sort(E)

out = []
for (L, W, nup) in [(6, 2, 1), (6, 4, 2), (7, 4, 1), (8, 4, 2), (8, 6, 3), (9, 6, 2)]:
    r = Ring(L, W, nup, t=1.0, u=1.0, mu=0.3, J=0.0, Jz=0.0)
    H = r.sparse().toarray()
    w = np.linalg.eigvalsh(H)
    p = predicted(L, W, nup, mu=0.3)
    err = float(np.max(np.abs(w - p))) if len(w) == len(p) else None
    out.append(dict(L=L, W=W, nup=nup, dim=len(w), npred=len(p), maxerr=err))
    print(L, W, nup, len(w), len(p), err)
json.dump(out, open('data/check_twist.json', 'w'), indent=1)
