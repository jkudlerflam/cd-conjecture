#!/usr/bin/env python3
"""Checks of ed.py against brute force on small rings: full spectra, translation operator,
and half-ring entropies computed from the explicit d^L state vector."""
import itertools, json, sys
import numpy as np
from ed import Ring, brute_force_check

out = {}
spec = []
for (L, W, nup) in [(4, 2, 1), (5, 2, 1), (5, 4, 2), (6, 4, 2), (6, 2, 0), (6, 4, 1), (7, 4, 2), (7, 6, 3)]:
    n, m, err = brute_force_check(L, W, nup)
    spec.append(dict(L=L, W=W, nup=nup, dim=n, maxerr=err))
    print('spectrum', L, W, nup, n, m, err)
out['spectrum'] = spec

# explicit state vector in the d^L letter basis
def to_letters(r, v):
    L, W = r.L, r.W
    d = 6
    psi = np.zeros(d ** L)
    Psi = v.reshape(r.nc, r.ns)
    for c in range(r.nc):
        pos = r.cp[c]; o = r.co[c]
        # arrows: b_i = right arrow of site i; b_{L-1} = o
        walls = set(pos.tolist())
        b = np.zeros(L, int); cur = o
        for i in range(L):
            if i in walls:
                cur = 1 - cur
            b[i] = cur
        assert b[L - 1] == o
        for s in range(r.ns):
            amp = Psi[c, s]
            if amp == 0:
                continue
            sp = r.sb[s]
            letters = []
            k = 0
            for i in range(L):
                left = b[i - 1]; right = b[i]
                if left == right:
                    letters.append(0 if left == 0 else 1)
                else:
                    up = (sp >> k) & 1; k += 1
                    if (left, right) == (1, 0):
                        letters.append(2 if up else 3)
                    else:
                        letters.append(4 if up else 5)
            idx = 0
            for a in letters:
                idx = idx * d + a
            psi[idx] += amp
    return psi

ent = []
rng = np.random.default_rng(3)
for (L, W, nup, ell) in [(6, 4, 2, 3), (7, 4, 2, 3), (7, 6, 3, 4), (8, 4, 2, 4)]:
    r = Ring(L, W, nup, t=1.0, u=0.6, mu=0.3, J=0.8, Jz=0.2)
    w, V = r.lowest(k=4)
    v = V[:, 0]
    psi = to_letters(r, v)
    M = psi.reshape(6 ** ell, 6 ** (L - ell))
    s = np.linalg.svd(M, compute_uv=False); p = s ** 2; p = p[p > 1e-16]
    S_bf = float(-(p * np.log(p)).sum())
    S_ed = r.half_entropy(v, ell)
    # translation: compare T applied in squeezed coordinates with a cyclic shift of letters
    Tv = to_letters(r, r.translate(v))
    shifted = np.moveaxis(psi.reshape([6] * L), -1, 0).ravel()   # site i -> i+1
    terr = float(np.max(np.abs(Tv - shifted)))
    ent.append(dict(L=L, W=W, nup=nup, ell=ell, S_bruteforce=S_bf, S_ed=S_ed, translate_err=terr))
    print('entropy', L, W, nup, ell, S_bf, S_ed, 'T err', terr)
out['entropy_translation'] = ent
json.dump(out, open('data/check_ed.json', 'w'), indent=1)
