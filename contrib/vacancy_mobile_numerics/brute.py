#!/usr/bin/env python3
"""Literal check of sqed.py: build H = sum_i h11(i, i+1) on a ring directly in the physical-letter basis,
from the 121 x 121 matrix of vacbond.vacancy_lift, restricted to (M occupied, cyclic type word, 2Q),
and compare its lowest levels with the union of sqed.py's momentum sectors.
Claude (Opus 5.5), 2026-10-28.
"""
import os, sys, json, itertools
import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg import eigsh

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import vacbond as VB  # noqa: E402
import sqed  # noqa: E402


def basis(L, word, Q2):
    M = len(word)
    rots = sorted(set(word[i:] + word[:i] for i in range(M)))
    states = []
    for x in itertools.combinations(range(L), M):
        for w in rots:
            choices = [VB.LETTERS[ty] for ty in w]
            for lab in itertools.product(*choices):
                if sum(VB.Q2_LETTER[a] for a in lab) != Q2:
                    continue
                s = [0] * L
                for j, xi in enumerate(x):
                    s[xi] = 1 + lab[j]
                states.append(tuple(s))
    return states


def hamiltonian(L, states, h11):
    idx = {s: i for i, s in enumerate(states)}
    d = 11
    cols = [np.nonzero(h11[:, k])[0] for k in range(d * d)]
    r, c, v = [], [], []
    for i, s in enumerate(states):
        for b in range(L):
            b1 = (b + 1) % L
            k = s[b] * d + s[b1]
            for o in cols[k]:
                ns = list(s); ns[b], ns[b1] = divmod(o, d)
                j = idx.get(tuple(ns))
                if j is None:
                    raise RuntimeError('leak out of the sector')
                r.append(j); c.append(i); v.append(h11[o, k])
    n = len(states)
    return sps.csr_matrix((v, (r, c)), shape=(n, n))


def check(L, word, Q2, t, J, mu_v, lam, nlev=10):
    h11 = VB.vacancy_lift(VB.demo_bond(lam), t=t, J=J, mu_v=mu_v)
    st = basis(L, word, Q2)
    H = hamiltonian(L, st, h11)
    herm = abs(H - H.T).max()
    if H.shape[0] <= 6000:
        wb = np.linalg.eigvalsh(H.toarray())[:nlev]
    else:   # ask for many more levels than compared, so that ARPACK does not drop degenerate copies
        wb = np.sort(eigsh(H, k=3 * nlev, which='SA', tol=1e-13, ncv=8 * nlev)[0])[:nlev]
    ring = sqed.Ring(L, word, Q2, t, J, mu_v, lam)
    ws = []
    for m in range(L // 2 + 1):
        w = ring.lowest(m, nev=2 * nlev, tol=1e-12)
        ws += list(w)      # the cos projector spans +k and -k, so each +-k pair appears twice
    ws = np.sort(ws)[:nlev]
    return dict(L=L, word=word, Q2=Q2, t=t, J=J, mu_v=mu_v, lam=lam, dim_brute=H.shape[0], dim_sqed=ring.N,
                herm=float(herm), max_diff=float(np.abs(wb - ws).max()), levels=[float(x) for x in wb])


if __name__ == '__main__':
    cases = [(6, 'BPP', 1, 1.0, 0.7, 0.3, 0.15), (7, 'BPPBPP', 0, 1.0, 0.25, 0.1, 0.15),
             (8, 'BPPBPP', 2, 1.0, 0.5, 0.0, 0.10), (7, 'BPPP', 0, 1.0, 0.4, 0.2, 0.15),
             (6, 'BPPBPP', 0, 0.0, 1.0, 0.0, 0.15)]
    out = []
    for cs in cases:
        r = check(*cs)
        print(json.dumps({k: v for k, v in r.items() if k != 'levels'}), flush=True)
        out.append(r)
    json.dump(out, open(os.path.join(HERE, 'data', 'brute_check.json'), 'w'), indent=1)
