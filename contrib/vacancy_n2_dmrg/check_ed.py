#!/usr/bin/env python3
"""Check qdmrg.py against exact diagonalization of open LV_2 chains (and the d = 8 free reference).
usage: check_ed.py   (about a minute).  Claude (Opus 5.5), 2026-10-28."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'vacancy_mobile_numerics')); sys.path.insert(0, HERE)
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import logical_vac as LV
import qdmrg as Q
import models as MD


def ed_open(L, n, M, N, h):
    states, d = LV.sector_states(L, n, M, N)
    # open chain: reuse build_H on bonds 0..L-2 by zeroing the wrap bond with a tiny wrapper
    dim = len(states)
    pw = d ** (L - 1 - np.arange(L, dtype=np.int64))
    digits = ((states[:, None] // pw[None, :]) % d).astype(np.int64)
    rows, cols, vals = [], [], []
    diag = np.zeros(dim)
    hoff = h.copy(); np.fill_diagonal(hoff, 0)
    for i in range(L - 1):
        j = i + 1
        p = digits[:, i] * d + digits[:, j]
        diag += h[p, p]
        for k in np.unique(p):
            sel = np.nonzero(p == k)[0]
            for o in np.nonzero(np.abs(hoff[:, k]) > 1e-14)[0]:
                oi, oj = divmod(o, d)
                new = states[sel] + (oi - digits[sel, i]) * pw[i] + (oj - digits[sel, j]) * pw[j]
                idx = np.searchsorted(states, new)
                assert np.all(states[idx] == new)
                rows.append(idx); cols.append(sel); vals.append(np.full(len(sel), hoff[o, k]))
    ar = np.arange(dim)
    H = sp.csr_matrix((np.concatenate(vals + [diag]), (np.concatenate(rows + [ar]), np.concatenate(cols + [ar]))), shape=(dim, dim))
    w, U = eigsh(H, k=2, which='SA', tol=1e-12)
    o = np.argsort(w)
    S = LV.interval_entropies(L, d, n, states, digits.astype(np.int8), U[:, o[0]])
    return w[o], S


out = []
for (L, M, N, vz) in [(8, 5, [2, 2], None), (9, 6, [2, 1], None), (9, 5, [2, 3], [0, 3])]:
    h = LV.model_bond(2, 1.0, 0.5, 0.0, [1.0, 0.5], [0.0, 0.0], vz)
    w, S = ed_open(L, 2, M, N, h)
    qloc, cfg = MD.lv_charges(2), MD.lv_config(L, M, N)
    dm, hist = Q.run(h, qloc, L, cfg, [64, 256], sweeps=3, log=lambda s: None)
    e = hist[-1]['E']
    Sd = hist[-1]['S'][1:L]
    rec = dict(L=L, M=M, N=N, vz=vz, E_ed=float(w[0]), gap_ed=float(w[1] - w[0]), E_dmrg=e,
               dE=e - float(w[0]), max_dS=float(np.max(np.abs(np.array(Sd) - np.array(S)))))
    print(rec, flush=True); out.append(rec)
json.dump(out, open(os.path.join(HERE, 'data', 'check_ed.json'), 'w'), indent=1)
