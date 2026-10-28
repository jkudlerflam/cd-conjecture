#!/usr/bin/env python3
"""Exact open-chain entanglement profile of the free target (calibration for the DMRG fits).
On an open chain the hopping phases are a gauge, so each flavour is an open tight-binding chain at
chemical potential mu_a; S(l) = sum_a S_a(l) from the correlation matrix.  The same Calabrese-Cardy
fits as analyze.py obc are applied.
usage: freeobc.py INC 64 96 128 256
Written by Claude (Opus 5.5) for contrib/routeA_incommensurate, 2026-10-28.
"""
import sys, json
import numpy as np
from scipy import linalg
import cell as CE


def profile(par, L):
    t, th, mu, _ = CE.target_params(par)
    S = np.zeros(L - 1)
    for a in range(len(t)):
        h = -t[a] * (np.eye(L, k=1) + np.eye(L, k=-1)) - mu[a] * np.eye(L)
        w, U = linalg.eigh(h)
        occ = U[:, w < 0]
        Cm = occ @ occ.conj().T
        for l in range(1, L):
            nu = linalg.eigvalsh(Cm[:l, :l])
            nu = nu[(nu > 1e-14) & (nu < 1 - 1e-14)]
            S[l - 1] += float(-(nu * np.log(nu) + (1 - nu) * np.log(1 - nu)).sum())
    return S


def cc_fits(S, L, w=8):
    l = np.arange(1, L)
    x = np.log((2 * L / np.pi) * np.sin(np.pi * l / L))
    m = (l >= L // 4) & (l <= 3 * L // 4)
    c_raw = 6 * np.polyfit(x[m], S[m], 1)[0]
    Sa = np.convolve(S, np.ones(w) / w, mode='valid')
    la = np.arange(len(Sa)) + (w + 1) / 2
    xa = np.log((2 * L / np.pi) * np.sin(np.pi * la / L))
    ma = (la >= L // 4) & (la <= 3 * L // 4)
    c_avg = 6 * np.polyfit(xa[ma], Sa[ma], 1)[0]
    return float(c_raw), float(c_avg)


if __name__ == '__main__':
    par = getattr(CE, sys.argv[1])
    out = []
    for L in map(int, sys.argv[2:]):
        S = profile(par, L)
        cr, ca = cc_fits(S, L)
        out.append(dict(L=L, c_raw=cr, c_avg=ca, S=S.tolist()))
        print(sys.argv[1], L, 'c_raw=%.3f c_avg=%.3f Smid=%.4f' % (cr, ca, S[L // 2 - 1]))
    json.dump(out, open(f'data/freeobc_{sys.argv[1]}.json', 'w'))
