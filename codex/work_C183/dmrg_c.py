"""Open-chain DMRG central-charge estimates for a B_2 point (C155, Claude Opus agent, 2026-11-07).

The open chain carries the uniform field f on every site (half from each adjacent bond, plus f/2 at the
two end sites), and no other boundary term.
Estimators from the entanglement profile S(l), l = 1..L-1:
  c_prof : fit S(l) = (c/6) X(l) + s0 + b (-1)^l exp(-...) ~ we use X(l) = log[(2L/pi) sin(pi l/L)] on the
           central window L/8 <= l <= 7L/8, with an explicit (-1)^l term (Calabrese-Cardy plus parity).
  c_sm   : the same fit on the three-point smoothed profile (S(l-1) + 2 S(l) + S(l+1))/4.
  c_two  : 6 [Sm_{L2}(L2/2) - Sm_{L1}(L1/2)] / log(L2/L1) from two chain lengths (smoothed midpoints).
"""
import numpy as np
from est import J_f, bond_from_JB, PAULI
import dmrg as D


def chain_bond(v):
    v = np.asarray(v, float); v = v / np.linalg.norm(v)
    J, f = J_f(v)
    h = bond_from_JB(J, f)
    edge = sum(f[a] * PAULI[a] for a in range(3)) / 2
    return h, edge


def profile_fits(S, L):
    l = np.arange(1, L)
    X = np.log(2 * L / np.pi * np.sin(np.pi * l / L))
    w = (l >= L // 8) & (l <= L - L // 8)
    A = np.stack([X[w] / 6, np.ones(w.sum()), (-1.0) ** l[w]], 1)
    coef, *_ = np.linalg.lstsq(A, S[w], rcond=None)
    Sm = S.copy(); Sm[1:-1] = (S[:-2] + 2 * S[1:-1] + S[2:]) / 4
    w2 = w.copy(); w2[0] = w2[-1] = False
    A2 = np.stack([X[w2] / 6, np.ones(w2.sum())], 1)
    coef2, *_ = np.linalg.lstsq(A2, Sm[w2], rcond=None)
    return float(coef[0]), float(coef2[0]), float(Sm[L // 2 - 1])


def run(v, Ls=(32, 64), chis=(16, 32, 64, 96, 96), seed=0, verbose=False):
    h, edge = chain_bond(v)
    out = dict(L={})
    mids = []
    for L in Ls:
        r = D.dmrg(h, L, chis=chis, seed=seed, edge_field=edge, verbose=verbose)
        S = r['S']
        cp, cs, mid = profile_fits(S, L)
        hist = r['hist']
        out['L'][L] = dict(E=r['E'], c_prof=cp, c_sm=cs, Smid=mid, trunc=hist[-1][2],
                           dE_last=float(hist[-1][1] - hist[-2][1]) if len(hist) > 1 else None,
                           S=[float(x) for x in S])
        mids.append(mid)
    if len(Ls) >= 2:
        out['c_two'] = float(6 * (mids[-1] - mids[0]) / np.log(Ls[-1] / Ls[0]))
    return out


if __name__ == '__main__':
    import time, sys
    tests = {'TFI c=1/2': np.array([-1, 0, 0, 0, 0, 0, 0, 0, -1.0]),
             'XX c=1': np.array([1, 1, 0, 0, 0, 0, 0, 0, 0.0]),
             'XXZ D=-0.5': np.array([1, 1, -.5, 0, 0, 0, 0, 0, 0.0])}
    Ls = tuple(int(a) for a in sys.argv[1:]) or (32, 64)
    for name, v in tests.items():
        t = time.time()
        o = run(v, Ls)
        print(name, {L: (round(d['c_prof'], 3), round(d['c_sm'], 3), f"{d['trunc']:.1e}", d['dE_last']) for L, d in o['L'].items()},
              'c_two', round(o.get('c_two', np.nan), 3), f'{time.time() - t:.0f}s', flush=True)
