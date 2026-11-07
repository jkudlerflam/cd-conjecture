"""Smallest majorant parameters (S, V', K) at finite Trotter number, from data/spec_*.npz.
log Z~(L) = log sum_j (u_j/u_0)^{L/2};  constraint log Z~ <= K + S F(L/(V' beta)) for all even L >= Lmin and all beta.
Prints, for each model and m, V'_min(S,K)/M and where it is attained, plus kappa_eff(K) = sup (L/beta)(log Z~ - K)/M.
"""
import os, sys, json, numpy as np
from scipy.special import logsumexp
from majorant import vmin_needed

HERE = os.path.dirname(os.path.abspath(__file__))
MODELS = {0.0: 0.5, 0.5: 0.5, 1.0: 0.5, 2.0: 0.75}   # Delta: M
VEL = {0.0: 1.0, 0.5: 3 * np.sqrt(3) / 4, 1.0: np.pi / 2, 2.0: None}
LS = np.concatenate([np.arange(2, 41, 2), np.arange(44, 401, 8), [512, 768, 1024, 2048, 4096]])


def logZt_from_u(u, Ls):
    u = np.sort(u)[::-1]
    lu = np.log(np.clip(u, 1e-300, None)) - np.log(u[0])
    return np.array([logsumexp(0.5 * L * lu) for L in Ls])


def load(Delta, m):
    """Delta: float (legacy XXZ files) or registry name"""
    from models import LEGACY
    if isinstance(Delta, str) and Delta in LEGACY:
        Delta = LEGACY[Delta]
    tag = ('D%g' % Delta) if not isinstance(Delta, str) else Delta
    p = os.path.join(HERE, 'data', 'spec_%s_m%d.npz' % (tag, m))
    if not os.path.exists(p):
        return None
    return np.load(p)


def table(Delta, m, Ls=LS):
    z = load(Delta, m)
    if z is None:
        return None
    betas = z['betas']
    g = np.array([logZt_from_u(z['u_%g' % b], Ls) for b in betas]).T   # (L, beta)
    return betas, g


def summary(Delta, mlist, Ss=(2, 4, 6), Ks=(np.log(2), 1.0, 2.0, 3.0), Lmin=2, dmax=None, M=None):
    if M is None:
        M = MODELS[Delta] if not isinstance(Delta, str) else __import__('models').get(Delta)[2]
    rows = []
    for m in mlist:
        t = table(Delta, m)
        if t is None:
            continue
        betas, g = t
        sel = LS >= Lmin
        bsel = np.ones(len(betas), bool) if dmax is None else (betas / m <= dmax)
        LL, BB = np.meshgrid(LS[sel], betas[bsel], indexing='ij')
        gg = g[sel][:, bsel]
        row = dict(m=m)
        for S in Ss:
            for K in Ks:
                V = vmin_needed(gg, LL, BB, S, K)
                i = np.unravel_index(np.argmax(V), V.shape)
                row['S%d_K%.2f' % (S, K)] = (round(V.max() / M, 4), int(LL[i]), float(BB[i]))
        for K in (0.0, np.log(2), 1.0):
            row['kappa_K%.2f' % K] = round(((LL / BB) * (gg - K)).max() / M, 4)
        row['maxlogZt'] = round(gg.max(), 4)
        rows.append(row)
    return rows


if __name__ == '__main__':
    mlist = range(2, 10)
    for Delta in MODELS:
        print('==== Delta =', Delta, ' M =', MODELS[Delta])
        for dmax in [None, 1.0]:
            print('  -- beta/m <=', dmax)
            for r in summary(Delta, mlist, dmax=dmax):
                print('   ', r)
