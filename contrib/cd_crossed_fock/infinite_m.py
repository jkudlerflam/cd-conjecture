"""m = infinity (continuous imaginary time) values of log Z~(L) = log Z(L,beta) + L eps0(beta) for small even L.
log Z(L,beta): exact diagonalization of the ring (positive-frame bond, symmetry sectors).
eps0(beta) = beta f(beta): Richardson extrapolation in delta^2 = (beta/m)^2 of the Lanczos values in data/eps0_*.json
(XX checked against the exact p_c). Writes data/infinite_m.json."""
import os, json, glob
import numpy as np
from scipy.special import logsumexp
from scipy.integrate import quad
from cflib import ring_bond_ops
from models import MODELS, get

HERE = os.path.dirname(os.path.abspath(__file__))
BETAS = [1.0, 2.0, 3.0, 4.0, 6.0, 8.0]


def eps0_table(name):
    rows = {}
    for p in glob.glob(os.path.join(HERE, 'data', 'eps0_%s_m*.json' % name)):
        m = int(p.split('_m')[-1].split('.')[0])
        for b, v in json.load(open(p)).items():
            rows.setdefault(float(b), {})[m] = v['eps0']
    return rows


def richardson(ms, vals, beta, npow=3):
    """fit eps0(m) = e_inf + sum_{k=1..npow} a_k delta^(2k) using the largest available m"""
    ms = np.array(ms, float); vals = np.array(vals)
    o = np.argsort(ms)[-(npow + 2):]
    d2 = (beta / ms[o]) ** 2
    A = np.vstack([np.ones_like(d2)] + [d2 ** k for k in range(1, npow + 1)]).T
    coef, res, *_ = np.linalg.lstsq(A, vals[o], rcond=None)
    # error estimate: difference between npow and npow-1 fits
    o2 = np.argsort(ms)[-(npow + 1):]
    d22 = (beta / ms[o2]) ** 2
    A2 = np.vstack([np.ones_like(d22)] + [d22 ** k for k in range(1, npow)]).T
    c2 = np.linalg.lstsq(A2, vals[o2], rcond=None)[0]
    return coef[0], abs(coef[0] - c2[0]), int(ms[o].max())


def sectors(L, mode):
    N = 2 ** L
    idx = np.arange(N)
    bits = (idx[:, None] >> np.arange(L)[None, ::-1]) & 1     # bit j = site j (C order)
    if mode == 'u1':
        stag = (bits * ((-1) ** np.arange(L))[None, :]).sum(1)
        return [idx[stag == q] for q in np.unique(stag)]
    if mode == 'z2':
        par = bits.sum(1) % 2
        return [idx[par == q] for q in (0, 1)]
    return [idx]


def ring_logZ(h, L, betas, mode):
    ops = ring_bond_ops(h, 2, L)
    H = sum(ops).tocsr()
    logs = [[] for _ in betas]
    for s in sectors(L, mode):
        ev = np.linalg.eigvalsh(H[s][:, s].toarray())
        for i, b in enumerate(betas):
            logs[i].append(logsumexp(-b * ev))
    return [float(logsumexp(l)) for l in logs]


def xx_pc(beta):
    return quad(lambda k: np.logaddexp(beta * np.cos(k) / 2, -beta * np.cos(k) / 2), 0, 2 * np.pi, limit=400)[0] / (2 * np.pi)


if __name__ == '__main__':
    out = {}
    for name in MODELS:
        h, mode, M, A, label = get(name)
        tab = eps0_table(name)
        if not tab:
            continue
        Ls = [2, 4, 6, 8, 10, 12] + ([14] if mode == 'u1' else [])
        lz = {L: ring_logZ(h, L, BETAS, mode) for L in Ls}
        res = {}
        for i, b in enumerate(BETAS):
            if b not in tab:
                continue
            ms = sorted(tab[b]); vals = [tab[b][m] for m in ms]
            e_inf, err, mmax = richardson(ms, vals, b)
            g = {L: lz[L][i] + L * e_inf for L in Ls}
            row = dict(eps0_inf=e_inf, eps0_err=err, m_max=mmax, logZt={str(L): g[L] for L in Ls},
                       LoverBeta_logZt_over_M={str(L): L / b * g[L] / M for L in Ls})
            if name == 'XX':
                row['eps0_exact'] = -xx_pc(b)
            res['%g' % b] = row
        out[name] = res
        print(name, {b: (round(r['eps0_inf'], 7), '%.1e' % r['eps0_err'], r.get('eps0_exact')) for b, r in res.items()})
        print('   (L/beta) logZt/M at L=2,4,8:', {b: [round(r['LoverBeta_logZt_over_M'][str(L)], 4) for L in (2, 4, 8)] for b, r in res.items()})
    json.dump(out, open(os.path.join(HERE, 'data', 'infinite_m.json'), 'w'), indent=1)
