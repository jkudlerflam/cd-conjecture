#!/usr/bin/env python3
"""Tables from data/*.json (run_chain.py outputs).

fes  : S (bond entropy), xi (transfer matrix, cells), e (energy per cell) per chi.
       c_k   = 6 x slope of S against ln xi over the last k chi values (k = 4, 6), c_all over all.
       p     = exponent in e(chi) - e_inf ~ xi^{-p} (z = 1 gives p = 2); e_inf exact for quadratic chains,
               fitted (three parameters, last five chi) otherwise.
       dens  = density per cell at the largest chi; nk = jumps of the total occupation n(k) (Fermi points).
fes2 : largest difference between the two cells of a two-cell unit cell, for the density and the bond energy.
obc  : Calabrese-Cardy fit S(l) = (c/6) ln[(2L/pi) sin(pi l/L)] + s0 over L/4 <= l <= 3L/4, raw and after a
       running average over w = 8 bonds; for quadratic chains the same fit of the exact entropies, and the
       exact entropies of the Hartree-Fock chain of the same bond (computed here) for the interacting runs.
usage: analyze.py fes|fes2|obc data/*.json   (prints one JSON line per file)
Written by Claude (Opus 5.5) for contrib/routeA_n2_chains, 2026-10-28.
"""
import os, sys, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import curve_fit


def exact_e_inf(spec, nk=20001):
    import run_chain as R, bands as BD, bonds as B
    K, meta = R.build(spec)
    c, h, U, Hq, H4, r = B.normal_order(K)
    if np.linalg.norm(K - Hq) > 1e-10:
        return None
    T, M = BD.blocks(h)
    ks = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    E = np.array([np.linalg.eigvalsh(BD.bloch(T, M, k)) for k in ks])
    return float(c + E[E < 0].sum() / nk)


def fes(path, exact=True):
    d = json.load(open(path))
    rows = [r for r in d['results']]
    S = np.array([np.mean(r['S']) for r in rows]); xi = np.array([r['xi'] for r in rows])
    e = np.array([r['e'] for r in rows]); chi = [r['chi'] for r in rows]
    lx = np.log(xi)
    out = dict(spec=d['meta']['spec'], chis=chi, S=np.round(S, 4).tolist(), xi=np.round(xi, 2).tolist(),
               err=[float('%.1e' % r['err']) for r in rows], iters=[r['iters'] for r in rows])
    if len(rows) >= 3:
        out['c_pair'] = np.round(6 * np.diff(S) / np.diff(lx), 3).tolist()
        for k in (4, 6):
            if len(rows) >= k:
                out[f'c_{k}'] = round(float(6 * np.polyfit(lx[-k:], S[-k:], 1)[0]), 3)
        out['c_all'] = round(float(6 * np.polyfit(lx, S, 1)[0]), 3)
    out['dens'] = round(float(np.sum(rows[-1]['dens'])), 5)
    out['nk_jumps'] = [(round(k, 3), int(s)) for k, s in rows[-1].get('nk_jumps', [])]
    einf = exact_e_inf(d['meta']['spec']) if exact else None
    if einf is not None:
        de = e - einf
        m = de > 0
        if m.sum() >= 3:
            out['p_exact'] = round(float(-np.polyfit(lx[m][-5:], np.log(de[m][-5:]), 1)[0]), 3)
        out['e_inf_exact'] = einf
    if len(rows) >= 5:
        try:
            f = lambda x, a, b, p: a + b * x ** (-p)
            popt, _ = curve_fit(f, xi[-5:], e[-5:], p0=(e[-1], e[0] - e[-1], 2.0), maxfev=20000)
            out['p_fit'] = round(float(popt[2]), 3)
            out['e_inf_fit'] = float(popt[0])
        except Exception as ex:
            out['p_fit'] = None
    meta = d['meta']
    out['hf_points'] = [(round(p[0], 3), round(p[2], 4)) for p in meta.get('hf_points', meta.get('bare_points', []))]
    return out


def fes2(path):
    d = json.load(open(path))
    out = dict(spec=d['meta']['spec'], rows=[])
    for r in d['results']:
        dens = np.array(r['dens'])            # (2 ops, 2 sites)
        tot = dens.sum(0)
        eb = np.array(r['ebond'])
        out['rows'].append(dict(chi=r['chi'], e=r['e'], dn=float(abs(tot[0] - tot[1])), deb=float(abs(eb[0] - eb[1])),
                                eb=float(eb.mean()), S=np.round(r['S'], 4).tolist(), xi=round(r['xi'], 2)))
    return out


def cc_fit(S, L, w=8):
    S = np.asarray(S)
    l = np.arange(1, L)
    x = np.log((2 * L / np.pi) * np.sin(np.pi * l / L))
    m = (l >= L // 4) & (l <= 3 * L // 4)
    c_raw = 6 * np.polyfit(x[m], S[m], 1)[0]
    Sa = np.convolve(S, np.ones(w) / w, mode='valid')
    la = np.arange(len(Sa)) + (w + 1) / 2
    xa = np.log((2 * L / np.pi) * np.sin(np.pi * la / L))
    ma = (la >= L // 4) & (la <= 3 * L // 4)
    c_avg = 6 * np.polyfit(xa[ma], Sa[ma], 1)[0]
    return round(float(c_raw), 3), round(float(c_avg), 3)


def obc(path):
    import run_chain as R
    d = json.load(open(path))
    L = d['L']
    out = dict(spec=d['meta']['spec'], L=L, rows=[])
    for r in d['results']:
        cr, ca = cc_fit(r['S'], L)
        out['rows'].append(dict(chi=r['chi'], trunc=float('%.1e' % r['trunc']), c_raw=cr, c_avg=ca, N=round(r['N'], 4),
                                Smid=round(r['S'][L // 2 - 1], 4)))
    if 'exact' in d:
        out['exact'] = cc_fit(d['exact']['S'], L)
        out['exact_N'] = d['exact']['N']
    spec = d['meta']['spec']
    if ':full' in spec:
        hspec = spec.replace(':full', ':hf')
        K, meta = R.build(hspec)
        ex = R.exact_open(K, L)
        out['hf_exact'] = cc_fit(ex['S'], L)
        out['hf_N'] = ex['N']
    return out


if __name__ == '__main__':
    mode = sys.argv[1]
    for p in sys.argv[2:]:
        try:
            o = dict(fes=fes, fes2=fes2, obc=obc)[mode](p)
        except Exception as ex:
            o = dict(path=p, error=repr(ex))
        print(json.dumps(o), flush=True)
