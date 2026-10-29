#!/usr/bin/env python3
"""Screening estimator on open chains: two-site DMRG (FDMRG of contrib/d2_followup_220_1282/fdmrg.py,
unchanged) at L = 32 and 64, bond dimension 24 then 48, end sites with half the uniform field (--edgefield: full field).

Open ends pin spirals and remove the rotor zero mode that inflates ring entropies of near-classical
spirals.  With x(l) = (1/6) ln[(2L/pi) sin(pi l/L)] and the window l in [L/4, 3L/4]:
  c_fit(L)  least squares S = c x + b + beta (E_b - mean E_b) over the window (Laflorencie et al.)
  c_dbl     6 [Sbar(64) - Sbar(32)] / ln 2, Sbar = window mean of S (for a CFT exactly (c/6) ln 2 apart)
usage: obcscreen.py --points FILE.json --start i --stop j [--sizes 32 64] [--chis 24 48]
       obcscreen.py --names ising xx ...
Written for contrib/d2_grid (Claude Opus 5.5 subagent, 2026-10-28).
"""
import os, sys, json, time, argparse
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'd2_followup_220_1282'))
import numpy as np
import strata
from vumps import mpo_from_term, bond, PAULI
from fdmrg import FDMRG


def window_fit(S, Eb, L):
    l = np.arange(1, L)
    sel = (l >= L // 4) & (l <= 3 * L // 4)
    x = np.log(2 * L / np.pi * np.sin(np.pi * l / L)) / 6
    A = np.stack([x[sel], np.ones(sel.sum())], 1)
    c_plain = float(np.linalg.lstsq(A, S[sel], rcond=None)[0][0])
    # bond energy of bond l (between sites l-1 and l, 0-based) is Eb[l-1]
    eb = Eb[l[sel] - 1] - Eb[l[sel] - 1].mean()
    A2 = np.stack([x[sel], np.ones(sel.sum()), eb], 1)
    c_eb = float(np.linalg.lstsq(A2, S[sel], rcond=None)[0][0])
    return c_plain, c_eb, float(S[sel].mean())


def run(h, sizes, chis, seed=1, edge=None):
    W = mpo_from_term(h, 2)
    out = {}
    t0 = time.time()
    for L in sizes:
        dm = FDMRG(W, L, seed=seed, chi0=8, edge=edge)
        for chi in chis:
            Eprev = None
            for sw in range(16):
                E, tr = dm.sweep(chi)
                if sw >= 2 and Eprev is not None and abs(E - Eprev) < 1e-9 * L:
                    break
                Eprev = E
        S = dm.S.copy()
        Eb = dm.bond_energies(h)
        cp, ce, Sbar = window_fit(S, Eb, L)
        out[L] = dict(E=E, trunc=tr, S=S.tolist(), c_plain=cp, c_eb=ce, Sbar=Sbar, t=time.time() - t0)
    Ls = sorted(out)
    c_dbl = float(6 * (out[Ls[-1]]['Sbar'] - out[Ls[-2]]['Sbar']) / np.log(Ls[-1] / Ls[-2])) if len(Ls) > 1 else np.nan
    return dict(sizes={str(L): out[L] for L in Ls}, c_dbl=c_dbl)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--points', default=None)
    ap.add_argument('--names', nargs='*', default=[])
    ap.add_argument('--start', type=int, default=0); ap.add_argument('--stop', type=int, default=10 ** 9)
    ap.add_argument('--sizes', type=int, nargs='+', default=[32, 64])
    ap.add_argument('--chis', type=int, nargs='+', default=[24, 48])
    ap.add_argument('--out', default=None)
    ap.add_argument('--edgefield', action='store_true', help='full uniform field on the end sites')
    a = ap.parse_args()
    jobs = []
    if a.points:
        for p in json.load(open(a.points))[a.start:a.stop]:
            J, f = strata.J_f(p['v'])
            jobs.append((p['name'], bond(J, f), f))
        out = a.out or os.path.join(HERE, 'data', 'obcscreen', os.path.basename(a.points).replace('.json', f'_{a.start}.jsonl'))
    else:
        from models_grid import model
        jobs = [(n, model(n)[0], model(n)[4]) for n in a.names]
        out = a.out or os.path.join(HERE, 'data', 'obcscreen', 'named.jsonl')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    done = {json.loads(l)['name'] for l in open(out) if l.strip()} if os.path.exists(out) else set()
    for name, h, f in jobs:
        if name in done:
            continue
        edge = 0.5 * np.einsum('a,aij->ij', np.asarray(f), PAULI) if (a.edgefield and f is not None) else None
        r = run(h, a.sizes, a.chis, edge=edge)
        r['edgefield'] = edge is not None
        r['name'] = name
        with open(out, 'a') as fo:
            fo.write(json.dumps(r) + '\n')
        L = a.sizes[-1]
        print(name, 'c_dbl=%.3f' % r['c_dbl'], 'c_plain(%d)=%.3f' % (L, r['sizes'][str(L)]['c_plain']),
              'c_eb=%.3f' % r['sizes'][str(L)]['c_eb'], 'trunc=%.1e' % r['sizes'][str(L)]['trunc'],
              't=%.0f' % r['sizes'][str(L)]['t'], flush=True)
