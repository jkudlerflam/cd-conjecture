#!/usr/bin/env python3
"""Assemble the summary table (markdown) from data/: FES (one-cell, or two-cell where the one-cell
uMPS is not injective), open-chain DMRG against the exact Hartree-Fock entropies, translation check.
usage: table.py > data/table.md        Written by Claude (Opus 5.5) for contrib/routeA_n2_chains, 2026-10-28."""
import os, sys, json, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from analyze import cc_fit

D = os.path.join(HERE, 'data')
FULL = ["c28:0.2:full:nu=0.45", "c28:0.5:full:nu=0.45", "c28:1.0:full:nu=0.45",
        "c29:0.2:full", "c29:0.5:full", "c29:0.5:full:nu=1.2", "c29:1.0:full:nu=1.2",
        "spatial:0.05:full", "spatial:0.1:full", "spatial:0.2:full",
        "compile:0.003:full", "compile:0.01:full", "compile:0.01:full:nu=0.9"]
CAL = ["xx:0:free", "c28:0:lead", "c29:0.2:lead", "spatial:0.05:free", "compile:0.003:free"]


def fname(spec, suffix):
    return os.path.join(D, spec.replace(':', '_') + suffix)


def fes_c(spec, k=4, chimin=24):
    """c from the last k accepted chi values (S of the cell averaged over its bonds), the chi range,
    the largest xi, and which run was used."""
    best = None
    for suf in ('_fes.json', '_fes2x.json', '_fes2.json'):
        p = fname(spec, suf)
        if not os.path.exists(p):
            continue
        d = json.load(open(p))
        rows = [r for r in d['results'] if r['chi'] >= chimin and np.isfinite(r['xi']) and r['xi'] < 1e5
                and r['err'] < 1e-3]
        if len(rows) < 3:
            continue
        rows = rows[-k:]
        S = np.array([np.mean(r['S']) for r in rows]); lx = np.log([r['xi'] for r in rows])
        c = 6 * np.polyfit(lx, S, 1)[0]
        cand = dict(c=float(c), chis=(rows[0]['chi'], rows[-1]['chi']), xi=float(rows[-1]['xi']), run=suf[1:-5],
                    n=len(rows))
        if best is None or cand['chis'][1] > best['chis'][1]:
            best = cand
    return best


def trans(spec):
    for suf in ('_fes2x.json', '_fes2.json'):
        p = fname(spec, suf)
        if os.path.exists(p):
            d = json.load(open(p))
            r = [x for x in d['results'] if x['chi'] == 64]
            if not r:
                continue
            r = r[0]
            dens = np.array(r['dens'])
            eb = np.array(r['ebond'])
            return dict(dn=float(np.abs(dens[:, 0] - dens[:, 1]).max()), deb=float(abs(eb[0] - eb[1])),
                        ok=bool(np.isfinite(r['xi']) and r['xi'] < 1e5 and r['err'] < 1e-3))
    return None


def obc(spec):
    p = fname(spec, '_obc.json')
    if not os.path.exists(p):
        return None
    d = json.load(open(p))
    if not d['results']:
        return None
    r = d['results'][-1]
    L = d['L']
    return dict(c=cc_fit(r['S'], L), chi=r['chi'], trunc=r['trunc'], L=L, N=r['N'])


def exact_hf(spec, L=64):
    import run_chain as R
    K, meta = R.build(spec.replace(':full', ':hf') if ':full' in spec else spec)
    ex = R.exact_open(K, L)
    return cc_fit(ex['S'], L), ex['N']


def fmt(x, nd=2):
    return '-' if x is None else (f'{x:.{nd}f}' if isinstance(x, float) else str(x))


if __name__ == '__main__':
    res = {json.loads(l)['spec']: json.loads(l) for l in open(os.path.join(D, 'resonances.jsonl'))}
    print('| chain | HF points, nu | c FES full (chi, xi) | c FES HF | c OBC full (chi) | c OBC HF exact | translation dn, deb (full / HF) |')
    print('|---|---|---|---|---|---|---|')
    for s in FULL:
        h = s.replace(':full', ':hf')
        rr = res.get(s, {})
        f1, f2 = fes_c(s), fes_c(h)
        o = obc(s)
        try:
            (ehf, _), Nh = exact_hf(s)
        except Exception:
            ehf = None
        t1, t2 = trans(s), trans(h)
        tstr = '-' if t1 is None else f"{t1['dn']:.0e}, {t1['deb']:.0e}" + ('' if t1['ok'] else ' (bad)')
        tstr += ' / ' + ('-' if t2 is None else f"{t2['dn']:.0e}, {t2['deb']:.0e}" + ('' if t2['ok'] else ' (bad)'))
        print(f"| {s} | {rr.get('n_points', '-')}, {rr.get('nu_hf', 0):.3f} | "
              f"{fmt(f1 and f1['c'])} ({f1 and f1['chis']}, {f1 and round(f1['xi'])}, {f1 and f1['run']}) | "
              f"{fmt(f2 and f2['c'])} ({f2 and f2['chis']}) | {fmt(o and o['c'][1])} ({o and o['chi']}) | {fmt(ehf)} | {tstr} |")
    print()
    print('| calibration | exact c | c FES (chi, xi) | c OBC DMRG (chi) | c OBC exact L = 64 |')
    print('|---|---|---|---|---|')
    for s in CAL:
        f1 = fes_c(s); o = obc(s)
        try:
            (ex, _), N = exact_hf(s)
        except Exception:
            ex = None
        print(f"| {s} | {1 if s.startswith('xx') else 2} | {fmt(f1 and f1['c'])} ({f1 and f1['chis']}, {f1 and round(f1['xi'])}) | "
              f"{fmt(o and o['c'][1])} ({o and o['chi']}) | {fmt(ex)} |")
