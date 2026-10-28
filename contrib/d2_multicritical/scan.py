#!/usr/bin/env python3
"""Two-dimensional ring-ED scan of a qubit family.

usage: scan.py FAMILY --fixed k=v ... --x NAME LO HI N --y NAME LO HI N --L 12 16
               [--k 3] [--procs P] [--out data/scan_xxx.jsonl] [--chunk i/n]

FAMILY is 'kit' (Kitaev variables tau, delta, V, h, D, hx, Jxz; the field enters as
-h Z) or 'xyz' (Jx, Jy, Jz, D, h, hx, Jxz, C15 convention), see params().
For every grid point and size it stores the lowest k levels of every (momentum,
parity) sector, the ground-state sector, and the ring entanglement profile S(l)
of the ground state with its chord-fit slope c_S (l >= L/4).
One JSON line per (point, L).  Resumes: points already in the file are skipped.
"""
import os
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '1')
import sys, json, argparse, time
import multiprocessing as mp
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ring

DEFAULTS = dict(kit=dict(tau=1.0, delta=0.0, V=0.0, h=0.0, D=0.0, hx=0.0, Jxz=0.0, hy=0.0),
                xyz=dict(Jx=0.0, Jy=0.0, Jz=0.0, D=0.0, h=0.0, hx=0.0, Jxz=0.0, hy=0.0))


def params(fam, **kw):
    p = dict(DEFAULTS[fam]); p.update(kw)
    if fam == 'kit':
        Jx = (-p['tau'] + p['delta']) / 2; Jy = (-p['tau'] - p['delta']) / 2; Jz = p['V'] / 4
        J, f = ring.z2_params(Jx, Jy, Jz, p['D'], p['h'], p['hx'], p['Jxz'])
    else:
        J, f = ring.z2_params(p['Jx'], p['Jy'], p['Jz'], p['D'], p['h'], p['hx'], p['Jxz'])
    f = f.copy(); f[1] += p['hy']
    return J, f


def bond_of(fam, **kw):
    return ring.bond(*params(fam, **kw))


def run_point(args):
    fam, kw, L, k = args
    h = bond_of(fam, **kw)
    t0 = time.time()
    sp, key, psi = ring.spectrum(L, h, k=k, want_gs=True)
    S = ring.entropies(psi, L)
    rec = dict(fam=fam, p=kw, L=L, gs=list(key), E0=float(min(np.nanmin(v) for v in sp.values())),
               lev={f'{m},{P}': [float(x) for x in v] for (m, P), v in sp.items()},
               S=S.tolist(), cS=ring.c_chord(S, L), t=time.time() - t0)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('fam')
    ap.add_argument('--fixed', nargs='*', default=[])
    ap.add_argument('--x', nargs=4, required=True)
    ap.add_argument('--y', nargs=4, default=None)
    ap.add_argument('--L', type=int, nargs='+', default=[12, 16])
    ap.add_argument('--k', type=int, default=3)
    ap.add_argument('--procs', type=int, default=os.cpu_count())
    ap.add_argument('--out', required=True)
    ap.add_argument('--chunk', default='0/1')
    a = ap.parse_args()
    fixed = {s.split('=')[0]: float(s.split('=')[1]) for s in a.fixed}
    xs = np.linspace(float(a.x[1]), float(a.x[2]), int(a.x[3]))
    ys = np.linspace(float(a.y[1]), float(a.y[2]), int(a.y[3])) if a.y else [None]
    pts = []
    for y in ys:
        for x in xs:
            kw = dict(fixed); kw[a.x[0]] = round(float(x), 10)
            if y is not None:
                kw[a.y[0]] = round(float(y), 10)
            pts.append(kw)
    ci, cn = map(int, a.chunk.split('/'))
    pts = pts[ci::cn]
    done = set()
    if os.path.exists(a.out):
        for line in open(a.out):
            r = json.loads(line); done.add((json.dumps(r['p'], sort_keys=True), r['L']))
    for L in a.L:
        ring.orbits(L)
    tasks = [(a.fam, kw, L, a.k) for L in a.L for kw in pts if (json.dumps(kw, sort_keys=True), L) not in done]
    print(f'{len(tasks)} tasks', flush=True)
    with mp.get_context('fork').Pool(a.procs) as pool, open(a.out, 'a') as fo:
        for rec in pool.imap_unordered(run_point, tasks, chunksize=1):
            fo.write(json.dumps(rec) + '\n'); fo.flush()


if __name__ == '__main__':
    main()
