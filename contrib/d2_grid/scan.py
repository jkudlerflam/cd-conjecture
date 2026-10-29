#!/usr/bin/env python3
"""Ring ED scan of one random great 2-sphere in a residual stratum (or the staggered-field grid).

For every point and L in --sizes: the two lowest levels in every momentum sector (ed.py of
contrib/d2_classification, unchanged), the ground-state block entropies S(l), and the chord fit
c_S(L) = 3 x slope of S(l) against ln[(L/pi) sin(pi l/L)] over l >= L/4.

usage: scan.py sphere STRATUM SEED --N 600 [--sizes 12 16] [--procs P] [--out ...]
       scan.py stag --u0 -1 --u1 1.5 --nu 26 --b0 0.05 --b1 2.5 --nb 25 [--procs P]
       scan.py points FILE.json (a list of {"name":..., "v":[9 numbers]}) [--sizes ...]
Written for contrib/d2_grid (Claude Opus 5.5 subagent, 2026-10-28).
"""
import os, sys, json, time, argparse
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'd2_classification'))
import numpy as np
import ed
import strata


def c_chord(S, L):
    lmin = int(np.ceil(L / 4))
    ls = np.arange(lmin, L // 2 + 1)
    x = np.log(L / np.pi * np.sin(np.pi * ls / L)) / 3
    return float(np.polyfit(x, np.asarray(S)[ls - 1], 1)[0])


def run_point(args):
    tag, v, sizes, k = args
    J, f = strata.J_f(v)
    h = ed.bond(J, f)
    rec = dict(tag=tag, v=[float(x) for x in v])
    for L in sizes:
        best = None
        E = []
        for m in range(L):
            H = ed.sector_hamiltonian(L, h, m)
            e, u = ed.lowest(H, k=k)
            E.append([float(x) for x in e])
            if best is None or e[0] < best[0] - 1e-12:
                best = (e[0], m, u[:, 0])
        psi = ed.expand(L, best[1], best[2])
        S = ed.entropies(psi, L)
        rec[f'L{L}'] = dict(E=E, m0=int(best[1]), S=[float(x) for x in S], cS=c_chord(S, L))
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode')
    ap.add_argument('rest', nargs='*')
    ap.add_argument('--N', type=int, default=600)
    ap.add_argument('--sizes', type=int, nargs='+', default=[12, 16])
    ap.add_argument('--k', type=int, default=2)
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--out', default=None)
    ap.add_argument('--start', type=int, default=0); ap.add_argument('--stop', type=int, default=10 ** 9)
    ap.add_argument('--u0', type=float, default=-1.0); ap.add_argument('--u1', type=float, default=1.6)
    ap.add_argument('--nu', type=int, default=27)
    ap.add_argument('--b0', type=float, default=0.05); ap.add_argument('--b1', type=float, default=2.5)
    ap.add_argument('--nb', type=int, default=25)
    a = ap.parse_args()
    tasks = []
    if a.mode == 'sphere':
        st, seed = a.rest[0], int(a.rest[1])
        X, P = strata.sphere_points(st, seed, a.N)
        tasks = [(f'{st}_{seed}_{i}', P[i], a.sizes, a.k) for i in range(a.N)]
        out = a.out or os.path.join(HERE, 'data', 'scan', f'{st}_{seed}.jsonl')
    elif a.mode == 'stag':
        for u in np.linspace(a.u0, a.u1, a.nu):
            for b in np.linspace(a.b0, a.b1, a.nb):
                tasks.append((f'stag_u{u:.4f}_b{b:.4f}', strata.staggered(u, b), a.sizes, a.k))
        out = a.out or os.path.join(HERE, 'data', 'scan', 'stag.jsonl')
    elif a.mode == 'points':
        pts = json.load(open(a.rest[0]))[a.start:a.stop]
        tasks = [(p['name'], np.array(p['v'], float), a.sizes, a.k) for p in pts]
        suf = '' if a.start == 0 and a.stop >= 10 ** 9 else f'_{a.start}'
        out = a.out or os.path.join(HERE, 'data', 'scan', os.path.basename(a.rest[0]).replace('.json', suf + '.jsonl'))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    done = set()
    if os.path.exists(out):
        for line in open(out):
            try:
                done.add(json.loads(line)['tag'])
            except Exception:
                pass
    tasks = [t for t in tasks if t[0] not in done]
    t0 = time.time()
    with open(out, 'a') as fo:
        if a.procs > 1:
            from multiprocessing import Pool
            with Pool(a.procs) as pool:
                for rec in pool.imap_unordered(run_point, tasks, chunksize=1):
                    fo.write(json.dumps(rec) + '\n'); fo.flush()
        else:
            for t in tasks:
                fo.write(json.dumps(run_point(t)) + '\n'); fo.flush()
    print(f'{len(tasks)} points in {time.time() - t0:.0f} s -> {out}', flush=True)


if __name__ == '__main__':
    main()
