#!/usr/bin/env python3
"""Momentum-resolved ring ED along a two-parameter cut (ed.py of contrib/d2_classification, unchanged).

For each point and L: the k lowest levels in every momentum sector m (k = 2 pi m / L), the ground-state
momentum, and (for L in --ent) the ground-state block entropies and the chord fit c_S.
usage: ringscan.py CUT [--n 17] [--sizes 12 14 16 18] [--ent 16 18] [--procs P] [--start i --stop j]
Written for contrib/d2_floating_ising (Claude Opus 5.5 subagent, 2026-10-29).
"""
import os, sys, json, time, argparse
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'd2_classification'))
import numpy as np
import ed
import models_fi as MF


def c_chord(S, L):
    ls = np.arange(int(np.ceil(L / 4)), L // 2 + 1)
    x = np.log(L / np.pi * np.sin(np.pi * ls / L)) / 3
    return float(np.polyfit(x, np.asarray(S)[ls - 1], 1)[0])


def run_point(args):
    tag, fam, p, sizes, ent, k = args
    J, f = MF.family(fam, p)
    h = ed.bond(J, f)
    rec = dict(tag=tag, fam=fam, p=p)
    for L in sizes:
        best = None; E = []
        for m in range(L):
            e, u = ed.lowest(ed.sector_hamiltonian(L, h, m), k=k)
            E.append([float(x) for x in e])
            if best is None or e[0] < best[0] - 1e-12:
                best = (e[0], m, u[:, 0])
        r = dict(E=E, m0=int(best[1]))
        if L in ent:
            S = ed.entropies(ed.expand(L, best[1], best[2]), L)
            r['S'] = [float(x) for x in S]; r['cS'] = c_chord(S, L)
        rec[f'L{L}'] = r
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cut')
    ap.add_argument('--n', type=int, default=17)
    ap.add_argument('--sizes', type=int, nargs='+', default=[12, 14, 16, 18])
    ap.add_argument('--ent', type=int, nargs='*', default=[16, 18])
    ap.add_argument('--k', type=int, default=4)
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--start', type=int, default=0); ap.add_argument('--stop', type=int, default=10 ** 9)
    ap.add_argument('--points', default=None, help='JSON list of {"tag","fam","p"} instead of a cut')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    if a.points:
        pts = [(q['tag'], q['fam'], q['p']) for q in json.load(open(a.points))]
    else:
        pts = MF.cut_points(a.cut, a.n)
    pts = pts[a.start:a.stop]
    out = a.out or os.path.join(HERE, 'data', 'ring', f'{a.cut}_{a.start}.jsonl')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    done = set()
    if os.path.exists(out):
        for line in open(out):
            try: done.add(json.loads(line)['tag'])
            except Exception: pass
    tasks = [(t, fam, p, a.sizes, set(a.ent), a.k) for t, fam, p in pts if t not in done]
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
