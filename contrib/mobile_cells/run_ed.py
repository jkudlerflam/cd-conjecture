#!/usr/bin/env python3
"""Ring spectra of the spinful mobile selector (ed.py).

  python3 run_ed.py static --L 9 10 ... --mu -1 --J 0.5 [--pmax 5]
      u = 0.  Enumerates every marker pattern (necklace of cell lengths in [2, pmax]), finds the
      lowest energy in each (S^z = 0 and S^z = 1), and counts the exact ground degeneracy as
      (number of distinct translates of each minimizing pattern) x (degeneracy inside it).
  python3 run_ed.py mobile NAME --L 9 10 ... --u 1 --mu -1 --J 0.5 [--Wwin 1] [--k 8]
      Lowest k levels in every wall-number sector W near the optimum, S^z = 0 and 1, with
      momenta, plus the half-ring entanglement of the ground state.
Output: data/ed/<name>.json.  Claude (Opus 5.5), 2026-10-27.
"""
import os, sys, json, time, argparse
import numpy as np
from ed import Ring

HERE = os.path.dirname(os.path.abspath(__file__))


def necklaces(L, pmin=2, pmax=5):
    out = set()
    def rec(rem, seq):
        if rem == 0:
            m = min(tuple(seq[i:] + seq[:i]) for i in range(len(seq)))
            out.add(m)
            return
        for p in range(pmin, min(pmax, rem) + 1):
            if rem - p == 0 or rem - p >= pmin:
                rec(rem - p, seq + [p])
    rec(L, [])
    return sorted(out)


def translates(L, cells):
    x = np.concatenate([[0], np.cumsum(cells)[:-1]])
    S = frozenset(int(v) for v in x)
    for p in range(1, L + 1):              # smallest shift that maps the marker set to itself
        if L % p == 0 and frozenset((v + p) % L for v in S) == S:
            return p
    return L


def static(a):
    par = dict(t=a.t, u=0.0, mu=a.mu, J=a.J, Jz=a.Jz)
    res = []
    for L in a.L:
        t0 = time.time()
        rows = []
        for cells in necklaces(L, 2, a.pmax):
            M = len(cells)
            ener = {}
            for dS in (0, 1):
                r = Ring(L, 2 * M, M + dS, cells=cells, **par)
                w, _ = r.lowest(k=min(4, r.dim))
                ener[dS] = w.tolist()
            rows.append(dict(cells=list(cells), M=M, translates=translates(L, cells), E=ener))
        Emin = min(min(r_['E'][0][0], r_['E'][1][0]) for r_ in rows)
        tol = 1e-8
        deg = 0
        ground = []
        for r_ in rows:
            n0 = sum(1 for e in r_['E'][0] if abs(e - Emin) < tol)
            n1 = 2 * sum(1 for e in r_['E'][1] if abs(e - Emin) < tol)      # S^z = +1 and -1
            if n0 + n1:
                deg += r_['translates'] * (n0 + n1)
                ground.append(dict(cells=r_['cells'], translates=r_['translates'], n_inside=n0 + n1))
        Es = sorted(set(round(min(r_['E'][0][0], r_['E'][1][0]), 9) for r_ in rows))
        rec = dict(L=L, E0=Emin, degeneracy=deg, ground=ground, next_pattern_gap=(Es[1] - Es[0]) if len(Es) > 1 else None,
                   n_patterns=len(rows), time_s=time.time() - t0)
        res.append(rec)
        print(f"static L={L} E0={Emin:.8f} deg={deg} ground={[(g['cells'], g['translates'], g['n_inside']) for g in ground]} "
              f"next pattern +{rec['next_pattern_gap']:.4f} ({len(rows)} patterns, {rec['time_s']:.0f}s)", flush=True)
    return res


def mobile(a):
    par = dict(t=a.t, u=a.u, mu=a.mu, J=a.J, Jz=a.Jz)
    res = []
    for L in a.L:
        t0 = time.time()
        # optimum W from a coarse look at S^z = 0 ground energies
        Ws = [W for W in range(2, L + 1, 2)]
        if a.density:
            wc = 2 * round(a.density * L / 2)
            Ws = [W for W in Ws if abs(W - wc) <= 2 * a.Wwin]
        levels = []
        gs = {}
        for W in Ws:
            for dS in (0, 1):
                nup = W // 2 + dS
                if nup > W:
                    continue
                r = Ring(L, W, nup, **par)
                if r.dim > a.maxdim:
                    print(f"  skip L={L} W={W} dS={dS} dim={r.dim}", flush=True)
                    continue
                w, V = r.lowest(k=min(a.k, r.dim))
                P = r.momenta(w, V)
                for e, p in zip(w, P):
                    levels.append(dict(W=W, Sz=dS, E=float(e), P=float(p)))
                if dS == 0:
                    gs[W] = (r, w, V[:, 0].copy())
        levels.sort(key=lambda d: d['E'])
        E0 = levels[0]['E']
        # entanglement of the lowest S^z = 0 state of the optimal W
        W0 = levels[0]['W']
        r, w, v0 = gs[W0]
        S = [r.half_entropy(v0, L // 2)]
        if L % 2:
            S.append(r.half_entropy(v0, L // 2 + 1))
        deg = sum(1 + (l['Sz'] > 0) for l in levels if abs(l['E'] - E0) < 1e-8)
        rec = dict(L=L, E0=E0, W0=W0, degeneracy=deg, S_half=S, gap_ground_unique=float(w[1] - w[0]),
                   levels=levels, time_s=time.time() - t0)
        res.append(rec)
        lv = [(l['W'], l['Sz'], round(L * (l['E'] - E0) / (2 * np.pi), 4), l['P']) for l in levels[:10]]
        print(f"mobile L={L} E0={E0:.8f} W0={W0} deg={deg} S={np.round(S, 5).tolist()} "
              f"lowest (W, Sz, L dE/2pi, P): {lv} ({rec['time_s']:.0f}s)", flush=True)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['static', 'mobile'])
    ap.add_argument('name')
    ap.add_argument('--L', type=int, nargs='+', required=True)
    for k, v in dict(t=1.0, u=1.0, mu=-1.0, J=0.5, Jz=0.0).items():
        ap.add_argument('--' + k, type=float, default=v)
    ap.add_argument('--pmax', type=int, default=5)
    ap.add_argument('--k', type=int, default=8)
    ap.add_argument('--density', type=float, default=None, help='expected wall density, to pick W sectors')
    ap.add_argument('--Wwin', type=int, default=1)
    ap.add_argument('--maxdim', type=float, default=6e6)
    a = ap.parse_args()
    res = static(a) if a.mode == 'static' else mobile(a)
    os.makedirs(os.path.join(HERE, 'data', 'ed'), exist_ok=True)
    json.dump(dict(args=vars(a), results=res), open(os.path.join(HERE, 'data', 'ed', a.name + '.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
