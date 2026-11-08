#!/usr/bin/env python3
"""Low-lying ring spectrum of the staggered chain with momentum resolution, and a ring c estimate.
G = T o (prod_j X_j) (one-site translation times global spin flip) is an exact symmetry of H' that maps
Sz -> -Sz; in the Sz = 0 sector its eigenvalues e^{ip} give a one-site momentum p (G^L = 1 on even rings
since X^L = 1). Levels are reported as x = (E - E0) / (2 pi v / L), v = sqrt(vK * v/K) from ed.luttinger.
c from E0(L)/L = e_inf - pi c v / (6 L^2) using pairs of sizes.
Written for codex/work_C153 (Claude Opus agent, 2026-11-07)."""
import numpy as np, sys, json
from ed import hamiltonian, lowest, luttinger
from scipy.sparse.linalg import eigsh


def G_apply(v, st, L):
    mask = (1 << L) - 1
    s2 = st ^ mask
    s2 = ((s2 << 1) | (s2 >> (L - 1))) & mask   # site j -> j+1
    idx = np.searchsorted(st, s2)
    out = np.zeros_like(v)
    out[idx] = v
    return out


def spectrum(L, u, b, k=10):
    H, st = hamiltonian(L, L // 2, u, b)
    w, V = eigsh(H, k=k, which='SA', tol=1e-11)
    o = np.argsort(w); w, V = w[o], V[:, o]
    # momenta: diagonalize G within (near-)degenerate groups
    ps = np.zeros(len(w)); i = 0
    while i < len(w):
        j = i + 1
        while j < len(w) and abs(w[j] - w[i]) < 1e-7: j += 1
        sub = V[:, i:j]
        Gm = sub.conj().T @ np.column_stack([G_apply(sub[:, a], st, L) for a in range(j - i)])
        ev = np.linalg.eigvals(Gm)
        ps[i:j] = np.sort(np.angle(ev))
        i = j
    return w, ps


if __name__ == '__main__':
    pts = [tuple(map(float, s.split(','))) for s in sys.argv[2:]]
    Ls = [int(x) for x in sys.argv[1].split(',')]
    res = []
    for (u, b) in pts:
        rec = dict(u=u, b=b, sizes=[])
        for L in Ls:
            lt = luttinger(L, u, b)
            w, ps = spectrum(L, u, b)
            e1 = lowest(L, L // 2 + 1, u, b, k=3)
            v = lt['v']; unit = 2 * np.pi * v / L
            rec['sizes'].append(dict(L=L, E0=float(w[0]), K=lt['K'], v=v, x0=[float(x) for x in (w - w[0]) / unit],
                                     p0=[float(p) for p in ps / (2 * np.pi / L)], x1=[float(x) for x in (e1 - w[0]) / unit]))
            print(f'(u,b)=({u},{b}) L={L} K={lt["K"]:.3f} v={v:.4f}')
            print('   Sz=0  x:', ' '.join(f'{x:.3f}[{p:+.1f}]' for x, p in zip(rec['sizes'][-1]['x0'], rec['sizes'][-1]['p0'])))
            print('   Sz=1  x:', ' '.join(f'{x:.3f}' for x in rec['sizes'][-1]['x1']), f' (1/(4K) = {1/(4*lt["K"]):.3f})', flush=True)
        s = rec['sizes']
        for a in range(len(s) - 1):
            L1, L2 = s[a]['L'], s[a + 1]['L']; e1, e2 = s[a]['E0'] / L1, s[a + 1]['E0'] / L2
            v = s[a + 1]['v']
            c = (e2 - e1) / (np.pi * v / 6 * (1 / L1 ** 2 - 1 / L2 ** 2))
            print(f'   c({L1},{L2}) = {c:.3f}')
            rec.setdefault('c_pairs', []).append([L1, L2, float(c)])
        res.append(rec)
    json.dump(res, open(sys.argv[0].replace('.py', '') + '_' + sys.argv[1].replace(',', '_') + '.json', 'w'))
