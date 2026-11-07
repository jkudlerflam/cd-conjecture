"""Additivity test for two-insertion crossed states.
For one-particle levels a, b (dominant single-insertion eigenstates), build
  |ab> = sum_{r1,r2} e^{2 pi i (n_a r1 + n_b r2)/m} Y_a(j_a + 2 r1) Y_b(j_b + 2 r2) |base>,  slices distinct,
with Y_a the combination of (operator, slice parity) that best reproduces level a, project it to the sector
(Q_a+Q_b, n_a+n_b), remove the vacuum and the single-insertion span W, and record its spectral weights.
Additivity means the weight sits at eps ~ eps_a + eps_b.
Usage: python3 two_insertion.py NAME BETA M [NLEV]"""
import sys, os, json
import numpy as np
from insertion import Sectors, single_insertion, one_particle_levels, OPS1, orth_basis, slice_op
from models import get

HERE = os.path.dirname(os.path.abspath(__file__))


def level_operators(S):
    """for each one-particle level: sector key, energy, list of (op, parity, coefficient)"""
    m = S.m
    vecs = {}
    for a, op in OPS1:
        for j in (0, 1):
            for n in range(m):
                w = sum(np.exp(2j * np.pi * n * r / m) * S.insertion_state([(j + 2 * r, op)]) for r in range(m))
                for key in S.data:
                    if key[1] != n:
                        continue
                    c = S.project(w, key)
                    if key == (0, 0):
                        c[0] = 0
                    if np.linalg.norm(c) > 1e-9 * max(1, np.linalg.norm(w)):
                        vecs.setdefault(key, []).append(((op, j), c))
    levels = []
    for key, lst in vecs.items():
        C = np.array([c for _, c in lst]).T
        Qb = orth_basis(C)
        f = (abs(Qb) ** 2).sum(1)
        for i in np.argsort(-f)[:Qb.shape[1]]:
            e = np.zeros(C.shape[0]); e[i] = 1
            alpha, *_ = np.linalg.lstsq(C, e, rcond=None)
            levels.append(dict(key=key, eps=float(S.eps(key)[i]), f=float(f[i]), ops=[(lst[k][0], alpha[k]) for k in range(len(lst))]))
    levels.sort(key=lambda d: d['eps'])
    return levels


def apply_level(S, lev, vec_fn_base):
    """returns list of (slice, op, weight) terms for the Fourier-transformed insertion of level lev"""
    m = S.m; n = lev['key'][1]
    terms = []
    for (op, j), a in lev['ops']:
        if abs(a) < 1e-12:
            continue
        for r in range(m):
            terms.append((j + 2 * r, op, a * np.exp(2j * np.pi * n * r / m)))
    return terms


def two_state(S, la, lb):
    ta = apply_level(S, la, None); tb = apply_level(S, lb, None)
    v = np.zeros(S.col.N, complex)
    base = S.base
    # group by slice of first insertion to reuse work
    for (ja, opa, wa) in ta:
        va = slice_op(S.col, base, ja, opa)
        for (jb, opb, wb) in tb:
            if jb == ja:
                continue
            v += wa * wb * slice_op(S.col, va, jb, opb)
    return S.T.Ke_half @ v


def run(S, nlev=6):
    mode, m = S.mode, S.m
    prepare_W(S)
    W = W_basis
    levels = level_operators(S)
    mod = 2 if mode == 'z2' else None
    out = []
    low = levels[:nlev]
    for i, la in enumerate(low):
        for lb in low[i:]:
            Q = la['key'][0] + lb['key'][0]
            if mod: Q %= mod
            n = (la['key'][1] + lb['key'][1]) % m
            key = (Q, n)
            if key not in S.data:
                continue
            v = two_state(S, la, lb)
            c = S.project(v, key)
            if key == (0, 0):
                c[0] = 0
            if key in W:          # remove the single-insertion span W of this sector
                B = W[key]
                c = c - B @ (B.conj().T @ c)
            nrm = np.linalg.norm(c) ** 2
            if nrm < 1e-20:
                continue
            p = abs(c) ** 2 / nrm
            e = S.eps(key)
            target = la['eps'] + lb['eps']
            jmax = int(np.argmax(p))
            out.append(dict(a=(int(la['key'][0]), int(la['key'][1]), round(la['eps'], 4)), b=(int(lb['key'][0]), int(lb['key'][1]), round(lb['eps'], 4)),
                            sector=(int(key[0]), int(key[1])), sum_eps=target, mean_eps=float((p * e).sum()),
                            dominant_eps=float(e[jmax]), dominant_weight=float(p[jmax]),
                            weight_below_90pct=float(p[e < 0.9 * target].sum()), lowest_eps_weight_gt_1pct=float(e[p > 0.01].min())))
    return S, levels, out


W_basis = {}


def prepare_W(S):
    global W_basis
    W_basis = {}
    m = S.m
    vecs = {}
    for a, op in OPS1:
        for j in (0, 1):
            v = S.insertion_state([(j, op)])
            for key in S.data:
                c = S.project(v, key)
                if key == (0, 0):
                    c[0] = 0
                if np.linalg.norm(c) > 1e-10 * np.linalg.norm(v):
                    vecs.setdefault(key, []).append(c)
    for key, lst in vecs.items():
        W_basis[key] = orth_basis(np.array(lst).T)


if __name__ == '__main__':
    name, beta, m = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
    nlev = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    h, mode, M, A, label = get(name)
    S = Sectors(name, beta, m, Qmax=3, h=h, mode=mode)
    S, levels, out = run(S, nlev)
    print(name, 'beta', beta, 'm', m, ' lowest one-particle levels (Q,n,eps,f):',
          [(int(l['key'][0]), int(l['key'][1]), round(l['eps'], 3), round(l['f'], 3)) for l in levels[:nlev]])
    for r in out:
        print('  a', r['a'], 'b', r['b'], '-> sector', r['sector'], ' eps_a+eps_b=%.3f  mean=%.3f  dominant=%.3f (w=%.2f)  w(<0.9 sum)=%.3f  lowest(w>1%%)=%.3f' % (
            r['sum_eps'], r['mean_eps'], r['dominant_eps'], r['dominant_weight'], r['weight_below_90pct'], r['lowest_eps_weight_gt_1pct']))
    json.dump(dict(levels=[dict(key=[int(x) for x in l['key']], eps=l['eps'], f=l['f']) for l in levels], pairs=out),
              open(os.path.join(HERE, 'data', 'two_%s_b%g_m%d.json' % (name, beta, m)), 'w'), indent=1)
