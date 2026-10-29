"""Small exact diagonalizations in the invariant zero-excitation sector.

Data bits are spectators at lambda=g=field=0. All off-grammar moves and
both vacancy colors are retained. No truncation by the penalty is made.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import itertools
import json
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

MOVES = {('R', 'R'): ('L', 'R'), ('L', 'R'): ('R', 'R'),
         ('R', 'M'): ('L', 'M'), ('L', 'M'): ('R', 'M')}


def out(a):
    return int(a == 'R' or a == '1')


def inc(a):
    return (int(a),) if a in ('0', '1') else ((0,) if a == 'L' else (0, 1))


def vacancy(a):
    return a in ('V', '0', '1')


def build(n, cells, nv, colored):
    material = (['L']+['R']*(n-1)+['M'])*cells
    seed = tuple(material+(['0'] if colored else ['V'])*nv)
    states, index = [seed], {seed: 0}
    rows, cols, vals = [], [], []
    prows, pcols, pvals = [], [], []
    penalties, badmat = [], []
    cur = 0
    while cur < len(states):
        w = states[cur]
        def emit(i, j, a, b, hop):
            z = list(w)
            z[i], z[j] = a, b
            z = tuple(z)
            if z not in index:
                index[z] = len(states)
                states.append(z)
            rr, cc, vv = (rows, cols, vals) if hop else (prows, pcols, pvals)
            rr.append(index[z]); cc.append(cur); vv.append(-1.)
        pen = 0
        for i, a in enumerate(w):
            j = (i+1) % len(w)
            b = w[j]
            if not vacancy(a) and not vacancy(b):
                if (a, b) in MOVES:
                    c, d = MOVES[a, b]
                    emit(i, j, c, d, False)
                pen += a == 'R' and b == 'L'
            elif colored:
                pen += out(a) not in inc(b)
                if vacancy(a) and not vacancy(b) and int(a) in inc(b):
                    emit(i, j, b, str(out(b)), True)
                elif not vacancy(a) and vacancy(b) and int(b) == out(a):
                    for q in inc(a):
                        emit(i, j, str(q), a, True)
            elif vacancy(a) != vacancy(b):
                emit(i, j, b, a, True)
        penalties.append(pen)
        mat = [a for a in w if not vacancy(a)]
        badmat.append(sum(a == 'R' and mat[(i+1) % len(mat)] == 'L'
                          for i, a in enumerate(mat)))
        cur += 1
        if len(states) > 30000:
            raise RuntimeError('Small-computation cap exceeded')
    shape = (len(states), len(states))
    T = sp.coo_matrix((vals, (rows, cols)), shape=shape).tocsc()
    W = sp.coo_matrix((pvals, (prows, pcols)), shape=shape).tocsc()
    assert (T-T.T).nnz == (W-W.T).nnz == 0
    pen, badmat = np.array(penalties), np.array(badmat)
    if colored:
        assert np.all(badmat <= pen)
    return T, W, pen, badmat, states


def run(T, W, pen, badmat, K, J=1.):
    H = T+J*W+J*K*sp.diags(pen.astype(float), format='csc')
    if H.shape[0] <= 500:
        es, vs = np.linalg.eigh(H.toarray())
        e, v = es[0], vs[:, 0]
    else:
        # A rigorous lower bound follows from the off-diagonal row sum.
        sigma = -float(np.max(np.asarray((abs(T)+J*abs(W)).sum(axis=1))))-1
        es, vs = eigsh(H, k=1, sigma=sigma, which='LM', tol=1e-11,
                       v0=np.ones(H.shape[0]), ncv=30)
        e, v = es[0], vs[:, 0]
    return dict(K=K, J=J, energy=float(e), residual=float(np.linalg.norm(H@v-e*v)),
                illegal_material_weight=float(np.linalg.norm(v[badmat > 0])**2),
                material_bad_expectation=float(np.vdot(v, badmat*v).real),
                local_bad_expectation=float(np.vdot(v, pen*v).real))


def projection_check(n=2, cells=1, nv=2):
    t, w, pen, _, states = build(n, cells, nv, True)
    a, b, p, bad, old = build(n, cells, nv, False)
    good = np.flatnonzero(pen == 0)
    goodold = np.flatnonzero(bad == 0)
    lookup = {tuple('V' if vacancy(x) else x for x in states[i]): i for i in good}
    order = [lookup[old[i]] for i in goodold]
    assert len(lookup) == len(good) == len(goodold)
    dt = t[order][:, order]-a[goodold][:, goodold]
    dw = w[order][:, order]-b[goodold][:, goodold]
    assert dt.nnz == dw.nnz == 0
    return dict(projected_dim=len(good), hopping_difference=0., material_difference=0.)


def main():
    result = dict(projection=projection_check(), runs=[])
    for n, cells, nv in [(2, 1, 2), (3, 2, 1)]:
        for colored in (False, True):
            t, w, pen, bad, states = build(n, cells, nv, colored)
            for K in (10., 100., 1000., 10000.):
                rec = dict(n=n, cells=cells, vacancies=nv, colored=colored, dim=len(states),
                           **run(t, w, pen, bad, K))
                result['runs'].append(rec)
                print(json.dumps(rec), flush=True)
    with open('vacancy_checks.json', 'w') as f:
        json.dump(result, f, indent=2)


if __name__ == '__main__':
    main()
