"""Literal 17-state bond checks and a full nonzero-port four-site closure."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import itertools
import json
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from bond17 import MAT, LET, IX, allowed, outgoing, incoming, params, bond


def main():
    h = bond(params(K=0., musel=0., field=0.), lam=.1)
    number = np.array([int(a[0] != 'V') for a in LET])
    marker = np.array([int(a[0] == 'M') for a in LET])
    charge = np.array([2*(a[2] or 0)+{'': 0, '0': 0, '1': 0, 'A': -1, 'B': 1}[a[1]] for a in LET])
    out = dict(hermiticity=float(np.max(abs(h-h.conj().T))), commutators={})
    for name, q in [('occupied_number', number), ('marker_number', marker), ('charge', charge)]:
        qb = np.add.outer(q, q).ravel()
        out['commutators'][name] = float(np.max(abs((qb[:, None]-qb[None, :])*h)))
    assert out['hermiticity'] == 0.
    assert max(out['commutators'].values()) == 0.
    assert all(allowed(a, b) == (a[0] in ('L', 'M') or b[0] == 'R' or b == ('M', '0', None)) for a in MAT for b in MAT)

    # Every illegal material interval requires a distinct extended bad bond.
    # This exhausts all 17^4 ring words, including all excited letters.
    checked = 0
    for w in itertools.product(range(17), repeat=4):
        material = [i for i in w if i < 15]
        if not material:
            continue
        nb = sum(not allowed(LET[i], LET[w[(j+1) % 4]]) for j, i in enumerate(w))
        bm = sum(not allowed(LET[i], LET[material[(j+1) % len(material)]])
                 for j, i in enumerate(material))
        assert bm <= nb
        checked += 1
    out['four_site_grammar_words_checked'] = checked

    # Nonzero g, lambda, td, and phases. Fixed N=3, one marker, Q=2.
    seed = tuple(IX[x] for x in [('M', '0', None), ('L', '', 0),
                                 ('R', '', 1), ('V', '1', None)])
    states, ix = [seed], {seed: 0}
    transitions = {i: [(j, h[j, i]) for j in np.flatnonzero(h[:, i])] for i in range(289)}
    rows, cols, vals, penalties, badmat = [], [], [], [], []
    cur = 0
    while cur < len(states):
        w = states[cur]
        for i, a in enumerate(w):
            j = (i+1) % 4
            for q, val in transitions[a*17+w[j]]:
                z = list(w)
                z[i], z[j] = divmod(q, 17)
                z = tuple(z)
                if z not in ix:
                    ix[z] = len(states)
                    states.append(z)
                rows.append(ix[z]); cols.append(cur); vals.append(val)
        penalties.append(sum(not allowed(LET[a], LET[w[(i+1) % 4]]) for i, a in enumerate(w)))
        m = [a for a in w if a < 15]
        badmat.append(sum(not allowed(LET[a], LET[m[(i+1) % 3]]) for i, a in enumerate(m)))
        cur += 1
        if len(states) > 20000:
            raise RuntimeError('Small-computation cap exceeded')
    b = sp.coo_matrix((vals, (rows, cols)), shape=(len(states), len(states))).tocsc()
    assert np.max(abs((b-b.getH()).data), initial=0.) < 1e-14
    pen, bad = np.array(penalties), np.array(badmat)
    assert np.all(bad <= pen)
    sigma = -float(np.max(np.asarray(abs(b).sum(axis=1))))-1
    out['full_bond_runs'] = []
    for K in (10., 100., 1000., 10000.):
        H = b+sp.diags(K*pen.astype(float), format='csc')
        es, vs = eigsh(H, k=1, sigma=sigma, which='LM', tol=1e-11,
                       v0=np.ones(len(states)), ncv=35)
        e, v = es[0], vs[:, 0]
        rec = dict(K=K, dim=len(states), energy=float(e),
                   residual=float(np.linalg.norm(H@v-e*v)),
                   illegal_material_weight=float(np.linalg.norm(v[bad > 0])**2),
                   local_bad_expectation=float(np.vdot(v, pen*v).real))
        out['full_bond_runs'].append(rec)
        print(json.dumps(rec), flush=True)
    with open('full_checks.json', 'w') as f:
        json.dump(out, f, indent=2)


if __name__ == '__main__':
    main()
