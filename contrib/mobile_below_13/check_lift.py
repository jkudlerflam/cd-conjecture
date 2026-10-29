"""Small exact diagonalization of the d=12 bond.
 (1) K^-2 suppression of illegal material (squeezed) weight in ground states with
     nonzero bits, vacancies and nonzero port; control = transparent uncoloured
     vacancy (C80/C75-F5 type), which plateaus.
Bit blindness is in check_bitblind.py.
"""
import os, json, itertools
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import bond12 as b
from closure import closure, nbad_material_pairs

HERE = os.path.dirname(os.path.abspath(__file__))
L0, L1, R0, R1, E, V = ('L', 0), ('L', 1), ('R', 0), ('R', 1), ('e', 0), b.V


def ground(H, pen, K):
    A = (H + sp.diags(K * pen, format='csc')).tocsc()
    if A.shape[0] < 1500:
        w, U = np.linalg.eigh(A.toarray()); return w[0], U[:, 0]
    sigma = -float(np.max(np.asarray(abs(H).sum(axis=1)))) - 1
    es, vs = eigsh(A, k=1, sigma=sigma, which='LM', tol=1e-13, ncv=40, v0=np.ones(A.shape[0]))
    v = vs[:, 0]
    assert np.linalg.norm(A @ v - es[0] * v) < 1e-8
    return es[0], v


out = dict(suppression=[])
cases = [('two_cells_nonzero_bits_one_vacancy', [L1, R0, R0, E, L1, R0, R0, E, V], .1),
         ('full_one_vacancy', [L1, R0, E, V], .1),
         ('full_two_vacancies', [L1, R0, E, V, V], .1),
         ('two_cells_two_vacancies_port', [L1, R0, E, L0, R0, E, V, V], .1)]
for name, seed, lam in cases:
    for t in ((1., .1) if len(seed) <= 5 else (1.,)):
        for transparent in (False, True):
            if transparent and len(seed) > 8:
                continue  # control closure too large; the 4-, 5- and 8-site controls suffice
            p = b.params(n=2, K=0.)
            depth = None if (len(seed) <= 5 or transparent) else 2
            H, pen, mb, st = closure(seed, [b.bond(p, lam=lam, t=t, transparent=transparent)] * len(seed),
                                     ill_depth=depth, penfun=(nbad_material_pairs if transparent else None),
                                     cap=600000)
            if not transparent:
                assert np.all(mb <= pen)
            for K in (10., 100., 1000., 10000.):
                e, v = ground(H, pen, K)
                rec = dict(case=name, t=t, control_transparent=transparent, K=K, dimension=len(st),
                           illegal_depth_truncation=depth,
                           energy=float(e), illegal_material_weight=float(np.sum(abs(v[mb > 0]) ** 2)),
                           illegal_total_weight=float(np.sum(abs(v[pen > 0]) ** 2)))
                out['suppression'].append(rec); print(json.dumps(rec), flush=True)
json.dump(out, open(os.path.join(HERE, 'lift_checks.json'), 'w'), indent=1)
