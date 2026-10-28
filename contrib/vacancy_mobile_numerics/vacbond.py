#!/usr/bin/env python3
"""The one-vacancy lift of the route-A demonstrator bond (d = 10 -> d = 11).

Local space C^11 = C|vac> (+) C^10, index 0 = vacancy, index 1 + a = demonstrator letter a
(a = 0, 1: port P;  a = 2..9: B, see ../routeA_demonstrator/cell.py).

    h11 = -t sum_a (|vac, a><a, vac| + h.c.) + J h_occ - (mu_v/2)(n (x) 1 + 1 (x) n),

with h_occ the demonstrator bond on occupied pairs only.  This is the construction of
../counterexample_campaign/doping/README.md (proposed there for the 23-letter bond), applied to the
ten-letter bond.  Conserved: occupied number, the cyclic order of site types (P/B) of the occupied
sites, the logical U(1) charge Q, the conjugation C, and reflection.
Claude (Opus 5.5), 2026-10-28.
"""
import os, sys, importlib.util
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('demo_cell', os.path.join(HERE, '..', 'routeA_demonstrator', 'cell.py'))
C = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(C)

DEMO = dict(D0=1.0, D1=1.0, kap=0.6, Jp=2.0, Jz=-1.0, mu=0.55, K=2.0)
LETTERS = {'P': np.arange(2), 'B': np.arange(2, 10)}          # demonstrator letter indices
Q2_LETTER = np.round(2 * C.charge_site()).astype(int)            # 2q of each demonstrator letter


def demo_bond(lam, **kw):
    p = dict(DEMO)
    p.update(kw)
    return C.bond(lam=lam, **p)


def vacancy_lift(h, t=1.0, J=1.0, mu_v=0.0):
    """Literal lift (same convention as ../counterexample_campaign/doping/bond24.py:vacancy_lift)."""
    q = int(round(np.sqrt(h.shape[0])))
    d = q + 1
    occ = np.array([a * d + b for a in range(1, d) for b in range(1, d)])
    out = np.zeros((d * d, d * d))
    out[np.ix_(occ, occ)] = J * h
    for a in range(1, d):
        out[a, a * d] = out[a * d, a] = -t        # |vac, a><a, vac| + h.c.
    n = np.r_[0, np.ones(q)]
    out -= mu_v / 2 * np.diag(np.add.outer(n, n).reshape(-1))
    return out


def type_bonds(lam, **kw):
    """Restriction of the demonstrator bond to each ordered pair of site types (exact: types are conserved)."""
    h = demo_bond(lam, **kw)
    out = {}
    for a in 'PB':
        for b in 'PB':
            idx = np.array([x * 10 + y for x in LETTERS[a] for y in LETTERS[b]])
            out[a + b] = h[np.ix_(idx, idx)]
    # check type conservation
    tp = np.array([0] * 2 + [1] * 8)
    T2 = np.add.outer(tp * 2, tp).reshape(-1)
    assert np.abs(h * (T2[:, None] != T2[None, :])).max() == 0
    return out, h


if __name__ == '__main__':
    import importlib.util, json
    spec = importlib.util.spec_from_file_location(
        'bond24', os.path.join(HERE, '..', 'counterexample_campaign', 'doping', 'bond24.py'))
    b24 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(b24)
    h = demo_bond(0.15)
    mine = vacancy_lift(h, t=1.3, J=0.27, mu_v=0.4)
    theirs = b24.vacancy_lift(h, hopping=1.3, internal_scale=0.27, chemical_potential=0.4)
    n = np.r_[0, np.ones(10)]
    N2 = np.add.outer(n, n).reshape(-1)
    q2 = np.r_[0, Q2_LETTER]
    Qt = np.add.outer(q2, q2).reshape(-1)
    res = dict(lift_vs_bond24=float(np.abs(mine - theirs).max()),
               hermitian=float(np.abs(mine - mine.T).max()),
               number_comm=float(np.abs(mine * (N2[:, None] - N2[None, :])).max()),
               charge_comm=float(np.abs(mine * (Qt[:, None] - Qt[None, :])).max()),
               norm_h=float(np.linalg.norm(h, 2)), nnz=int(np.count_nonzero(mine)))
    print(json.dumps(res, indent=1))
