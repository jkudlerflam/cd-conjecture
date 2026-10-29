"""Referee rebuild of the C85 14-letter three-colour vacancy bond from
codex/work_C85/PROOF.md Section 1 (o, I table) and eqs. (3)-(4), on top of my11.py.
Does not import bond14.py."""
import numpy as np
import my11

MAT = my11.LETTERS
LET = MAT + ['V0', 'V1', 'V2']
IX = {l: i for i, l in enumerate(LET)}
d = len(LET)
O = {'L0': 0, 'L1': 2, 'R0': 1, 'R1': 1, 'DA': 1, 'DB': 1, 'XA': 1, 'XB': 1, 'E': 0, 'MA': 0, 'MB': 0,
     'V0': 0, 'V1': 1, 'V2': 2}
I = {'L0': (0, 2), 'L1': (0, 2), 'R0': (0, 1, 2), 'R1': (1,), 'DA': (0, 2), 'DB': (0, 2), 'XA': (0, 2),
     'XB': (0, 2), 'E': (0, 1), 'MA': (0,), 'MB': (0,), 'V0': (0,), 'V1': (1,), 'V2': (2,)}


def ext_allowed(a, b):
    return O[a] in I[b]


def is_vac(l):
    return l[0] == 'V'


def bond(pm, lam, J=1., t=1., mu=0., coloured=True):
    """pm: my11 parameter dict (its K is the material penalty K). Returns (h, letters)."""
    ls = LET if coloured else MAT + ['V']
    dd = len(ls); ix = {l: i for i, l in enumerate(ls)}
    h = np.zeros((dd * dd, dd * dd), complex)
    occ = [dd * i + j for i in range(11) for j in range(11)]
    h[np.ix_(occ, occ)] = J * my11.bond(pm, lam)
    for a in MAT:
        if coloured:
            for q in I[a]:
                src = ix[f'V{q}'] * dd + ix[a]; dst = ix[a] * dd + ix[f'V{O[a]}']
                h[dst, src] += -t; h[src, dst] += -t
        else:
            src = ix['V'] * dd + ix[a]; dst = ix[a] * dd + ix['V']
            h[dst, src] += -t; h[src, dst] += -t
    for a in ls:
        for b in ls:
            k = ix[a] * dd + ix[b]
            if coloured and (is_vac(a) or is_vac(b)) and not ext_allowed(a, b):
                h[k, k] += J * pm['K']
            h[k, k] -= mu * ((not is_vac(a)) + (not is_vac(b))) / 2
    return h, ls
