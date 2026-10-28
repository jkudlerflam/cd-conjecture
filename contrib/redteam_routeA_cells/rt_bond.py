"""Red-team reimplementation of the 23-letter folded-bus bond of
contrib/routeA_spatial/ (written from the README's prose, not from ycell.py).

The bond is assembled as ONE explicit 529 x 529 Hermitian matrix h on C^23 x C^23:
kinetic moves, junction switch, loop link, port, arrow penalty K, and every
on-site term split half/half over the two copies.  numpy only.
"""
import itertools
import numpy as np

SPEC_D = ['v', '0', 'p', 'm']      # stub token species (v = vacuum)
SPEC_E = ['0', 'p', 'm']           # express token species (excited only)
QCH = {'v': 0, '0': 0, 'p': 1, 'm': -1}

LETTERS = ([('l', None, a) for a in (0, 1)] + [('r', None, a) for a in (0, 1)]
           + [('D', s, a) for s in SPEC_D for a in (0, 1)]
           + [('E', s, a) for s in SPEC_E for a in (0, 1)]
           + [('M', None, None)] + [('M', s, None) for s in SPEC_D])
IX = {x: i for i, x in enumerate(LETTERS)}
D0 = len(LETTERS)


def arrows(x):
    k, s, a = x
    if k == 'l':
        return ('L', 'L')
    if k == 'r':
        return ('R', 'R')
    if k in 'DE':
        return ('L', 'R')
    return ('R', 'L') if s is None else ('L', 'L')


def legal_pair(x, y):
    return arrows(x)[1] == arrows(y)[0]


def charge(x):
    k, s, a = x
    return (a or 0) + (QCH[s] if s is not None else 0)


def is_exc(x):
    return x[1] in ('0', 'p', 'm')


def is_marker(x):
    return x[0] == 'M'


WORK = dict(t=1.0, tx=1.0, tsw=1.0, td=1e-3, Delta=0.0, delta=0.0, Dx=4.0, DM=4.0,
            g=0.05, ph_stub={'0': 0.0, 'p': 0.0, 'm': 0.0},
            ph_x={'0': 0.4, 'p': 0.0, 'm': 0.0}, ph_d={'0': 0.9, 'p': 0.0, 'm': 0.0},
            ph_sw={'0': 0.0, 'p': 0.0, 'm': 0.0}, phi=np.pi / 2, eps=0.0, musel=0.0, K=8.0)


def P(**kw):
    out = {k: (dict(v) if isinstance(v, dict) else v) for k, v in WORK.items()}
    out.update(kw)
    return out


def zsign(s, a):
    return (-1) ** a if s == '0' else 1


def onsite_diag(x, p):
    k, s, a = x
    e = 0.0
    if k in 'lrDE':
        e += p['eps'] * a
    if s in ('0', 'p', 'm'):
        e += p['Delta'] + p['delta'] * QCH[s]
        if k == 'E':
            e += p['Dx']
        if k == 'M':
            e += p['DM']
    if k == 'M':
        e += p['musel']
    return e


def onsite_matrix(p):
    """Single-site operator: diagonal energies plus the stub conversions."""
    h1 = np.diag([onsite_diag(x, p) for x in LETTERS]).astype(complex)
    for (x, y) in [(('D', 'p', 0), ('D', '0', 1)), (('D', 'm', 1), ('D', '0', 0))]:
        h1[IX[x], IX[y]] += p['g']
        h1[IX[y], IX[x]] += np.conj(p['g'])
    return h1


def F_ops():
    """Port token operators on one site: F0, F+, F- (acting on markers and stub tokens)."""
    F = {f: np.zeros((D0, D0)) for f in ('0', 'p', 'm')}
    mp = {'0': [('v', '0'), ('0', 'v')], 'p': [('v', 'p'), ('m', 'v')], 'm': [('v', 'm'), ('p', 'v')]}
    for f, pairs in mp.items():
        for (s_in, s_out) in pairs:
            F[f][IX[('M', s_out, None)], IX[('M', s_in, None)]] = 1.0
            for a in (0, 1):
                F[f][IX[('D', s_out, a)], IX[('D', s_in, a)]] = 1.0
    return F


def bond(p, lam=0.0, kinetic=True, port=True, penalty=True, onsite=True):
    """The literal two-site bond h on C^23 x C^23, basis index i*23+j for (left i, right j)."""
    h = np.zeros((D0 * D0, D0 * D0), complex)

    def mv(x, y, x2, y2, c):
        i, j = IX[x] * D0 + IX[y], IX[x2] * D0 + IX[y2]
        h[j, i] += c
        h[i, j] += np.conj(c)
    if kinetic:
        for a in (0, 1):
            for b in (0, 1):
                for s in SPEC_D:
                    ph = 0.0 if s == 'v' else p['ph_stub'][s]
                    mv(('D', s, a), ('r', None, b), ('l', None, a), ('D', s, b),
                       -p['t'] * np.exp(1j * ph) * zsign(s, a))
                for s in SPEC_E:
                    mv(('E', s, a), ('r', None, b), ('l', None, a), ('E', s, b),
                       -p['tx'] * np.exp(1j * p['ph_x'][s]) * zsign(s, a))
            mv(('D', 'v', a), ('M', None, None), ('l', None, a), ('M', 'v', None), -p['t'])
            for s in SPEC_E:
                mv(('E', s, a), ('M', None, None), ('l', None, a), ('M', s, None),
                   -p['tx'] * np.exp(1j * p['ph_x'][s]) * zsign(s, a))
                if p['td'] != 0:
                    mv(('D', s, a), ('M', None, None), ('l', None, a), ('M', s, None),
                       -p['td'] * np.exp(1j * p['ph_d'][s]) * zsign(s, a))
        for m in [('M', None, None)] + [('M', s, None) for s in SPEC_D]:
            for s in SPEC_E:
                for a in (0, 1):
                    mv(m, ('D', s, a), m, ('E', s, a), -p['tsw'] * np.exp(1j * p['ph_sw'][s]))
    if port and lam != 0:
        F = F_ops()
        Pm = np.diag([1.0 if is_marker(x) else 0.0 for x in LETTERS])
        Pd = np.diag([1.0 if x[0] == 'D' else 0.0 for x in LETTERS])
        V = (np.kron(Pm @ F['0'], Pd @ F['0'])
             + np.exp(1j * p['phi']) * np.kron(Pm @ F['p'], Pd @ F['m'])
             + np.exp(-1j * p['phi']) * np.kron(Pm @ F['m'], Pd @ F['p']))
        h += lam * V
    if penalty:
        for x in LETTERS:
            for y in LETTERS:
                if not legal_pair(x, y):
                    i = IX[x] * D0 + IX[y]
                    h[i, i] += p['K']
    if onsite:
        h1 = onsite_matrix(p)
        I = np.eye(D0)
        h += 0.5 * (np.kron(h1, I) + np.kron(I, h1))
    return h


def checks(p):
    out = {}
    out['d0'] = D0
    out['letters_by_kind'] = {k: sum(1 for x in LETTERS if x[0] == k) for k in 'lrDEM'}
    h = bond(p, lam=0.3)
    out['hermitian_err'] = float(np.abs(h - h.conj().T).max())
    q1 = np.diag([charge(x) for x in LETTERS]).astype(float)
    Q2 = np.kron(q1, np.eye(D0)) + np.kron(np.eye(D0), q1)
    out['U1_commutator'] = float(np.abs(h @ Q2 - Q2 @ h).max())
    # marker number conserved on each site
    m1 = np.diag([1.0 if is_marker(x) else 0.0 for x in LETTERS])
    for side, Mx in (('left', np.kron(m1, np.eye(D0))), ('right', np.kron(np.eye(D0), m1))):
        out['marker_conserved_' + side] = float(np.abs(h @ Mx - Mx @ h).max())
    # outer arrows preserved and legal <-> legal by every off-diagonal element
    bad = 0
    hk = bond(p, lam=0.3, penalty=False, onsite=False)
    for i, j in zip(*np.nonzero(np.abs(hk) > 0)):
        x, y = LETTERS[i // D0], LETTERS[i % D0]
        x2, y2 = LETTERS[j // D0], LETTERS[j % D0]
        if arrows(x)[0] != arrows(x2)[0] or arrows(y)[1] != arrows(y2)[1]:
            bad += 1
        if legal_pair(x, y) != legal_pair(x2, y2):
            bad += 1
    out['arrow_violations'] = bad
    # per-site excitation label (vacuum vs excited) conserved by the lam=0 bond
    e1 = np.diag([1.0 if is_exc(x) else 0.0 for x in LETTERS])
    h0 = bond(p, lam=0.0)
    out['exc_conserved_lam0'] = float(np.abs(h0 @ (np.kron(e1, np.eye(D0)) + np.kron(np.eye(D0), e1))
                                             - (np.kron(e1, np.eye(D0)) + np.kron(np.eye(D0), e1)) @ h0).max())
    # is the number of excited tokens at each site individually conserved? (it is not: tokens move)
    # letters actually touched by some term (none redundant)
    used = set()
    for i, j in zip(*np.nonzero(np.abs(bond(p, lam=0.3, penalty=False, onsite=False)) > 0)):
        used.update([i // D0, i % D0, j // D0, j % D0])
    out['letters_unused_by_moves'] = [LETTERS[i] for i in range(D0) if i not in used]
    return out


if __name__ == '__main__':
    import json
    print(json.dumps(checks(P()), default=str, indent=1))
