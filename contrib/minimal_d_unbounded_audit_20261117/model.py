"""Literal eleven-state boundary-biased shifted-bit variant.

The old C84 bond is recovered at bias=0, q=w=1, field_side='both'.
No projection is performed in bond(); full finite-K cell matrices are
available separately. One inactive data endpoint is always retained.
"""
from itertools import product
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp

LETTERS = ([('L', a) for a in (0, 1)]+[('R', a) for a in (0, 1)]
           +[(r, s) for r in ('d', 'x') for s in ('A', 'B')]
           +[('e', 0), ('m', 'A'), ('m', 'B')])
INDEX = {x: i for i, x in enumerate(LETTERS)}
MARKERS = LETTERS[8:]


def marker(x):
    return x[0] in ('e', 'm')


def excited(x):
    return x[0] in ('d', 'x', 'm')


def charge(x):
    if x[0] in ('L', 'R'):
        return 2*x[1]
    return 0 if x[0] == 'e' else (-1 if x[1] == 'A' else 1)


def allowed(left, right):
    if left == ('L', 0) or marker(left):
        return right != ('R', 1)
    if left == ('L', 1):
        return right in LETTERS[:8] and right != ('R', 1)
    return right in [('R', 0), ('R', 1), ('e', 0)]


def qclock(p, q):
    Q = np.diag([q*q]+[1+q*q]*p)
    return Q-q*(np.eye(p+1, k=1)+np.eye(p+1, k=-1))


def parameters(n, **changes):
    # Diagnostic floating defaults only. For large n or tiny q, supply a
    # separately high-precision musel; double precision loses E_(n+1).
    p = dict(n=n, q=.25, w=.5, D=8., Delta=3., td=.0, g=.05,
             K=100., field=0., field_side='left', bias=1.,
             flux=np.pi/(2*(n+1)), theta=-np.pi/4)
    p.update(changes)
    p.setdefault('musel', 1+p['q']**2-(n+3)*la.eigvalsh(qclock(n+1, p['q']))[0])
    return p


def parts(p):
    size = 11*11
    kinetic = np.zeros((size, size), complex)
    port = np.zeros_like(kinetic)
    onsite = np.zeros(11)

    def add(source, target, value, dest=kinetic):
        a = 11*INDEX[source[0]]+INDEX[source[1]]
        b = 11*INDEX[target[0]]+INDEX[target[1]]
        dest[b, a] += value; dest[a, b] += np.conj(value)

    for bit in (0, 1):
        add((('R', 0), ('R', bit)), (('L', bit), ('R', 0)), -p['q'])
        for s in ('A', 'B'):
            sign = (-1)**bit if s == 'A' else 1
            phi = p['flux'] if s == 'A' else 0.
            add((('d', s), ('R', bit)), (('L', bit), ('d', s)), -p['w']*sign)
            add((('x', s), ('R', bit)), (('L', bit), ('x', s)), -np.exp(1j*phi)*sign)
    add((('R', 0), ('e', 0)), (('L', 0), ('e', 0)), -p['q'])
    for s in ('A', 'B'):
        phi = p['flux'] if s == 'A' else 0.
        psi = (p['n']+1)*phi+np.pi/2
        add((('x', s), ('e', 0)), (('L', 0), ('m', s)), -np.exp(1j*phi))
        add((('d', s), ('e', 0)), (('L', 0), ('m', s)), -p['td']*np.exp(1j*psi))
        for m in MARKERS:
            add((m, ('d', s)), (m, ('x', s)), -1.)
    add((('d', 'A'), ('R', 1)), (('d', 'B'), ('R', 0)), p['g'])
    add((('e', 0), ('R', 0)), (('m', 'A'), ('d', 'B')),
        np.exp(1j*p['theta']), port)
    add((('e', 0), ('R', 0)), (('m', 'B'), ('d', 'A')),
        np.exp(-1j*p['theta']), port)
    for m in MARKERS:
        i = 11*INDEX[m]+INDEX[('R', 0)]
        kinetic[i, i] -= p['bias']
    for i, x in enumerate(LETTERS):
        onsite[i] += p['musel']*marker(x)
        if excited(x):
            onsite[i] += p['Delta']+p['D']*(x[0] in ('x', 'm'))
        if x == ('L', 1) or (p['field_side'] == 'both' and x == ('R', 1)):
            onsite[i] += p['field']
    penalty = np.diag([p['K']*(not allowed(a, b)) for a in LETTERS for b in LETTERS])
    return kinetic, penalty, onsite, port


def bond(p, lam=0.):
    kinetic, penalty, onsite, port = parts(p)
    diagonal = np.diag((onsite[:, None]+onsite[None, :]).ravel()/2)
    return kinetic+penalty+diagonal+lam*port


def words(p):
    if p == 0:
        return [(m,) for m in MARKERS]
    out = []
    for bits in product((0, 1), repeat=p-1):
        for j in range(p):
            for token in [('R', 0), ('d', 'A'), ('d', 'B'), ('x', 'A'), ('x', 'B')]:
                out.append(tuple([('L', a) for a in bits[:j]]+[token]
                                 +[('R', a) for a in bits[j:]]+[('e', 0)]))
        for m in MARKERS:
            out.append(tuple([('L', a) for a in bits]+[('L', 0), m]))
    return out


def vacuum_words(bits):
    p = len(bits)+1
    return [tuple([('L', a) for a in bits[:j]]+[('R', 0)]
                  +[('R', a) for a in bits[j:]]+[('e', 0)]) for j in range(p)] + [
                      tuple([('L', a) for a in bits]+[('L', 0), ('e', 0)])]


def cell_matrix(params, p, full=True):
    states = (list(product(LETTERS[:8], repeat=p)) if full else None)
    states = ([tuple(data)+(m,) for data in states for m in MARKERS]
              if full else words(p))
    index = {w: j for j, w in enumerate(states)}
    kinetic, penalty, onsite, _ = parts(params)
    local = kinetic+penalty
    columns = [[(int(row), local[row, col]) for row in np.flatnonzero(local[:, col])]
               for col in range(121)]
    rows, cols, vals = [], [], []
    fixed = ('e', 0)
    for col, word in enumerate(states):
        rows.append(col); cols.append(col); vals.append(sum(onsite[INDEX[x]] for x in word))
        for j in range(-1, p):
            left = fixed if j == -1 else word[j]
            right = word[j+1]
            for row, value in columns[11*INDEX[left]+INDEX[right]]:
                changed = list(word)
                if j == -1:
                    assert LETTERS[row//11] == fixed
                else:
                    changed[j] = LETTERS[row//11]
                changed[j+1] = LETTERS[row % 11]
                target = tuple(changed)
                if target in index:
                    rows.append(index[target]); cols.append(col); vals.append(value)
                elif full:
                    raise AssertionError((word, target))
    H = sp.coo_matrix((vals, (rows, cols)), shape=(len(states), len(states))).tocsr()
    bad = np.array([sum(not allowed(a, b) for a, b in zip((fixed,)+w, w)) for w in states])
    return states, H, bad
