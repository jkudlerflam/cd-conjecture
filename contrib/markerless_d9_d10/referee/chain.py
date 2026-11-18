"""Literal many-site Hamiltonians from a two-site bond matrix (referee code).

A system is a set of words (tuples of letter indices) closed under the nonzero
off-diagonal entries of the bond, on an open chain or a ring.  Nothing about
cells, tokens or sectors is assumed when the matrix is built: every matrix
element is  sum_i <x_i' x_{i+1}'| h |x_i x_{i+1}>  taken from the 2-site matrix.
"""
import numpy as np
import scipy.sparse as sp


def bond_table(h, d, tol=1e-14):
    """dict (r,s) -> list of ((p,q), <p q|h|r s>) for nonzero entries"""
    tab = {}
    rows, cols = np.nonzero(np.abs(h) > tol)
    for row, col in zip(rows, cols):
        tab.setdefault((col // d, col % d), []).append(((row // d, row % d), h[row, col]))
    return tab


def bonds_of(N, periodic):
    b = [(i, i + 1) for i in range(N - 1)]
    if periodic and N > 1:
        b.append((N - 1, 0))
    return b


def closure(seeds, h_generic, d, periodic):
    """all words reachable from the seeds by nonzero off-diagonal bond entries"""
    tab = bond_table(h_generic, d)
    seen = {}
    stack = []
    for s in seeds:
        s = tuple(s)
        if s not in seen:
            seen[s] = len(seen)
            stack.append(s)
    N = len(stack[0])
    bl = bonds_of(N, periodic)
    while stack:
        wd = stack.pop()
        for (i, j) in bl:
            for (pq, val) in tab.get((wd[i], wd[j]), ()):
                if pq == (wd[i], wd[j]):
                    continue
                nw = list(wd)
                nw[i], nw[j] = pq
                nw = tuple(nw)
                if nw not in seen:
                    seen[nw] = len(seen)
                    stack.append(nw)
    words = [None] * len(seen)
    for wd, k in seen.items():
        words[k] = wd
    return words, seen


def hamiltonian(words, index, h, d, periodic):
    tab = bond_table(h, d)
    N = len(words[0])
    bl = bonds_of(N, periodic)
    rows, cols, vals = [], [], []
    leak = 0
    for k, wd in enumerate(words):
        for (i, j) in bl:
            for (pq, val) in tab.get((wd[i], wd[j]), ()):
                if pq == (wd[i], wd[j]):
                    rows.append(k); cols.append(k); vals.append(val)
                else:
                    nw = list(wd)
                    nw[i], nw[j] = pq
                    kk = index.get(tuple(nw))
                    if kk is None:
                        leak += 1
                        continue
                    rows.append(kk); cols.append(k); vals.append(val)
    H = sp.coo_matrix((vals, (rows, cols)), shape=(len(words), len(words)), dtype=complex).tocsc()
    return H, leak
