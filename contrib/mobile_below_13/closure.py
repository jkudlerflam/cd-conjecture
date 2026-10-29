"""Ring closures for small exact diagonalization (numpy/scipy only)."""
import numpy as np
import scipy.sparse as sp
import bond12 as b

d = 12
VI = b.IX[b.V]
ALW = np.array([[b.allowed(x, y) for y in b.LET] for x in b.LET])


def nbad(w, open_chain=False):
    L = len(w)
    return sum(not ALW[w[i], w[(i + 1) % L]] for i in range(L - open_chain))


def nbad_open(w):
    return nbad(w, True)


def matbad(w):
    m = [x for x in w if x != VI]
    return sum(not ALW[m[i], m[(i + 1) % len(m)]] for i in range(len(m))) if m else 0


def nbad_material_pairs(w):
    L = len(w)
    return sum((w[i] != VI and w[(i + 1) % L] != VI) and not ALW[w[i], w[(i + 1) % L]] for i in range(L))


def closure(seed, bonds, legal_only=False, cap=200000, ill_depth=None, penfun=None):
    """seed: tuple of letters (tuples) or indices.  bonds: list of L bond
    matrices (144x144, bond i acts on sites i,i+1 mod L) WITHOUT penalty.
    legal_only: restrict to penfun==0 states (K = infinity projection).
    ill_depth: keep illegal (penfun>0) states only within this many steps of
    the penfun==0 set (variational truncation).  Returns the Hermitian
    restriction H (csc), penalty-count vector, squeezed-illegal vector, states."""
    seed = tuple(b.IX[x] if not isinstance(x, (int, np.integer)) else int(x) for x in seed)
    L = len(seed)
    penfun = penfun or nbad
    trans = [None if h is None else {c: [(int(r), h[r, c]) for r in np.flatnonzero(h[:, c])]
                                     for c in range(d * d)} for h in bonds]

    def moves(w):
        for i in range(L):
            if trans[i] is None:  # open boundary
                continue
            j = (i + 1) % L
            for r, v in trans[i][w[i] * d + w[j]]:
                z = list(w); z[i], z[j] = divmod(r, d)
                yield tuple(z), v
    # phase 1: reachable set with depth bookkeeping (Dijkstra-like on 0/1 weights)
    from collections import deque
    pen0 = penfun(seed)
    depth = {seed: 0 if pen0 == 0 else 1}
    dq = deque([seed])
    while dq:
        w = dq.popleft()
        dw = depth[w]
        if ill_depth is not None and dw >= ill_depth:
            continue
        for z, _ in moves(w):
            pz = penfun(z)
            if legal_only and pz > 0:
                continue
            dz = 0 if pz == 0 else dw + 1
            if z not in depth or dz < depth[z]:
                depth[z] = dz
                if dz == 0:
                    dq.appendleft(z)
                else:
                    dq.append(z)
        if len(depth) > cap:
            raise RuntimeError('closure cap')
    states = sorted(depth, key=lambda x: (depth[x], x))
    states.remove(seed); states = [seed] + states
    index = {w: i for i, w in enumerate(states)}
    rows, cols, vals = [], [], []
    for c, w in enumerate(states):
        for z, v in moves(w):
            if z in index:
                rows.append(index[z]); cols.append(c); vals.append(v)
    n = len(states)
    H = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsc()
    assert (abs(H - H.getH()).max() if n > 1 else 0) < 1e-12
    pen = np.array([penfun(w) for w in states], float)
    mb = np.array([matbad(w) for w in states], float)
    return H, pen, mb, states


def squeeze(w):
    return tuple(x for x in w if x != VI)
