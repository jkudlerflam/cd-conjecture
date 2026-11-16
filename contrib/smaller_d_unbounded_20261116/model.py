"""Twenty material states and two vacancy colors; finite diagnostics only.

The material bond has two independent charged buses, an explicit vacuum
head, and a biased path with a terminal killing term. No term is projected
out of the literal bond. See README.md for the conditional infrared status.
"""
from itertools import product

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
from scipy.sparse.linalg import splu


MATERIAL = (
    [("L", a) for a in (0, 1)]
    + [("R", a) for a in (0, 1)]
    + [("v",)]
    + [(rail, r, s) for rail in ("d", "x") for r in (0, 1) for s in (-1, 1)]
    + [("e", a) for a in (0, 1)]
    + [("m",)]
    + [("m", r, s) for r in (0, 1) for s in (-1, 1)]
)
LETTERS = MATERIAL + [("hole", a) for a in (0, 1)]
INDEX = {x: i for i, x in enumerate(LETTERS)}
assert len(MATERIAL) == 20 and len(LETTERS) == 22


def arrows(x):
    if x[0] == "L" or x[0] == "m":
        return 0, 0
    if x[0] == "R":
        return 1, 1
    if x[0] in ("v", "d", "x"):
        return 0, 1
    if x[0] == "e":
        return 1, 0
    if x[0] == "hole":
        return x[1], x[1]
    raise ValueError(x)


def allowed(x, y):
    return arrows(x)[1] == arrows(y)[0]


def charge(x):
    if x[0] in ("L", "R", "e"):
        return 2 * x[1]
    if len(x) == 3:
        return x[2]
    return 0


def marker(x):
    return x[0] in ("e", "m")


def excited(x):
    return len(x) == 3 and x[0] in ("d", "x", "m")


def vacuum_path(n, q):
    """Positions 0,...,n, with the same killing term also at n=0."""
    h = np.diag([q * q] + [1 + q * q] * n)
    if n:
        h -= q * (np.eye(n + 1, k=1) + np.eye(n + 1, k=-1))
    return h


def parameters(n=2, **changes):
    p = dict(n=n, q=0.2, w=0.7, v=0.6, eta=0.35, D=8.0,
             a=4.0, g=0.0, K=100.0, field=0.0)
    p.update(changes)
    e, vec = la.eigh(vacuum_path(n, p["q"]))
    p["E"] = float(e[0])
    p["vacuum"] = vec[:, 0] * np.sign(vec[0, 0])
    p["Delta"] = p["a"] + p["E"]
    p["selector"] = -(n + 2) * p["E"]
    return p


def material_parts(p):
    """Kinetic, conversion (unit coefficient), port, and onsite matrices."""
    d = len(MATERIAL)
    kinetic = np.zeros((d * d, d * d), complex)
    conversion = np.zeros_like(kinetic)
    port = np.zeros_like(kinetic)
    onsite = np.zeros(d)

    def term(mat, source, target, amplitude, positive=False):
        i = INDEX[source[0]] * d + INDEX[source[1]]
        j = INDEX[target[0]] * d + INDEX[target[1]]
        mat[j, i] += amplitude
        mat[i, j] += np.conj(amplitude)
        if positive:
            mat[i, i] += p["q"] ** 2
            mat[j, j] += 1.0

    for bit in (0, 1):
        term(kinetic, (("v",), ("R", bit)),
             (("L", bit), ("v",)), -p["q"], positive=True)
        term(kinetic, (("v",), ("e", bit)),
             (("L", bit), ("m",)), -p["q"], positive=True)
        for r in (0, 1):
            w = p["w"] * (p["v"] if r else 1)
            for s in (-1, 1):
                z = (-1) ** bit if s == -1 else 1
                phase = np.exp(1j * np.pi / (2 * p["n"])) if s == -1 else 1
                term(kinetic, (("d", r, s), ("R", bit)),
                     (("L", bit), ("d", r, s)), -w * z)
                term(kinetic, (("x", r, s), ("R", bit)),
                     (("L", bit), ("x", r, s)), -phase * z)
                term(kinetic, (("x", r, s), ("e", bit)),
                     (("L", bit), ("m", r, s)), -phase * z)
    for m in [x for x in MATERIAL if marker(x)]:
        for r in (0, 1):
            for s in (-1, 1):
                term(kinetic, (m, ("d", r, s)), (m, ("x", r, s)), -1)
    for r in (0, 1):
        for kind in ("R", "e"):
            term(conversion, (("d", r, -1), (kind, 1)),
                 (("d", r, 1), (kind, 0)), 1)
        weight = np.sqrt(p["eta"]) if r else 1
        # In the tree gauge the two source coefficients are (a,b)=(1,-1)
        # and (1,i), respectively. Their hopping phases differ by i.
        phase = -1 if r else -1j
        term(port, (("m",), ("v",)),
             (("m", r, -1), ("d", r, 1)), weight * phase)
        term(port, (("m",), ("v",)),
             (("m", r, 1), ("d", r, -1)), weight)
    for i, x in enumerate(MATERIAL):
        if x == ("m",):
            onsite[i] += p["q"] ** 2
        if excited(x):
            onsite[i] += p["Delta"] + p["D"] * (x[0] in ("x", "m"))
        onsite[i] += p["selector"] * marker(x)
        onsite[i] += p["field"] * (x == ("L", 1))
    return kinetic, conversion, port, onsite


def material_bond(p, lam=0.0):
    k, c, v, o = material_parts(p)
    d = len(MATERIAL)
    diagonal = [(o[i] + o[j]) / 2 + p["K"] * (not allowed(x, y))
                for i, x in enumerate(MATERIAL) for j, y in enumerate(MATERIAL)]
    return k + p["g"] * c + lam * v + np.diag(diagonal)


def mobile_bond(p, lam=0.0, J=1.0, t=1.0, mu=0.0):
    """Colored vacancy lift; onsite material terms retain the contact factor."""
    d = len(LETTERS)
    h = np.zeros((d * d, d * d), complex)
    ii = [d * i + j for i in range(20) for j in range(20)]
    h[np.ix_(ii, ii)] = J * material_bond(p, lam)
    for i, x in enumerate(MATERIAL):
        a, b = arrows(x)
        src = INDEX[("hole", a)] * d + i
        dst = i * d + INDEX[("hole", b)]
        h[dst, src] = h[src, dst] = -t
    for i, x in enumerate(LETTERS):
        for j, y in enumerate(LETTERS):
            if i >= 20 or j >= 20:
                h[d * i + j, d * i + j] += J * p["K"] * (not allowed(x, y))
            h[d * i + j, d * i + j] -= mu * ((i < 20) + (j < 20)) / 2
    return h


def cell_words(n):
    """All legal cells, in a fixed marker-position sector, for every payload."""
    out, info = [], []
    for bits in product((0, 1), repeat=n):
        for j in range(n):
            for token in [("v",)] + [(r, bus, s) for r in ("d", "x")
                                     for bus in (0, 1) for s in (-1, 1)]:
                word = ([('L', a) for a in bits[:j]] + [token]
                        + [('R', a) for a in bits[j:-1]] + [('e', bits[-1])])
                out.append(tuple(word))
                info.append((bits, token, j))
        for m in [("m",)] + [("m", bus, s) for bus in (0, 1) for s in (-1, 1)]:
            out.append(tuple([('L', a) for a in bits] + [m]))
            info.append((bits, m, n))
    assert len(set(out)) == len(out) == (9 * n + 5) * 2 ** n
    return out, info


def cell_matrices(p):
    n = p['n']
    words, info = cell_words(n)
    wi = {w: i for i, w in enumerate(words)}
    kinetic, conversion, port, onsite = material_parts(p)
    d, size = 20, len(words)
    fixed = ("e", 0)
    matrices = []
    for local, diagonal in ((kinetic, onsite), (conversion, np.zeros(d))):
        columns = {i: [(j, local[j, i]) for j in np.flatnonzero(local[:, i])]
                   for i in range(d * d)}
        rows, cols, vals = [], [], []
        for col, word in enumerate(words):
            rows.append(col); cols.append(col)
            vals.append(sum(diagonal[INDEX[x]] for x in word))
            for pos in range(-1, n):
                left = fixed if pos == -1 else word[pos]
                right = word[pos + 1]
                for row, value in columns[INDEX[left] * d + INDEX[right]]:
                    changed = list(word)
                    if pos == -1:
                        assert MATERIAL[row // d] == fixed
                    else:
                        changed[pos] = MATERIAL[row // d]
                    changed[pos + 1] = MATERIAL[row % d]
                    changed = tuple(changed)
                    assert changed in wi, (word, changed, pos)
                    rows.append(wi[changed]); cols.append(col); vals.append(value)
        matrices.append(sp.coo_matrix((vals, (rows, cols)), shape=(size, size)).tocsc())
    code = np.zeros((size, 2 ** n))
    bitlist = list(product((0, 1), repeat=n))
    for row, (bits, token, position) in enumerate(info):
        if token in (("v",), ("m",)):
            code[row, bitlist.index(bits)] = p['vacuum'][position]
    return words, info, matrices[0], matrices[1], code, port


def pair_coefficient(p):
    """Literal finite-cell coefficient of lambda^2 g^2 on the full pair code."""
    words, info, H, C, G, port = cell_matrices(p)
    size = len(words)
    wi = {w: i for i, w in enumerate(words)}
    ex = np.array([excited(token) for _, token, _ in info])
    selected = np.array([i * size + j for i in range(size) for j in range(size)
                         if ex[i] and ex[j] and info[i][1][1] == info[j][1][1]])
    pairs = {int(x): i for i, x in enumerate(selected)}
    A = (sp.kron(H, sp.eye(size)) + sp.kron(sp.eye(size), H)).tocsc()
    A = A[selected][:, selected] - 2 * (p['E'] + p['selector']) * sp.eye(len(selected), format='csc')
    W = (sp.kron(C, sp.eye(size)) + sp.kron(sp.eye(size), C)).tocsc()[selected][:, selected]
    S = np.zeros((len(selected), 4 ** p['n']), complex)
    for il in np.flatnonzero(np.any(G != 0, axis=1)):
        for ir in np.flatnonzero(np.any(G != 0, axis=1)):
            left, right = words[il], words[ir]
            col = INDEX[left[-1]] * 20 + INDEX[right[0]]
            for row in np.flatnonzero(port[:, col]):
                wl = left[:-1] + (MATERIAL[row // 20],)
                wr = (MATERIAL[row % 20],) + right[1:]
                idx = pairs[wi[wl] * size + wi[wr]]
                S[idx] += port[row, col] * np.kron(G[il], G[ir])
    lu = splu(A)
    X = lu.solve(S)
    Y = W @ X
    H22 = -Y.conj().T @ lu.solve(Y)
    H20 = -S.conj().T @ X
    return H22, H20, dict(cell_dimension=size, resolvent_dimension=len(selected),
                         code_residual=float(la.norm(H @ G - (p['E'] + p['selector']) * G)))


def fock_operators(modes):
    dim = 2 ** modes
    ops = []
    for i in range(modes):
        f = np.zeros((dim, dim))
        bit = 1 << (modes - 1 - i)
        for col in range(dim):
            if col & bit:
                parity = (col >> (modes - i)).bit_count()
                f[col ^ bit, col] = (-1) ** parity
        ops.append(f)
    return ops


def gaussian_decomposition(h, modes):
    one = [1 << (modes - 1 - i) for i in range(modes)]
    constant = h[0, 0]
    M = h[np.ix_(one, one)] - constant * np.eye(modes)
    f = fock_operators(modes)
    reconstructed = constant * np.eye(2 ** modes, dtype=complex)
    for i in range(modes):
        for j in range(modes):
            reconstructed += M[i, j] * f[i].T @ f[j]
    return M, float(la.norm(h - reconstructed) / max(la.norm(h), 1e-300))


def scalar_bus(n, w, D, a=4.0):
    """Tree in the order D_0,...,D_(n-1),X_0,...,X_(n-1),M."""
    H = np.diag([a] * n + [a + D] * (n + 1))
    for j in range(n - 1):
        H[j, j + 1] = H[j + 1, j] = -w
        H[n + j, n + j + 1] = H[n + j + 1, n + j] = -1
    H[0, n] = H[n, 0] = -1
    H[2 * n - 1, 2 * n] = H[2 * n, 2 * n - 1] = -1
    return H


def scalar_exchange(n, w, D, a=4.0):
    H = scalar_bus(n, w, D, a)
    m = len(H)
    R = la.inv(np.kron(H, np.eye(m)) + np.kron(np.eye(m), H))
    x = R[:, (m - 1) * m]
    U, V = np.zeros((m * m, n)), np.zeros((m * m, n))
    for j in range(n):
        U[j * m:(j + 1) * m, j] = x[j * m:(j + 1) * m]
        V[j::m, j] = x[j::m]
    return U.T @ R @ V


def gram(n, w, a=4.0):
    """Finite rational identity for the return-filter Gram matrix (C39)."""
    A = a * np.eye(n) - w * (np.eye(n, k=1) + np.eye(n, k=-1))
    R = la.inv(np.kron(A, np.eye(n)) + np.kron(np.eye(n), A))
    z = R[:, 0]
    G = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            G[i, j] = 2 * sum(z[i * n + y] * R[i * n + y, x * n + j] * z[x * n + j]
                             for x in range(n) for y in range(n))
    return G
