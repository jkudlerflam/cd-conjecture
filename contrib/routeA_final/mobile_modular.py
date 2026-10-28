"""Finite-field pair-product certificate for C32's mobile 28-letter core.

This ports Astra C37's certificate (codex/work_C37/modular.py and
certify_pairs.py, frozen 23-letter core) to the mobile core. The algebra is
C37's, line for line: one conversion source per species pair, a
two-carrier inverse, Jordan-Wigner signs, port coefficients -iP, i, -i, and
conjugation implemented as i -> -i. Only the carrier graph changes:

    vertices 0 = A, 1..n = d_0..d_{n-1}, n+1..2n = x_0..x_{n-1}, 2n+1 = B,
    positive edge matrices (each edge adds |c| to both endpoint diagonals),
    on-site m, plus D on A, x and B,
    neutral loop d_{n-1} -> B with amplitude -t_d e^{i(psi_0 - n phi)},
    psi_0 = pi, n phi = pi/2, so the gauged neutral loop is -i t_d,
    charged loop -t_d (psi_c = 0),
    left cell's source at B, right cell's source at A,
    conversions at data site j = carrier vertex j + 1.

The same contraction runs over complex floats (backend 'float') and over
F_p with i -> sqrt(-1) mod p (backend 'mod'). The float run is compared with
mobile_leading.py, which was validated against a literal extraction in C32's
complete pair space; the mod run is the certificate. A third check reruns
the machinery on C37's frozen graph and compares the determinant polynomial
with C37's saved pairs_n*.json.

What a pass at n proves (C37 PROOFS sections 4-5 and RESULT, transferred):
the leading quadratic symbol h(z) = T z + M + e + T^dag / z at the rational
point (m, D, t_d, e) = (1, 1000, 10^-2000, 10^-1000), and at every point with
the same reduction mod p (D -> D + p j, exponents shifted by multiples of
p - 1), has a squarefree pair-product polynomial: all 2n roots are distinct and
no two unordered pairs of roots have the same product. Combined with C32's
crossing theorem at D >= D_n (PROOFS section 6), which places all 2n roots on
the unit circle when loop and field are small enough, this gives 2n simple
Fermi points, distinct sea centres modulo pi and nonzero quartic margins
(pair transfer and umklapp included) for the mobile core at that n.

numpy and scipy only, plus integer arithmetic. Usage:
    python3 -B mobile_modular.py check          # float cross-checks
    python3 -B mobile_modular.py 1 2 ... 23     # certificates
Output: mobile_pairs_n*.json, one line per n on stdout.
"""
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

P = 65521
HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------- arithmetic
# (SparseLU, dot, det, poly, gcd, pair_polynomial copied from C37 with no change
#  of algorithm; see codex/work_C37/modular.py and certify_pairs.py)

def dot(a, b, p=P):
    a = np.asarray(a, dtype=np.int64)
    b = np.asarray(b, dtype=np.int64)
    assert a.shape[-1] * (p - 1) ** 2 < 2 ** 63
    return (a @ b) % p


class SparseLU:
    def __init__(self, A, p=P):
        self.p = p
        self.N = A.shape[0]
        graph = A.copy().astype(float)
        graph.data[:] = 1
        graph = graph + sp.eye(self.N, format='csc') * 10
        lu = splu(graph, permc_spec='MMD_AT_PLUS_A', diag_pivot_thresh=0)
        self.order = np.argsort(lu.perm_c)
        B = A[self.order, :][:, self.order].tocsr()
        rows = [{int(j): int(v) % p for j, v in zip(B.indices[B.indptr[i]:B.indptr[i + 1]],
                                                     B.data[B.indptr[i]:B.indptr[i + 1]]) if int(v) % p}
                for i in range(self.N)]
        cols = [set() for _ in rows]
        for i, row in enumerate(rows):
            for j in row:
                if i > j:
                    cols[j].add(i)
        self.L = [[] for _ in rows]
        self.U = []
        self.di = []
        for k in range(self.N):
            pivot = rows[k].get(k, 0)
            if not pivot:
                raise ArithmeticError(('zero pivot', k))
            di = pow(pivot, -1, p)
            self.di.append(di)
            rr = [(j, v) for j, v in rows[k].items() if j > k and v]
            self.U.append((np.array([j for j, v in rr], dtype=int), np.array([v for j, v in rr], dtype=np.int64)))
            for i in list(cols[k]):
                val = rows[i].pop(k, 0)
                if not val:
                    continue
                fac = val * di % p
                self.L[i].append((k, fac))
                for j, v in rr:
                    old = rows[i].get(j, 0)
                    new = (old - fac * v) % p
                    if new:
                        rows[i][j] = new
                        if i > j:
                            cols[j].add(i)
                    elif old:
                        del rows[i][j]
                        if i > j:
                            cols[j].discard(i)
            rows[k].clear()
            cols[k].clear()
        self.L = [(np.array([j for j, v in r], dtype=int), np.array([v for j, v in r], dtype=np.int64))
                  for r in self.L]

    def solve(self, B):
        vec = B.ndim == 1
        if vec:
            B = B[:, None]
        p = self.p
        X = np.array(B[self.order], dtype=np.int64, copy=True)
        for i, (js, vs) in enumerate(self.L):
            if len(js):
                X[i] = (X[i] - dot(vs, X[js], p)) % p
        for i in range(self.N - 1, -1, -1):
            js, vs = self.U[i]
            if len(js):
                X[i] = (X[i] - dot(vs, X[js], p)) % p
            X[i] = X[i] * self.di[i] % p
        out = np.empty_like(X)
        out[self.order] = X
        return out[:, 0] if vec else out


def det(A, p):
    A = np.array(A, dtype=np.int64, copy=True)
    n = len(A)
    v = 1
    for k in range(n):
        r = next((j for j in range(k, n) if A[j, k]), None)
        if r is None:
            return 0
        if r != k:
            A[[r, k]] = A[[k, r]]
            v = -v
        d = int(A[k, k])
        v = v * d % p
        if k + 1 < n:
            f = A[k + 1:, k] * pow(d, -1, p) % p
            A[k + 1:, k + 1:] = (A[k + 1:, k + 1:] - f[:, None] * A[k, k + 1:]) % p
    return v % p


def poly(T, Tc, M, eps, p):
    n = len(T)
    B = (M + eps * np.eye(n, dtype=np.int64)) % p
    vals = np.array([det((T * z * z + B * z + Tc) % p, p) for z in range(2 * n + 1)], dtype=np.int64)
    out = np.zeros(2 * n + 1, dtype=np.int64)
    fall = np.array([1], dtype=np.int64)
    fact = 1
    for k in range(2 * n + 1):
        out[:k + 1] = (out[:k + 1] + fall * int(vals[0]) * pow(fact, -1, p)) % p
        vals = np.diff(vals) % p
        if k < 2 * n:
            fall = np.convolve(fall, np.array([-k, 1], dtype=np.int64)) % p
            fact = fact * (k + 1) % p
    return out


def trim(a):
    a = list(map(int, a))
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def gcd(a, b, p):
    a = trim(a)
    b = trim(b)
    while any(b):
        inv = pow(b[-1], -1, p)
        while len(a) >= len(b) and any(a):
            c = a[-1] * inv % p
            shift = len(a) - len(b)
            for j, x in enumerate(b):
                a[j + shift] = (a[j + shift] - c * x) % p
            a = trim(a)
        a, b = b, a
    return [(x * pow(a[-1], -1, p)) % p for x in a]


def pair_polynomial(f, p):
    m = len(f) - 1
    N = m * (m + 1) // 2
    if not f[-1] or not f[0]:
        raise ArithmeticError('degree drop or zero root')
    c = (np.array(f, dtype=np.int64) * pow(int(f[-1]), -1, p)) % p
    s = [m]
    for k in range(1, 2 * N + 1):
        val = sum(int(c[m - j]) * s[k - j] for j in range(1, min(k - 1, m) + 1))
        if k <= m:
            val += k * int(c[m - k])
        s.append(-val % p)
    t = [0] + [((s[k] * s[k] + s[2 * k]) * pow(2, -1, p)) % p for k in range(1, N + 1)]
    q = [1]
    for k in range(1, N + 1):
        q.append(-sum(q[k - j] * t[j] for j in range(1, k + 1)) * pow(k, -1, p) % p)
    return list(reversed(q))


# ---------------------------------------------------------------- graphs

def mobile_graph(n, typ, sgn, m, D, td, ii):
    """Edges and diagonal of the mobile carrier. typ 0 neutral, 1 charged.
    sgn = -1 gives the conjugate (i -> -i) graph. Returns dense (S, S) object
    array of field elements given as Python numbers (ints for mod, complex for
    float); ii is the image of i."""
    S = 2 * n + 2
    H = [[0] * S for _ in range(S)]

    def edge(i, j, c, cc, absc):
        H[i][i] += absc
        H[j][j] += absc
        H[j][i] += c          # <j|H|i>
        H[i][j] += cc         # <i|H|j> = conj(c)

    edge(0, 1, -1, -1, 1)
    for j in range(n - 1):
        edge(1 + j, 2 + j, -1, -1, 1)
    edge(1, n + 1, -1, -1, 1)
    for j in range(n - 1):
        edge(n + 1 + j, n + 2 + j, -1, -1, 1)
    edge(2 * n, 2 * n + 1, -1, -1, 1)
    if td:
        if typ == 0:
            edge(n, 2 * n + 1, -td * sgn * ii, td * sgn * ii, td)
        else:
            edge(n, 2 * n + 1, -td, -td, td)
    for v in range(S):
        H[v][v] += m
    for v in [0] + list(range(n + 1, 2 * n + 2)):
        H[v][v] += D
    return H, 0, 2 * n + 1, list(range(1, n + 1))   # graph, right source, left source, data vertices


def frozen_graph(n, typ, sgn, a, D, td, ii):
    """C37's frozen carrier (modular.py), for the regression check."""
    s = 2 * n + 1
    H = [[0] * s for _ in range(s)]
    for i in range(n):
        H[i][i] = a
    for i in range(n, s):
        H[i][i] = a + D
    for j, k in [(j, j + 1) for j in range(n - 1)] + [(n + j, n + j + 1) for j in range(n)] + [(0, n)]:
        H[j][k] = H[k][j] = -1
    c = td * (sgn * ii if typ == 0 else 1)
    cc = td * (-sgn * ii if typ == 0 else 1)
    H[s - 1][n - 1] -= c
    H[n - 1][s - 1] -= cc
    return H, 0, s - 1, list(range(n))


# ---------------------------------------------------------------- contraction

def matrices(n, graph, params, backend, p=P):
    """C37's matrices(), generalized to a carrier graph and a backend.
    backend 'mod': entries in F_p; 'float': complex floats."""
    if backend == 'mod':
        ii = next(j for j in range(2, p) if j * j % p == p - 1)
        red = lambda x: int(x) % p  # noqa: E731
    else:
        ii = 1j
        red = complex  # noqa: E731
    G = {}
    for typ in (0, 1):
        for sgn in (1, -1):
            H, rsrc, lsrc, data = graph(n, typ, sgn, *params, ii)
            G[typ, sgn] = H
    s = len(G[0, 1])
    if backend == 'mod':
        mats = {k: sp.csc_matrix(np.array([[red(x) for x in row] for row in H], dtype=np.int64))
                for k, H in G.items()}
        I = sp.eye(s, dtype=np.int64, format='csc')
    else:
        mats = {k: sp.csc_matrix(np.array([[complex(x) for x in row] for row in H]))
                for k, H in G.items()}
        I = sp.identity(s, format='csc')

    def pair(x, y, minus=False):
        sg = -1 if minus else 1
        A = (sp.kron(mats[int(x != 0), sg], I) + sp.kron(I, mats[int(y != 0), sg])).tocsc()
        if backend == 'mod':
            A.data %= p
        return A

    def solver(A):
        if backend == 'mod':
            return SparseLU(A, p)
        return splu(A)

    def mdot(a, b):
        return dot(a, b, p) if backend == 'mod' else a @ b

    dtype = np.int64 if backend == 'mod' else complex
    b = np.zeros(s * s, dtype=dtype)
    b[lsrc * s + rsrc] = 1
    sources = {}
    for key in [(0, 0), (1, 2), (2, 1)]:
        x, y = key
        for sign in [1, -1]:
            A = pair(x, y, sign == -1)
            sol = solver(A).solve(b)
            if backend == 'mod':
                assert np.array_equal((A @ sol) % p, b)
            sources[key, sign] = sol.reshape(s, s)
    basis = {}
    groups = {}
    for x, y in [(0, 0), (1, 2), (2, 1)]:
        for cell, typ in enumerate([x, y]):
            for ns in ([1, 2] if typ == 0 else [0]):
                target = (ns, y) if cell == 0 else (x, ns)
                for j in range(n):
                    key = (x, y, cell, ns, j)
                    idx = len(groups.setdefault(target, []))
                    groups[target].append(key)
                    basis[key] = (target, idx)
    gram = {}
    for target, keys in groups.items():
        vp = np.zeros((s * s, len(keys)), dtype=dtype)
        vm = np.zeros_like(vp)
        for col, (x, y, cell, ns, j) in enumerate(keys):
            dv = data[j]
            for sign, V in [(1, vp), (-1, vm)]:
                R = sources[(x, y), sign]
                z = V[:, col].reshape(s, s)
                if cell == 0:
                    z[dv, :] = R[dv, :]
                else:
                    z[:, dv] = R[:, dv]
        A = pair(*target)
        w = solver(A).solve(vp)
        if backend == 'mod':
            assert np.array_equal((A @ w) % p, vp % p)
            gram[target] = -dot(vm.T, w, p) % p
        else:
            gram[target] = -(vm.T @ w)

    def terms(bits):
        out = []
        for x, y in [(0, 0), (1, 2), (2, 1)]:
            coeff = (-ii * ((-1) ** bin(bits & ((1 << n) - 1)).count('1'))) if x == 0 else (ii if x == 1 else -ii)
            for cell, typ in enumerate([x, y]):
                local = (bits >> (cell * n)) & ((1 << n) - 1)
                for j in range(n):
                    bit = (local >> j) & 1
                    if typ == 0:
                        ns = 1 if bit else 2
                    elif typ == 1:
                        if bit:
                            continue
                        ns = 0
                    else:
                        if not bit:
                            continue
                        ns = 0
                    sg = (-1) ** bin(local & ((1 << j) - 1)).count('1')
                    target, idx = basis[(x, y, cell, ns, j)]
                    c = coeff * sg
                    out.append(((target, bits ^ (1 << (cell * n + j))), idx, c % p if backend == 'mod' else c))
        return out

    ts = [terms(bb) for bb in [0] + [1 << j for j in range(2 * n)]]
    out = np.zeros((2 * n + 1, 2 * n + 1), dtype=dtype)
    for i, xx in enumerate(ts):
        for j, yy in enumerate(ts):
            total = 0
            for (key, k, c) in xx:
                for (key2, l, d) in yy:
                    if key == key2:
                        g = gram[key[0]][k, l]
                        total -= c * d * (int(g) if backend == 'mod' else g)
            out[i, j] = total % p if backend == 'mod' else total
    if backend == 'mod':
        h = (out[1:, 1:] - out[0, 0] * np.eye(2 * n, dtype=np.int64)) % p
        return h[:n, n:], h[n:, :n], (h[:n, :n] + h[n:, n:]) % p, ii
    h = out[1:, 1:] - out[0, 0] * np.eye(2 * n)
    return h[:n, n:], h[n:, :n], h[:n, :n] + h[n:, n:], ii


# ---------------------------------------------------------------- checks

def float_check():
    """Same contraction in complex floats, compared with mobile_leading.py
    (literal-validated) at psi_0 = pi, psi_c = 0. The contraction omits the
    visibility 1/(n+2)^2, so it is restored before comparing. The onsite
    block of this generic contraction is a difference of one-particle and
    vacuum diagonals, so in floats it loses accuracy as D grows; the small-D
    rows test the onsite algebra, the large-D rows the hopping."""
    from mobile_leading import leading
    res = []
    for n, D, td in [(1, 2.0, 0.05), (2, 2.0, 0.05), (3, 2.0, 0.02), (4, 2.0, 0.01),
                     (3, 5.0, 1e-4), (4, 10.0, 1e-7), (5, 10.0, 1e-9)]:
        T, Tc, M, _ = matrices(n, mobile_graph, (1.0, D, td), 'float')
        vis = 1.0 / (n + 2) ** 2
        Tl, Ml = leading(n, m=1.0, D=D, td=td, psi0=np.pi, psic=0.0)
        sc = np.max(abs(Tl))
        res.append(dict(n=n, D=D, td=td,
                        T_rel_err=float(np.max(abs(vis * T - Tl)) / sc),
                        Tc_vs_Tdag=float(np.max(abs(vis * Tc - Tl.conj().T)) / sc),
                        M_err_over_T=float(np.max(abs(vis * M - Ml)) / sc),
                        M_over_T=float(np.max(abs(Ml)) / sc)))
    return res


def frozen_regression(ns=(2, 3, 4)):
    """Rerun the machinery on C37's frozen graph at C37's parameters and
    compare the determinant polynomial with C37's saved certificates."""
    res = []
    for n in ns:
        p = P
        tdm = pow(pow(10, 2000, p), -1, p)
        epsm = pow(pow(10, 1000, p), -1, p)
        T, Tc, M, _ = matrices(n, frozen_graph, (4, 1000, tdm), 'mod')
        f = poly(T, Tc, M, epsm, p).tolist()
        path = os.path.join(HERE, '..', '..', 'codex', 'work_C37', f'pairs_n{n}.json')
        ref = json.load(open(path))['f']
        res.append(dict(n=n, matches_C37=bool(f == ref)))
    return res


def certify(n, m=1, D=1000, td_exp=2000, eps_exp=1000, p=P):
    start = time.monotonic()
    tdm = pow(pow(10, td_exp, p), -1, p)
    epsm = pow(pow(10, eps_exp, p), -1, p)
    T, Tc, M, ii = matrices(n, mobile_graph, (m, D, tdm), 'mod', p)
    f = poly(T, Tc, M, epsm, p)
    checks = []
    for z in [2 * n + 1, 2 * n + 3, 2 * n + 7]:
        value = sum(int(c) * pow(z, j, p) for j, c in enumerate(f)) % p
        direct = det((T * z * z + (M + epsm * np.eye(n, dtype=np.int64)) * z + Tc) % p, p)
        assert value == direct
        checks.append([z, value])
    q = pair_polynomial(f, p)
    g = gcd(q, [(j * q[j]) % p for j in range(1, len(q))], p)
    r = dict(core='C32 mobile 28-letter', n=n, m=m, D=D, td_power=-td_exp, field_power=-eps_exp,
             psi0='pi', psic=0, express_phase='pi/(2n)', port='(1, i, -i)', p=p, image_of_i=ii,
             carrier_dimension=2 * n + 2, pair_resolvent_dimension=(2 * n + 2) ** 2,
             f=f.tolist(), pair_product_polynomial=q, gcd=g, extra_determinant_checks=checks,
             status='no_pair_product_collision' if len(g) == 1 else 'inconclusive_modular_collision',
             seconds=round(time.monotonic() - start, 2))
    with open(os.path.join(HERE, f'mobile_pairs_n{n}.json'), 'w') as fh:
        json.dump(r, fh, indent=1)
        fh.write('\n')
    print(json.dumps({k: v for k, v in r.items() if k not in ('f', 'pair_product_polynomial', 'gcd')}),
          flush=True)
    return r


def certify_frozen(n, a=4, D=1000, td_exp=2000, eps_exp=1000, p=P):
    """C37's frozen 23-letter certificate at C37's parameters, for n beyond
    C37's range (Form V's internal chain). Same machinery as certify()."""
    start = time.monotonic()
    tdm = pow(pow(10, td_exp, p), -1, p)
    epsm = pow(pow(10, eps_exp, p), -1, p)
    T, Tc, M, ii = matrices(n, frozen_graph, (a, D, tdm), 'mod', p)
    f = poly(T, Tc, M, epsm, p)
    q = pair_polynomial(f, p)
    g = gcd(q, [(j * q[j]) % p for j in range(1, len(q))], p)
    r = dict(core='frozen 23-letter (C37 graph and parameters)', n=n, a=a, D=D, td_power=-td_exp,
             field_power=-eps_exp, p=p, image_of_i=ii, carrier_dimension=2 * n + 1,
             pair_resolvent_dimension=(2 * n + 1) ** 2, f=f.tolist(), pair_product_polynomial=q, gcd=g,
             status='no_pair_product_collision' if len(g) == 1 else 'inconclusive_modular_collision',
             seconds=round(time.monotonic() - start, 2))
    with open(os.path.join(HERE, f'frozen_pairs_n{n}.json'), 'w') as fh:
        json.dump(r, fh, indent=1)
        fh.write('\n')
    print(json.dumps({k: v for k, v in r.items() if k not in ('f', 'pair_product_polynomial', 'gcd')}),
          flush=True)
    return r


def controls(n, p=P):
    """The test must detect the known collisions: no loop (tree, centres
    coincide) and no field."""
    out = {}
    for name, td_zero, eps_zero in [('no_loop', True, False), ('no_field', False, True)]:
        tdm = 0 if td_zero else pow(pow(10, 2000, p), -1, p)
        epsm = 0 if eps_zero else pow(pow(10, 1000, p), -1, p)
        T, Tc, M, _ = matrices(n, mobile_graph, (1, 1000, tdm), 'mod', p)
        f = poly(T, Tc, M, epsm, p)
        try:
            q = pair_polynomial(f, p)
            g = gcd(q, [(j * q[j]) % p for j in range(1, len(q))], p)
            out[name] = 'collision detected' if len(g) > 1 else 'no collision'
        except ArithmeticError as e:
            out[name] = 'degenerate: ' + str(e)
    return out


if __name__ == '__main__':
    if sys.argv[1:] == ['check']:
        print(json.dumps({'float_vs_structured': float_check(),
                          'frozen_regression': frozen_regression(),
                          'controls': {n: controls(n) for n in (2, 3, 4)}}, indent=1))
    elif sys.argv[1:2] == ['frozen']:
        for n in map(int, sys.argv[2:]):
            certify_frozen(n)
    else:
        for n in map(int, sys.argv[1:]):
            certify(n)
