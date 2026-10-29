"""(2) Legal cells, code, zero-port payload splitting, first order, second-order
Gaussian structure and small ED for the 12-letter bond, from my12.py only,
plus comparison of H2 with C75's 18-letter parent (bond18.pair, local copy
in work_C80). Argument: path to a scratch copy of codex/work_C80."""
import sys, os, json, itertools, time
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu, eigsh
import my12 as M

W = sys.argv[1]
sys.path.insert(0, W)
import bond18 as parent

out = {}
t0 = time.time()

# (a) legal words ====
rec = []
for n in range(1, 8):
    ws = M.legal_words(n)
    dec = [M.decode(w) for w in ws]
    ok_round = all(M.encode_word(*x) == w for x, w in zip(dec, ws))
    expected = set((j, s, b) for b in itertools.product((0, 1), repeat=n)
                   for j in range(n + 1) for s in (['v', 'DA', 'DB', 'XA', 'XB'] if j < n else ['v', 'MA', 'MB']))
    rec.append(dict(n=n, legal=len(ws), formula=(5 * n + 3) * 2 ** n, bijective=(set(dec) == expected and len(set(dec)) == len(dec)),
                    roundtrip=ok_round, max_excited=max(sum(M.excited(l) for l in w) for w in ws)))
out['legal_words'] = rec
print(rec, flush=True)

# marker-free zero-penalty words (open strings of data letters) up to length 6
mf = {}
for L in range(2, 7):
    zs = [w for w in itertools.product(M.DATA, repeat=L) if M.nbad(w) == 0]
    ring = [w for w in zs if M.allowed(w[-1], w[0])]
    mf[L] = dict(open_zero_penalty=len(zs), ring_zero_penalty=len(ring),
                 ring_types=sorted(set(''.join(sorted(set(l[0] for l in w))) for w in ring)))
out['marker_free'] = mf


# (b) zero-port payload splitting at order 1/K ====
def vacuum_vec(n):
    v = np.sqrt(2 / (n + 2)) * np.sin(np.pi * np.arange(1, n + 2) / (n + 2))
    return v


def splitting_matrix(n):
    """C_ab = <G_a| P H Q H P |G_b> on one cell (frozen left marker E0)."""
    p = M.default_params(n)
    mv = M.rule_table(M.move_rules(p))
    v = vacuum_vec(n)
    nb = 2 ** n
    cols = []  # dict bad_word -> amplitude per bit string
    bad_counts = set()
    for b in range(nb):
        bits = tuple((b >> k) & 1 for k in range(n))
        vec = {}
        for j in range(n + 1):
            w = M.encode_word(j, 'v', bits)
            full = ('E0',) + w
            for i in range(n + 1):
                for a2, b2, val in mv.get((full[i], full[i + 1]), []):
                    y = list(full); y[i] = a2; y[i + 1] = b2
                    if y[0] != 'E0':
                        continue
                    y = tuple(y[1:])
                    k = M.nbad(y, 'E0')
                    if k:
                        bad_counts.add(k)
                        vec[y] = vec.get(y, 0) + val * v[j]
        cols.append(vec)
    keys = sorted(set().union(*[c.keys() for c in cols]))
    kx = {k: i for i, k in enumerate(keys)}
    A = np.zeros((len(keys), nb), complex)
    for b, c in enumerate(cols):
        for k, val in c.items():
            A[kx[k], b] = val
    return A.conj().T @ A, bad_counts


def walsh(diag, n):
    """coefficients of diag(E) in Z-string basis, max |coef| by weight."""
    nb = 2 ** n
    res = {}
    for S in range(nb):
        chi = np.array([(-1) ** bin(b & S).count('1') for b in range(nb)])
        c = float(np.dot(diag, chi) / nb)
        w = bin(S).count('1')
        res[w] = max(res.get(w, 0.), abs(c))
    return res


rec = []
for n in range(1, 12):
    C, bc = splitting_matrix(n)
    off = float(np.max(abs(C - np.diag(np.diag(C)))))
    dg = np.real(np.diag(C))
    wl = walsh(dg, n) if n <= 10 else None
    rec.append(dict(n=n, K_times_spread=float(np.ptp(dg)), K_times_mean_shift=float(np.mean(dg)),
                    offdiag_max=off, bad_pairs_in_intermediates=sorted(bc),
                    walsh_max_by_weight=wl))
    print(rec[-1], flush=True)
out['payload_splitting_order_1_over_K'] = rec


# ED cross-check of the 1/K coefficient on the full zero-excitation cell space
def zero_exc_cell(n, K):
    p = M.default_params(n, K=K)
    seeds = [M.encode_word(j, 'v', bits) for bits in itertools.product((0, 1), repeat=n) for j in range(n + 1)]
    bonds = [(-1, 0)] + [(i, i + 1) for i in range(n)]
    st, H, V = M.closure(seeds, p, bonds, left='E0')
    return st, H


rec = []
for n in (2, 3, 4):
    C, _ = splitting_matrix(n)
    pred = np.sort(np.real(np.diag(C)))
    for K in (1e3, 1e4, 1e5):
        st, H = zero_exc_cell(n, K)
        ev = np.linalg.eigvalsh(H.toarray())[:2 ** n]
        Ev = -2 * np.cos(np.pi / (n + 2))
        rec.append(dict(n=n, K=K, dim=len(st), K_times_spread=float(K * np.ptp(ev)),
                        predicted_K_times_spread=float(np.ptp(pred)),
                        max_dev_K_times_shift=float(np.max(abs(np.sort(K * (Ev - ev)) - np.sort(pred))))))
        print(rec[-1], flush=True)
out['payload_splitting_ED'] = rec


# (c)-(e) two cells: first order, H2, Gaussian, comparison with parent ====
def code_basis(n, states_ix, N):
    v = vacuum_vec(n)
    G = np.zeros((N, 4 ** n))
    for A in range(2 ** n):
        ba = tuple((A >> k) & 1 for k in range(n))
        for B in range(2 ** n):
            bb = tuple((B >> k) & 1 for k in range(n))
            for j in range(n + 1):
                for k in range(n + 1):
                    w = M.encode_word(j, 'v', ba) + M.encode_word(k, 'v', bb)
                    G[states_ix[w], A + 2 ** n * B] += v[j] * v[k]
    return G


def two_cell(n, p, project):
    seeds = [M.encode_word(j, 'v', a) + M.encode_word(k, 'v', b)
             for a in itertools.product((0, 1), repeat=n) for b in itertools.product((0, 1), repeat=n)
             for j in range(n + 1) for k in range(n + 1)]
    bonds, ports = M.cell_bonds(n, 2)
    st, H, V = M.closure(seeds, p, bonds, left='E0', port_bonds=ports, project=project, cap=600000)
    ix = {s: i for i, s in enumerate(st)}
    exc = np.array([(sum(M.excited(l) for l in s[:n + 1]), sum(M.excited(l) for l in s[n + 1:])) for s in st])
    return st, ix, H, V, exc


def H2_infinite(n, p):
    st, ix, H, V, exc = two_cell(n, dict(p, K=np.inf), True)
    N = len(st)
    G = code_basis(n, ix, N)
    E0 = 2 * (p['musel'] - 2 * np.cos(np.pi / (n + 2)))
    hv = H @ G
    res0 = float(np.max(abs(hv - E0 * G)))
    first = float(np.max(abs(G.conj().T @ (V @ G))))
    src = np.flatnonzero((exc[:, 0] == 1) & (exc[:, 1] == 1))
    S = (V @ G)[src]
    lost = float(np.linalg.norm((V @ G)) ** 2 - np.linalg.norm(S) ** 2)
    A = H[src][:, src] - E0 * sp.eye(len(src))
    H2 = -S.conj().T @ splu(A.tocsc()).solve(S)
    return H2, dict(dim=N, code_eigen_residual=res0, first_order_max=first, source_outside_11=lost)


def gaussian_residual(Hm, n):
    m = 2 * n
    inds = [1 << j for j in range(m)]
    h = Hm[np.ix_(inds, inds)] - Hm[0, 0] * np.eye(m)
    cs = []
    for j in range(m):
        c = np.zeros((2 ** m, 2 ** m))
        for x in range(2 ** m):
            if (x >> j) & 1:
                c[x ^ (1 << j), x] = (-1) ** bin(x & ((1 << j) - 1)).count('1')
        cs.append(c)
    ref = Hm[0, 0] * np.eye(len(Hm), dtype=complex)
    for i in range(m):
        for j in range(m):
            ref = ref + h[i, j] * cs[i].T @ cs[j]
    return float(np.linalg.norm(Hm - ref, 2) / np.linalg.norm(Hm, 2))


rec = []
for n in (1, 2, 3):
    base = dict(D=4., g=0., td=1e-2, musel=0.3)
    pm = M.default_params(n, **base)
    H20, info = H2_infinite(n, pm)
    f = {}
    for g in (2e-2, 1e-2):
        H2g, _ = H2_infinite(n, dict(pm, g=g))
        f[g] = (H2g - H20) / g ** 2
    H22 = (4 * f[1e-2] - f[2e-2]) / 3
    # parent (C75 18-letter) at the same point
    pp = parent.params(n, D=4., g=1e-2, td=1e-2, musel=0.3, K=20.)
    ref = parent.pair(pp)
    H2g, _ = H2_infinite(n, dict(pm, g=1e-2))
    r = dict(n=n, **info,
             H20_nonscalar=float(np.max(abs(H20 - H20[0, 0] * np.eye(len(H20))))),
             H22_gaussian_residual=gaussian_residual(H22, n),
             H22_norm=float(np.linalg.norm(H22, 2)),
             H22_hermiticity=float(np.max(abs(H22 - H22.conj().T))),
             parent_H2_max_diff=float(np.max(abs(H2g - ref['H2']))),
             parent_H2_norm=float(np.max(abs(ref['H2']))))
    rec.append(r); print(r, flush=True)
out['second_order_infinite_K'] = rec


# (f) ED at finite K vs quasi-degenerate second order ====
def ed_check(n, K, lams):
    p = M.default_params(n, D=4., g=.04, td=1e-2, musel=0., K=K)
    st, ix, H, V, exc = two_cell(n, p, False)
    N = len(st)
    vac = np.flatnonzero((exc[:, 0] == 0) & (exc[:, 1] == 0))
    src = np.flatnonzero((exc[:, 0] == 1) & (exc[:, 1] == 1))
    Hv = H[vac][:, vac].toarray()
    ev, U = np.linalg.eigh(Hv)
    nv = 4 ** n
    E = ev[:nv]; G = np.zeros((N, nv), complex); G[vac] = U[:, :nv]
    first = float(np.max(abs(G.conj().T @ (V @ G))))
    S = (V @ G)[src]; A = H[src][:, src]
    X = np.zeros_like(S)
    for a in range(nv):
        X[:, a] = splu((A - E[a] * sp.eye(len(src))).tocsc()).solve(S[:, a])
    B = S.conj().T @ X
    H2 = -(B + B.conj().T) / 2
    rows = []
    for lam in lams:
        Ht = (H + lam * V).tocsc()
        vals, vec = eigsh(Ht, k=nv, sigma=float(E[0]) - 0.3, which='LM', ncv=min(N - 1, max(60, 4 * nv)), tol=1e-13)
        pred = np.linalg.eigvalsh(np.diag(E) + lam ** 2 * H2)
        err = float(np.max(abs(np.sort(vals) - pred)))
        rows.append(dict(n=n, K=K, lam=lam, dim=N, zero_port_splitting=float(np.ptp(E)), first_order_max=first,
                         err=err, err_over_lam4=err / lam ** 4,
                         residual=float(np.max(np.linalg.norm(Ht @ vec - vec * vals, axis=0)))))
        print(rows[-1], flush=True)
    return rows


rec = []
for n in (1, 2):
    rec += ed_check(n, 100., (0.2, 0.1, 0.05))
try:
    rec += ed_check(3, 100., (0.1, 0.05))
except Exception as e:
    rec.append(dict(n=3, error=str(e)))
out['ed_finite_K'] = rec
out['seconds'] = time.time() - t0
json.dump(out, open('check_cells12.json', 'w'), indent=1)
