"""Shared library for the cd_free_route checks (numpy and scipy only).

Bonds, ring exact diagonalization, contiguous marginals, gauge-optimized M,
and the iterated Petz periodization channel."""
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
import scipy.sparse.linalg as spla
from scipy.optimize import linprog

sx = np.array([[0, 1], [1, 0]], dtype=complex) / 2
sy = np.array([[0, -1j], [1j, 0]], dtype=complex) / 2
sz = np.array([[1, 0], [0, -1]], dtype=complex) / 2


def xxz_bond(delta, J=1.0):
    return J * (np.kron(sx, sx) + np.kron(sy, sy) + delta * np.kron(sz, sz)).real


def heis_bond(J=1.0):
    return xxz_bond(1.0, J)


def swap_bond(d):
    P = np.zeros((d * d, d * d))
    for a in range(d):
        for b in range(d):
            P[b * d + a, a * d + b] = 1.0
    return P


def bond_dim(h):
    return int(round(np.sqrt(h.shape[0])))


def ring_hamiltonian(h, L):
    """Sparse H = sum_x h_{x,x+1} on the ring (C^d)^L, site 0 most significant."""
    d = bond_dim(h)
    hs = sp.csr_matrix(h)
    N = d ** L
    H = sp.csr_matrix((N, N))
    for x in range(L - 1):
        H = H + sp.kron(sp.kron(sp.identity(d ** x), hs), sp.identity(d ** (L - x - 2)))
    # closing bond (L-1, 0): conjugate by cyclic shift
    hh = h.reshape(d, d, d, d)  # out (a,b) in (c,e) with a on site L-1, b on site 0
    rows, cols, vals = [], [], []
    idx = np.arange(N).reshape(d, d ** (L - 2), d)  # site0, middle, siteL-1
    for a in range(d):
        for b in range(d):
            for c in range(d):
                for e in range(d):
                    v = hh[a, b, c, e]
                    if abs(v) < 1e-15:
                        continue
                    # input: site L-1 = c, site 0 = e ; output: site L-1 = a, site 0 = b
                    rows.append(idx[b, :, a].ravel())
                    cols.append(idx[e, :, c].ravel())
                    vals.append(np.full(d ** (L - 2), v))
    if rows:
        H = H + sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N, N))
    return H


def ground_states(h, L, k=4, tol=1e-9):
    H = ring_hamiltonian(h, L)
    N = H.shape[0]
    if N <= 1024:
        w, v = sla.eigh(H.toarray(), driver='evr')
    else:
        w, v = spla.eigsh(H, k=k, which='SA', tol=1e-12)
        o = np.argsort(w); w, v = w[o], v[:, o]
    E0 = w[0]
    g = [i for i in range(len(w)) if w[i] - E0 < tol * max(1, abs(E0))]
    return w, v[:, g], E0


def ti_marginal(vecs, d, L, l):
    """Translation-averaged reduced density matrix of the equal mixture of the
    columns of vecs on l contiguous sites (sites 0..l-1 after averaging)."""
    rho = np.zeros((d ** l, d ** l), dtype=complex)
    n = vecs.shape[1]
    for j in range(n):
        psi = vecs[:, j].reshape([d] * L)
        for s in range(L):
            p = np.transpose(psi, [(i + s) % L for i in range(L)])
            m = p.reshape(d ** l, d ** (L - l))
            rho += m @ m.conj().T
    rho /= (n * L)
    return (rho + rho.conj().T) / 2


def psd_power(A, p, cut=1e-12):
    w, v = sla.eigh((A + A.conj().T) / 2, driver='evr')
    keep = w > cut * max(w.max(), 1e-300)
    wp = np.zeros_like(w)
    wp[keep] = w[keep] ** p
    return (v * wp) @ v.conj().T


def ptrace_first(rho, d, n):
    """Trace out the first site of an n-site density matrix."""
    r = rho.reshape(d, d ** (n - 1), d, d ** (n - 1))
    return np.einsum('aiaj->ij', r)


def ptrace_last(rho, d, n):
    r = rho.reshape(d ** (n - 1), d, d ** (n - 1), d)
    return np.einsum('iaja->ij', r)


def petz_periodize(rho_l, d, l, steps=200, tol=1e-13, cut=1e-12):
    """Iterated Petz extension with the ring's l-site marginal rho_l.

    Window tau on l sites [A, B] (B = l-1 sites). One step:
      tau' = (id_A (x) R_{B->BC})(tau),  R(X) = rho_l^{1/2} (rho_B^{-1/2} X rho_B^{-1/2} (x) 1) rho_l^{1/2},
    the bond (A, B_1) is then frozen; record its marginal and trace out A.
    Returns the list of recorded two-site marginals."""
    rB = ptrace_last(rho_l, d, l)
    rBm = psd_power(rB, -0.5, cut)
    rl = psd_power(rho_l, 0.5, cut)
    DB = d ** (l - 1)
    KB = np.kron(np.eye(d), rBm)
    Kl = np.kron(np.eye(d), rl)
    Id = np.eye(d)
    tau = rho_l.copy()
    out = []
    n = l + 1
    for it in range(steps):
        X = KB @ tau @ KB
        tp = Kl @ np.kron(X, Id) @ Kl
        r = tp.reshape(d * d, d ** (n - 2), d * d, d ** (n - 2))
        m2 = np.einsum('aibi->ab', r)
        out.append((m2 + m2.conj().T) / 2)
        tau_new = ptrace_first(tp, d, n)
        diff = np.abs(tau_new - tau).max()
        tau = (tau_new + tau_new.conj().T) / 2
        if it > 5 and diff < tol:
            break
    return out


def M_of_bond(h):
    """Gauge-optimized half-spread M = min over a, o of (lmax - lmin)/2 of h + o x 1 - 1 x o.
    Computed by a crude local minimization (Nelder-Mead over traceless o)."""
    from scipy.optimize import minimize
    d = bond_dim(h)
    basis = []
    for i in range(d):
        for j in range(d):
            E = np.zeros((d, d), dtype=complex)
            if i < j:
                E[i, j] = E[j, i] = 1
            elif i > j:
                E[i, j] = -1j; E[j, i] = 1j
            elif i < d - 1:
                E[i, i] = 1; E[i + 1, i + 1] = -1
            else:
                continue
            basis.append(E)
    I = np.eye(d)

    def f(c):
        o = sum(ci * B for ci, B in zip(c, basis))
        w = np.linalg.eigvalsh(h + np.kron(o, I) - np.kron(I, o))
        return (w[-1] - w[0]) / 2
    best = f(np.zeros(len(basis)))
    for trial in range(3):
        x0 = np.zeros(len(basis)) if trial == 0 else 0.3 * np.random.default_rng(trial).standard_normal(len(basis))
        r = minimize(f, x0, method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-12, 'maxiter': 20000})
        best = min(best, r.fun)
    return best


def petz_defect(h, vecs, L, l, cut=1e-12):
    """One-step Petz periodization. Returns e_trial = Tr[h sigma_2*] where
    sigma_2* is the (A, B1) marginal of (id_A (x) Tr_C R_{B->BC})(rho_l)."""
    d = bond_dim(h)
    rl = ti_marginal(vecs, d, L, l)
    out = petz_periodize(rl, d, l, steps=1, cut=cut)
    return np.trace(h @ out[0]).real


def jw_ops(n):
    """Annihilation operators c_j on n qubits (JW, qubit 0 most significant)."""
    Zm = np.diag([1., -1.]); sm = np.array([[0., 1.], [0., 0.]]); I2 = np.eye(2)
    ops = []
    for j in range(n):
        m = np.array([[1.]])
        for i in range(n):
            m = np.kron(m, Zm if i < j else (sm if i == j else I2))
        ops.append(m)
    return ops


def c106_bond(t, eps, a):
    """d = 4 Jordan-Wigner bond of the C106 two-band chain
    H = sum_x f_x^dag(-i t sz) f_{x+1} + h.c. + f_x^dag(eps sx - a eps) f_x,
    on-site term split equally between the two bonds."""
    c = jw_ops(4)  # modes: (x,1),(x,2),(x+1,1),(x+1,2)
    szm = np.diag([1., -1.]); sxm = np.array([[0., 1.], [1., 0.]])
    T = -1j * t * szm
    h = np.zeros((16, 16), dtype=complex)
    for i in range(2):
        for j in range(2):
            h += T[i, j] * c[i].conj().T @ c[2 + j]
    h = h + h.conj().T
    O = eps * sxm - a * eps * np.eye(2)
    for base in (0, 2):
        for i in range(2):
            for j in range(2):
                h += 0.5 * O[i, j] * c[base + i].conj().T @ c[base + j]
    return (h + h.conj().T) / 2


def ladder_bond(eps, J):
    """d = 4 fast/slow ladder: XX on qubit a, eps XX on qubit b, on-site J Sz_a Sz_b split over the two bonds."""
    I2 = np.eye(2)
    A_ = lambda o: np.kron(o, I2)
    B_ = lambda o: np.kron(I2, o)
    h = np.zeros((16, 16), dtype=complex)
    for o in (sx, sy):
        h += np.kron(A_(o), A_(o)) + eps * np.kron(B_(o), B_(o))
    zz = A_(sz) @ B_(sz)
    h += 0.5 * J * (np.kron(zz, np.eye(4)) + np.kron(np.eye(4), zz))
    return (h + h.conj().T) / 2
