"""Local bonds, generators and Cartan charges for the families searched here.

Every family returns a dict with
  d        local dimension
  hs       list of d^2 x d^2 Hermitian bond components (bond = sum coeff_i hs[i])
  charges  d x r integer array of additive Cartan charges (diagonal in the basis)
  gens     Hermitian d x d generators of the on-site symmetry (for Casimir labels)
  names    names of the components
"""
import numpy as np
from scipy import linalg
from ced import spin, swap, irrep_projectors


def _gs(cands, ip=lambda a, b: np.trace(a.conj().T @ b).real, tol=1e-9):
    out = []
    for t in cands:
        t = t.astype(complex)
        for g in out:
            t = t - ip(g, t) * g
        n = ip(t, t)
        if n > tol:
            out.append(t / np.sqrt(n))
    return out


def hermitian_basis(d):
    out = []
    for i in range(d):
        for j in range(i + 1, d):
            t = np.zeros((d, d), complex); t[i, j] = t[j, i] = 1; out.append(t)
            t = np.zeros((d, d), complex); t[i, j] = -1j; t[j, i] = 1j; out.append(t)
    for k in range(1, d):
        t = np.zeros((d, d), complex)
        t[:k, :k] = np.eye(k); t[k, k] = -k
        out.append(t)
    return out


def normalize_gens(gens, casimir_fund=None):
    """Orthonormalize w.r.t. Tr(ab) and return; caller fixes overall scale."""
    return _gs([(g + g.conj().T) / 2 for g in gens])


def sp4():
    """Sp(4) on C^4 in the spin-3/2 basis (m=3/2,1/2,-1/2,-3/2), J_{i,3-i}=(-1)^i."""
    J = np.zeros((4, 4))
    for i in range(4):
        J[i, 3 - i] = (-1) ** i
    gens = _gs([(t + J @ t.T @ J) / 2 for t in hermitian_basis(4)])
    assert len(gens) == 10
    gens = [g / np.sqrt(2) for g in gens]          # Tr(T_a T_b) = delta/2
    ps = irrep_projectors(gens)
    ranks = [p[1] for p in ps]
    P = {r: p[2] for p in ps for r in [p[1]]}
    assert sorted(ranks) == [1, 5, 10], ranks
    charges = np.array([[1, 0], [0, 1], [0, -1], [-1, 0]])
    return dict(d=4, hs=[P[5], P[10]], names=['P5', 'P10'], charges=charges, gens=gens,
                proj=P, J=J)


def spin_s(s):
    """SU(2) spin s: projectors P_J on total pair spin J, Sz charge (2m)."""
    sx, sy, sz = spin(s)
    d = int(round(2 * s + 1))
    gens = [sx, sy, sz]
    X = sum(np.kron(a, a) for a in gens).real
    vals = [(J * (J + 1) - 2 * s * (s + 1)) / 2 for J in range(d)]
    ps = []
    for j, v in enumerate(vals):
        p = np.eye(d * d)
        for k, w in enumerate(vals):
            if j != k:
                p = p @ (X - w * np.eye(d * d)) / (v - w)
        p[np.abs(p) < 1e-13] = 0
        ps.append(p)
    charges = np.round(2 * np.diag(sz)).astype(int)[:, None]
    return dict(d=d, X=X, P=ps, gens=[g.astype(complex) for g in gens], charges=charges, sz=sz)


def spin_orbital():
    """SU(2)xSU(2) on C^2 (x) C^2: components (A B, A, B), A = S.S, B = T.T."""
    sx, sy, sz = spin(0.5)
    I2 = np.eye(2)
    S = [np.kron(a, I2) for a in (sx, sy, sz)]
    T = [np.kron(I2, a) for a in (sx, sy, sz)]
    A = sum(np.kron(a, a) for a in S).real
    B = sum(np.kron(b, b) for b in T).real
    charges = np.array([[1, 1], [1, -1], [-1, 1], [-1, -1]])
    return dict(d=4, hs=[A @ B, A, B], names=['AB', 'A', 'B'], charges=charges,
                gens=[g.astype(complex) for g in S + T], S=S, T=T, A=A, B=B)


def son_vector(N):
    """SO(N) on C^N in a basis where the Cartan is diagonal.
    Components (P, K) with P the swap and K = N * (singlet projector)."""
    # real antisymmetric generators L_ab in the real basis
    gens_real = []
    for a in range(N):
        for b in range(a + 1, N):
            t = np.zeros((N, N), complex)
            t[a, b] = -1j; t[b, a] = 1j
            gens_real.append(t)
    # change of basis to weight basis: pairs (e_{2j} +- i e_{2j+1})/sqrt2, and e_{N-1} if N odd
    U = np.zeros((N, N), complex)
    col = 0
    charges = []
    r = N // 2
    for j in range(r):
        U[2 * j, col] = 1 / np.sqrt(2); U[2 * j + 1, col] = 1j / np.sqrt(2)
        q = [0] * r; q[j] = 1; charges.append(q); col += 1
    for j in range(r):
        U[2 * j, col] = 1 / np.sqrt(2); U[2 * j + 1, col] = -1j / np.sqrt(2)
        q = [0] * r; q[j] = -1; charges.append(q); col += 1
    if N % 2:
        U[N - 1, col] = 1; charges.append([0] * r)
    gens = [U.conj().T @ g @ U for g in gens_real]
    # check Cartan diagonal: generator L_{2j,2j+1}
    P = swap(N)
    # singlet vector sum_a e_a e_a in new basis: |s> = sum_a (U^T e_a?) compute directly
    s = np.zeros(N * N, complex)
    for a in range(N):
        va = U.conj().T[:, a]      # e_a in new basis coordinates
        s += np.kron(va, va)
    K = np.outer(s, s.conj())
    K[np.abs(K) < 1e-14] = 0
    return dict(d=N, hs=[P.astype(complex), K], names=['P', 'K'], charges=np.array(charges),
                gens=gens, U=U)


def g2():
    """G2 on C^7 as the stabilizer in so(7) of the octonionic 3-form, in a weight basis."""
    trip = [(0, 1, 2), (0, 3, 4), (0, 6, 5), (1, 3, 5), (1, 4, 6), (2, 3, 6), (2, 5, 4)]
    phi = np.zeros((7, 7, 7))
    for (a, b, c) in trip:
        for (x, y, z), sgn in [((a, b, c), 1), ((b, c, a), 1), ((c, a, b), 1),
                               ((b, a, c), -1), ((a, c, b), -1), ((c, b, a), -1)]:
            phi[x, y, z] = sgn
    # basis of so(7): E_ab antisymmetric real
    basis = []
    for a in range(7):
        for b in range(a + 1, 7):
            t = np.zeros((7, 7)); t[a, b] = 1; t[b, a] = -1; basis.append(t)
    # action on 3-form: (X.phi)_{abc} = X_ad phi_dbc + X_bd phi_adc + X_cd phi_abd
    M = []
    for X in basis:
        act = (np.einsum('ad,dbc->abc', X, phi) + np.einsum('bd,adc->abc', X, phi)
               + np.einsum('cd,abd->abc', X, phi))
        M.append(act.ravel())
    M = np.array(M).T
    ns = linalg.null_space(M)
    assert ns.shape[1] == 14, ns.shape
    Xs = [sum(c * b for c, b in zip(ns[:, i], basis)) for i in range(14)]
    gens = [1j * X for X in Xs]   # Hermitian
    # Cartan: find a commuting pair. Use su(3) stabilizer of e_0 inside g2 -> its Cartan.
    # Generic element approach: pick random H in g2, its centralizer in g2 is 2-dim Cartan.
    rng = np.random.default_rng(7)
    H0 = sum(rng.normal() * g for g in gens)
    C = np.array([(H0 @ g - g @ H0).ravel() for g in gens]).T
    cen = linalg.null_space(np.vstack([C.real, C.imag]))
    assert cen.shape[1] == 2, cen.shape
    H1 = sum(c * g for c, g in zip(cen[:, 0], gens))
    H2 = sum(c * g for c, g in zip(cen[:, 1], gens))
    # simultaneous diagonalization
    w, V = linalg.eigh(H1 + np.pi * H2)
    a1 = np.real(np.diag(V.conj().T @ H1 @ V))
    a2 = np.real(np.diag(V.conj().T @ H2 @ V))
    wts = np.stack([a1, a2], 1)
    # integer lattice: express weights in basis of two nonzero short weights
    nz = [i for i in range(7) if np.linalg.norm(wts[i]) > 1e-8]
    b1 = wts[nz[0]]
    b2 = None
    for i in nz[1:]:
        if abs(np.cross(b1, wts[i])) > 1e-6:
            b2 = wts[i]; break
    Bm = np.stack([b1, b2], 1)
    coords = np.linalg.solve(Bm, wts.T).T
    assert np.allclose(coords, np.round(coords), atol=1e-8), coords
    charges = np.round(coords).astype(int)
    gens_w = [V.conj().T @ g @ V for g in gens]
    gens_w = _gs(gens_w)
    gens_w = [g / np.sqrt(2) for g in gens_w]
    ps = irrep_projectors(gens_w)
    P = {p[1]: p[2] for p in ps}
    assert sorted(P) == [1, 7, 14, 27], sorted(P)
    return dict(d=7, hs=[P[7], P[14], P[27]], names=['P7', 'P14', 'P27'], charges=charges,
                gens=gens_w, proj=P, casimirs={p[1]: p[0] for p in ps})


def uls(d):
    return dict(d=d, hs=[swap(d)], names=['P'], charges=np.eye(d, dtype=int)[:, :-1],
                gens=None)
