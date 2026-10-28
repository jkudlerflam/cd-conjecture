#!/usr/bin/env python3
"""Route-A cells: one big site carrying an exact identical-cell Hamiltonian H0 with a right
port and a left port; neighbouring cells are coupled by  lam * sum_a FR_a(j) FL_a(j+1).

load_c8()   the 48-dimensional N = 3, q = 2 cell of codex/work_C8/two_port_N3.npz
            (physical order R, L, Raux, Laux, logical; second-order bond -7177/8192 + 27/8192 P_3).
build_twoport(K)  an identical two-qubit-port cell (C8 geometry, C9 section 2 freedom) whose
            second-order logical bond is a0 + tau^2 K for a prescribed Hermitian K.

Every cell is returned as a dict with H0 (D x D), G (D x N), FR, FL (3 x D x D), all in the
eigenbasis of H0 (so H0 is diagonal) and with G's columns spanning the kernel of H0.
numpy/scipy only.  Written by Claude (Opus 5.5) for contrib/routeA_numerics, 2026-10-27.
"""
import os
import numpy as np
from scipy import linalg

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], complex)


def swap(N):
    return np.eye(N * N).reshape(N, N, N, N).transpose(1, 0, 2, 3).reshape(N * N, N * N)


def to_eigenbasis(H0, G, FR, FL):
    e, V = linalg.eigh(H0)
    e[np.abs(e) < 1e-12] = 0.0
    N = G.shape[1]
    # make the first N eigenvectors exactly span ran(G), with G itself as the basis there
    Gp = V.conj().T @ G
    assert np.linalg.norm(Gp[N:]) < 1e-10, 'kernel of H0 is not ran G'
    V = V.copy()
    V[:, :N] = G                                   # same span; G is an isometry
    rot = lambda O: V.conj().T @ O @ V
    return dict(H0=np.diag(e).astype(complex), e=e, V=V, N=N, D=len(e),
                G=np.eye(len(e), N, dtype=complex),
                FR=np.array([rot(a) for a in FR]), FL=np.array([rot(a) for a in FL]))


def load_c8():
    d = np.load(os.path.join(ROOT, 'codex', 'work_C8', 'two_port_N3.npz'))
    c = to_eigenbasis(d['H0'], d['G'], d['FR'], d['FL'])
    c['H2_c8'] = d['H2']
    c['name'] = 'c8_N3'
    return c


def hermitian_basis(N):
    """Orthonormal (tr A B = delta) Hermitian basis of N x N matrices, identity/sqrt(N) first."""
    out = [np.eye(N, dtype=complex) / np.sqrt(N)]
    for i in range(N):
        for j in range(i + 1, N):
            a = np.zeros((N, N), complex); a[i, j] = a[j, i] = 1 / np.sqrt(2); out.append(a)
            a = np.zeros((N, N), complex); a[i, j] = -1j / np.sqrt(2); a[j, i] = 1j / np.sqrt(2); out.append(a)
    for k in range(1, N):
        a = np.zeros((N, N), complex)
        a[np.arange(k), np.arange(k)] = 1
        a[k, k] = -k
        out.append(a / np.linalg.norm(a))
    return out


def second_order(cell):
    """Direct second-order bond from every excited eigenvector (independent of the weight formula)."""
    e, G, FR, FL, N = cell['e'], cell['G'], cell['FR'], cell['FL'], cell['N']
    ex = np.arange(N, cell['D'])
    tr = [a[ex][:, :N] for a in FR]          # <exc| F_R |g>
    tl = [a[ex][:, :N] for a in FL]
    amp = sum(np.kron(a, b) for a, b in zip(tr, tl))           # (exc x exc, N x N)
    den = (e[ex][:, None] + e[ex][None, :]).reshape(-1)
    return -amp.conj().T @ (amp / den[:, None])


def build_twoport(K_target, N, tau=None, completion='qr', seed=0, extra=None):
    """Identical two-qubit-port cell whose second-order logical bond is
         H2 = a0 * 1 + tau^2 * K_conn,     K_conn = K_target minus its trace part,
    for any Hermitian K_target on C^N (x) C^N with N^2 <= 8.  This is the C8 section-4.2 geometry
    (R, L, Raux, Laux, logical; G = Bell pairs (x) logical) with G1's syndrome trick, so only the
    two excited energies 1 and 3 are needed (C9 section 2 allows any logical bond; its own k = m + 2
    energy recipe is numerically ill conditioned, the Cauchy matrix has eigenvalues down to 1e-17).
        M_R,e = 1/2 + tau X_e sum_mu A_p(mu) (x) T_mu,
        M_L,e = 1/2 - tau X_e sum_mu A_p(mu)^T (x) (sum_nu J_mu,nu T_nu),     X = sqrt6 (1, -1).
    T_mu runs over an orthonormal Hermitian basis of C^N including T_0 = 1/sqrt N, A_p over eight
    orthonormal traceless Hermitian 3 x 3 syndrome matrices.  K_target = sum J_mu,nu T_mu (x) T_nu.
    extra: optional N x N Hermitian on-cell logical term, added as G extra G^dagger to H0."""
    r = 3
    T = hermitian_basis(N)
    A = hermitian_basis(r)[1:]
    J = np.array([[np.real(np.trace(np.kron(a, b) @ K_target)) for b in T] for a in T])
    J[0, 0] = 0.0
    if np.abs(J[0]).max() < 1e-14:          # no (1 (x) T) terms: T_0 not needed on the right port
        T, J = T[1:], J[1:, 1:]
    m = len(T)
    assert m <= len(A), 'N^2 > 8 with on-site terms of the right cell'
    X = np.sqrt(6.0) * np.array([1.0, -1.0])
    E = np.array([1.0, 3.0])
    w = np.array([0.5, 0.5])
    DR = sum(np.kron(A[mu], T[mu]) for mu in range(m))
    DL = sum(np.kron(A[mu].T, sum(J[mu, nu] * T[nu] for nu in range(m))) for mu in range(m))
    if tau is None:
        nrm = max(np.linalg.norm(DR, 2), np.linalg.norm(DL, 2))
        tau = 0.75 * 0.5 / (np.sqrt(6.0) * nrm)          # every weight >= 1/8
    MR = [w[e] * np.eye(r * N) + tau * X[e] * DR for e in range(2)]
    ML = [w[e] * np.eye(r * N) - tau * X[e] * DL for e in range(2)]
    D = 16 * N
    G = np.zeros((D, N), complex)
    for a in range(2):
        for b in range(2):
            for i in range(N):
                G[((((a * 2 + b) * 2 + a) * 2 + b) * N + i), i] = 0.5
    I2 = np.eye(2)
    tens = lambda xs: np.kron(np.kron(np.kron(np.kron(xs[0], xs[1]), xs[2]), xs[3]), xs[4])
    FR = [tens([s, I2, I2, I2, np.eye(N)]) for s in PAULI]
    FL = [tens([I2, s, I2, I2, np.eye(N)]) for s in PAULI]
    S = np.concatenate([a @ G for a in FR + FL], axis=1)
    gs = np.concatenate([G, S], axis=1)
    assert np.allclose(gs.conj().T @ gs, np.eye(gs.shape[1]), atol=1e-12)
    rng = np.random.default_rng(seed)
    if completion == 'qr':
        frame, _ = np.linalg.qr(gs, mode='complete')
    else:
        Xr = rng.normal(size=(D, D)) + 1j * rng.normal(size=(D, D)); Xr[:, :gs.shape[1]] = gs
        frame, _ = np.linalg.qr(Xr)
    frame[:, :gs.shape[1]] = gs
    FQ = frame[:, N:]
    s6 = 6 * N
    W = np.zeros((D - N, s6), complex)
    for e in range(2):
        M = np.zeros((s6, s6), complex); M[:r * N, :r * N] = MR[e]; M[r * N:, r * N:] = ML[e]
        ev, u = linalg.eigh(M)
        assert ev.min() > 0
        W[e * s6:(e + 1) * s6] = (u * np.sqrt(ev)) @ u.conj().T
    if completion == 'qr':
        U, _ = np.linalg.qr(W, mode='complete')
    else:
        Xr = rng.normal(size=(D - N, D - N)) + 1j * rng.normal(size=(D - N, D - N)); Xr[:, :s6] = W
        U, _ = np.linalg.qr(Xr)
    U[:, :s6] = W
    en = np.full(D - N, 3.0); en[:s6] = 1.0
    exc = FQ @ U.conj().T
    H0 = (exc * en) @ exc.conj().T
    H0 = (H0 + H0.conj().T) / 2
    a0 = -sum(1.0 / (E[e] + E[f]) * w[e] * w[f] for e in range(2) for f in range(2)) * r
    Kc = K_target - np.trace(K_target) / N ** 2 * np.eye(N * N)
    pred = a0 * np.eye(N * N) + tau ** 2 * Kc
    c = to_eigenbasis(H0, G, FR, FL)
    c.update(name='twoport', pred2=pred, tau=tau, K_target=K_target)
    if extra is not None:
        c['H0'] = c['H0'].copy(); c['H0'][:N, :N] += extra
        c['extra'] = extra
    return c


if __name__ == '__main__':
    c = load_c8()
    N = c['N']
    print('C8 cell: D =', c['D'], ' spectrum multiplicities',
          {float(x): int(n) for x, n in zip(*np.unique(np.round(c['e'], 9), return_counts=True))})
    h2 = second_order(c)
    print('|H2 - (-7177/8192 + 27/8192 P)| =', np.linalg.norm(h2 - (-7177 / 8192 * np.eye(9) + 27 / 8192 * swap(3))))
    print('first-order compressions:', max(np.linalg.norm(a[:N, :N]) for a in list(c['FR']) + list(c['FL'])))
