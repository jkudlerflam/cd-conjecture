"""Bond operators of the benchmark chains, exact diagonalization on periodic rings, and the
lattice quantities entering the bounds: ||h|| (half the spread of spec h), the energy-current
norm ||j|| with j = i[h_12, h_23], and the Casimir coefficient from E_0(L).

Matrix-free: H acts on a (d,)*L tensor bond by bond, so no Hamiltonian is stored.
"""
import numpy as np
import scipy.sparse.linalg as sla

w3 = np.exp(2j * np.pi / 3)
X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1.0 + 0j, -1])
I2 = np.eye(2)


def spin_ops(S):
    m = np.arange(S, -S - 1, -1)
    Sz = np.diag(m).astype(complex)
    Sp = np.zeros((len(m), len(m)), complex)
    for i in range(1, len(m)):
        Sp[i - 1, i] = np.sqrt(S * (S + 1) - m[i] * (m[i] + 1))
    Sx = (Sp + Sp.conj().T) / 2; Sy = (Sp - Sp.conj().T) / 2j
    return Sx, Sy, Sz


def swap(d):
    P = np.zeros((d * d, d * d))
    for a in range(d):
        for b in range(d):
            P[b * d + a, a * d + b] = 1
    return P.astype(complex)


def models():
    """name: (h, d, c, v_exact or None)."""
    M = {}
    M["Ising"] = (-np.kron(Z, Z) - 0.5 * (np.kron(X, I2) + np.kron(I2, X)), 2, 0.5, 2.0)
    M["XX"] = (0.5 * (np.kron(X, X) + np.kron(Y, Y)), 2, 1.0, 2.0)
    Sx, Sy, Sz = spin_ops(0.5)
    g = np.pi / 3
    M["XXZ_D0.5"] = (np.kron(Sx, Sx) + np.kron(Sy, Sy) + 0.5 * np.kron(Sz, Sz), 2, 1.0, (np.pi / 2) * np.sin(g) / g)
    sig = np.diag([1, w3, w3 ** 2]); tau = np.roll(np.eye(3), 1, axis=0).astype(complex)
    I3 = np.eye(3)
    hp = -(np.kron(sig, sig.conj().T) + np.kron(sig.conj().T, sig)) - 0.5 * (np.kron(tau + tau.conj().T, I3) + np.kron(I3, tau + tau.conj().T))
    M["Potts3"] = (hp, 3, 0.8, 3 * np.sqrt(3) / 2)
    S1 = spin_ops(1.0)
    SS = sum(np.kron(a, a) for a in S1)
    M["TB_spin1"] = (SS - SS @ SS, 3, 1.5, None)
    for d in (2, 3, 4, 5):
        M[f"ULS_d{d}"] = (swap(d), d, d - 1.0, 2 * np.pi / d)
    return M


def half_spread(h):
    e = np.linalg.eigvalsh(h)
    return (e[-1] - e[0]) / 2


def current_norm(h, d):
    """||i[h_12, h_23]|| on three sites."""
    Id = np.eye(d)
    a = np.kron(h, Id); b = np.kron(Id, h)
    return np.max(np.abs(np.linalg.eigvalsh(1j * (a @ b - b @ a))))


def apply_H(psi, h4, L, d):
    out = np.zeros_like(psi)
    for i in range(L):
        j = (i + 1) % L
        t = np.tensordot(h4, psi, axes=([2, 3], [i, j]))  # new axes (i, j) in front
        out += np.moveaxis(t, [0, 1], [i, j])
    return out


def ground_energy(h, d, L, k=1):
    h4 = h.reshape(d, d, d, d)
    shape = (d,) * L
    op = sla.LinearOperator((d ** L, d ** L), matvec=lambda v: apply_H(v.reshape(shape), h4, L, d).ravel(), dtype=complex)
    ev = sla.eigsh(op, k=k, which="SA", tol=1e-10, maxiter=5000)[0]
    return np.sort(ev.real)


def full_spectrum(h, d, L):
    """Dense spectrum (only for small d**L)."""
    D = d ** L
    h4 = h.reshape(d, d, d, d)
    H = np.zeros((D, D), complex)
    eye = np.eye(D, dtype=complex)
    for col in range(D):
        H[:, col] = apply_H(eye[col].reshape((d,) * L), h4, L, d).ravel()
    return np.linalg.eigvalsh(H)
