#!/usr/bin/env python3
"""Two-body screen for nearest-neighbour qubit bonds (contrib/d2_two_fluid).

Particle language: |0> = spin down (empty), |1> = spin up (particle); left site of the
bond is the left kron factor.  The bond is

  h = -tau |10><01| - conj(tau) |01><10|                    hopping, tau = t e^{i phi}
      + V |11><11| - mu (n x 1 + 1 x n)/2                    interaction, chemical potential
      + gL |11><01| + gR |11><10| + h.c.                     conditional creation (charge 1)
      + Delta |11><00| + h.c.                                pairing (charge 2)
      + fL |10><00| + fR |01><00| + h.c.                     unconditional creation (charge 1)

Every Hermitian 4x4 bond can be written this way up to a constant and the telescoping
term o x 1 - 1 x o (`to_pauli` returns the (J, B) of specs/_d2_background.md).

With Delta = 0 and fL + fR = 0 the all-down state is an exact eigenstate (the "exact-vacuum"
family, gauge fL = fR = 0).  Its lowest excitations above the vacuum are then obtained, in
the Tamm-Dancoff truncation N <= 2, from a tridiagonal problem at each total momentum K:

  basis |S>  = one particle with momentum K,
        |r>  = two particles at separation r >= 1, total momentum K,
  <S|H|S> = -2t cos(K+phi) - mu,      <r|H|r> = -2mu + V delta_{r1},
  <1|H|S> = A(K) = gL e^{iK/2} + gR e^{-iK/2},   <r+1|H|r> = -2t cos(K/2+phi).

With gL = gR = 0 this is the exact two-up-spin problem of the U(1) chain.
"""
import numpy as np
from scipy.linalg import eigvalsh_tridiagonal

PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], complex)
I2 = np.eye(2)


def ket(a, b):
    v = np.zeros(4, complex); v[2 * a + b] = 1; return v


def bond(t=1.0, phi=0.0, V=0.0, mu=0.0, gL=0.0, gR=0.0, Delta=0.0, fL=0.0, fR=0.0):
    tau = t * np.exp(1j * phi)
    P = lambda a, b, c, d: np.outer(ket(a, b), ket(c, d).conj())
    h = -tau * P(1, 0, 0, 1) + V * P(1, 1, 1, 1)
    h += gL * P(1, 1, 0, 1) + gR * P(1, 1, 1, 0) + Delta * P(1, 1, 0, 0)
    h += fL * P(1, 0, 0, 0) + fR * P(0, 1, 0, 0)
    h = h + h.conj().T - (h.diagonal().real * np.eye(4))   # diagonal entries were doubled
    n = np.diag([0., 1.])
    h += -mu * (np.kron(n, I2) + np.kron(I2, n)) / 2
    return h


def to_pauli(h):
    """J_ab and B_a with h = sum J_ab s^a s^b + sum B_a (s^a x 1 + 1 x s^a)/2 + const + telescoping."""
    J = np.array([[np.trace(h @ np.kron(PAULI[a], PAULI[b])).real / 4 for b in range(3)] for a in range(3)])
    bl = np.array([np.trace(h @ np.kron(PAULI[a], I2)).real / 4 for a in range(3)])
    br = np.array([np.trace(h @ np.kron(I2, PAULI[a])).real / 4 for a in range(3)])
    return J, bl + br


def from_pauli(J, B):
    J = np.asarray(J, float); B = np.asarray(B, float)
    h = np.zeros((4, 4), complex)
    for a in range(3):
        h += B[a] * 0.5 * (np.kron(PAULI[a], I2) + np.kron(I2, PAULI[a]))
        for b in range(3):
            h += J[a, b] * np.kron(PAULI[a], PAULI[b])
    return h


def tda_levels(K, t, phi, V, mu, gL, gR, R=600, nlev=2):
    """Lowest nlev eigenvalues of the N<=2 tridiagonal problem at total momentum K and the
    lower edge of the two-particle continuum at that K."""
    A = abs(gL * np.exp(1j * K / 2) + gR * np.exp(-1j * K / 2))
    c = abs(2 * t * np.cos(K / 2 + phi))
    d = np.full(R + 1, -2.0 * mu); d[0] = -2 * t * np.cos(K + phi) - mu; d[1] += V
    e = np.full(R, c); e[0] = A
    w = eigvalsh_tridiagonal(d, e, select='i', select_range=(0, nlev - 1))
    return w, -2 * mu - 2 * c


def band_scan(t, phi, V, mu, gL, gR, nK=721, R=600, nlev=2):
    Ks = np.linspace(-np.pi, np.pi, nK)
    W = np.zeros((nK, nlev)); edge = np.zeros(nK)
    for i, K in enumerate(Ks):
        W[i], edge[i] = tda_levels(K, t, phi, V, mu, gL, gR, R, nlev)
    return Ks, W, edge


def _pivots(x, K, t, phi, V, mu, gL, gR):
    """Exact (semi-infinite, no cutoff) LDL pivots of T - x for x below the continuum edge.
    The uniform tail r >= 2 is positive definite there and contributes no negative pivot;
    its Schur complement onto r = 1 is (x' - sqrt(x'^2 - 4c^2))/2 with x' = x + 2mu."""
    A2 = np.abs(gL * np.exp(1j * K / 2) + gR * np.exp(-1j * K / 2)) ** 2
    c = 2 * t * np.cos(K / 2 + phi)
    xp = x + 2 * mu
    p1 = V - (xp - np.sqrt(np.maximum(xp ** 2 - 4 * c ** 2, 0))) / 2
    p1 = np.where(np.abs(p1) < 1e-300, 1e-300, p1)
    p0 = (-2 * t * np.cos(K + phi) - mu - x) - A2 / p1
    return p0, p1


def isolated_levels(K, t, phi, V, mu, gL, gR, iters=60):
    """Isolated eigenvalues below the two-particle continuum of the N <= 2 problem, exactly
    (infinite chain), vectorised over the array K.  Returns (E1, E2, edge); E1, E2 are nan
    where fewer isolated levels exist.  Uses Sturm counting with bisection."""
    K = np.asarray(K, float)
    c = 2 * t * np.cos(K / 2 + phi)
    edge = -2 * mu - 2 * np.abs(c)
    lo = edge - 4 * (abs(t) + abs(V) + abs(gL) + abs(gR) + abs(mu)) - 10

    def count(x):
        p0, p1 = _pivots(x, K, t, phi, V, mu, gL, gR)
        return (p0 < 0).astype(int) + (p1 < 0).astype(int)

    ntot = count(edge - 1e-12)
    out = []
    for m in (1, 2):
        a = lo.copy(); b = edge - 1e-12
        for _ in range(iters):
            mid = (a + b) / 2
            ok = count(mid) >= m
            b = np.where(ok, mid, b); a = np.where(ok, a, mid)
        out.append(np.where(ntot >= m, (a + b) / 2, np.nan))
    return out[0], out[1], edge


def pair_band_u1(K, t, phi, V, mu):
    """Exact U(1) adjacent-pair bound state E(K) = V + c^2/V - 2mu, c = 2t cos(K/2+phi),
    normalizable for V < 0 and |V| > |c|; nan otherwise."""
    c = 2 * t * np.cos(K / 2 + phi)
    E = V + c ** 2 / V - 2 * mu
    return np.where((V < 0) & (abs(V) > abs(c)), E, np.nan)
