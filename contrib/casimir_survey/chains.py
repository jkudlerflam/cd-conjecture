#!/usr/bin/env python3
"""Bonds h (d^2 x d^2, site 1 the most significant factor) for the survey."""
import numpy as np

X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1.0 + 0j, -1])
I2 = np.eye(2)


def spin(S):
    m = np.arange(S, -S - 1, -1)
    sp = np.zeros((len(m), len(m)), complex)
    for i in range(1, len(m)):
        sp[i - 1, i] = np.sqrt(S * (S + 1) - m[i] * (m[i] + 1))
    return (sp + sp.conj().T) / 2, (sp - sp.conj().T) / 2j, np.diag(m).astype(complex)


def swap(d):
    return np.eye(d * d).reshape(d, d, d, d).transpose(1, 0, 2, 3).reshape(d * d, d * d).astype(complex)


def graded_swap(par):
    """P_g |a b> = (-1)^{p_a p_b} |b a>, par = list of parities (0 boson, 1 fermion)."""
    d = len(par)
    P = np.zeros((d * d, d * d), complex)
    for a in range(d):
        for b in range(d):
            P[b * d + a, a * d + b] = (-1) ** (par[a] * par[b])
    return P


def number_op(d, which):
    n = np.zeros((d, d)); n[which, which] = 1
    return n


def onsite(o, d):
    """(o (x) 1 + 1 (x) o)/2: an on-site term distributed over the two bonds of each site."""
    return (np.kron(o, np.eye(d)) + np.kron(np.eye(d), o)) / 2


def ss(S):
    s = spin(S)
    return sum(np.kron(a, a) for a in s)


def spin_projectors(S):
    d = int(round(2 * S + 1)); x = ss(S).real
    vals = [(J * (J + 1) - 2 * S * (S + 1)) / 2 for J in range(d)]
    Ps = []
    for j, v in enumerate(vals):
        p = np.eye(d * d)
        for k, w in enumerate(vals):
            if k != j:
                p = p @ (x - w * np.eye(d * d)) / (v - w)
        Ps.append(p.astype(complex))
    return Ps


def takhtajan_babujian(S):
    """h = sum_J H_J P_J (H_J harmonic numbers), the integrable spin-S chain, v = pi/2."""
    Ps = spin_projectors(S)
    Hn = [sum(1.0 / k for k in range(1, J + 1)) for J in range(len(Ps))]
    return sum(c * p for c, p in zip(Hn, Ps))


def blbq(theta):
    x = ss(1.0)
    return np.cos(theta) * x + np.sin(theta) * (x @ x)


def models():
    """Calibration chains: name -> (h, d, c, A_exact or None)."""
    M = {}
    M['XX'] = (0.5 * (np.kron(X, X) + np.kron(Y, Y)), 2, 1.0, 2.0)
    M['Heisenberg'] = (swap(2), 2, 1.0, np.pi)
    M['Ising'] = (-np.kron(Z, Z) - 0.5 * (np.kron(X, I2) + np.kron(I2, X)), 2, 0.5, 1.0)
    M['SU3'] = (swap(3), 3, 2.0, 2 * 2 * np.pi / 3)
    M['SU4'] = (swap(4), 4, 3.0, 3 * 2 * np.pi / 4)
    M['TB1'] = (takhtajan_babujian(1.0), 3, 1.5, 1.5 * np.pi / 2)
    return M


def sym_rep(N, k):
    """SU(N) on Sym^k(C^N), bosonic occupation basis. Returns (basis, E) with E[a][b] = b_a^dag b_b."""
    import itertools
    basis = [c for c in itertools.product(range(k + 1), repeat=N) if sum(c) == k]
    basis.sort(reverse=True)
    idx = {c: i for i, c in enumerate(basis)}
    d = len(basis)
    E = [[np.zeros((d, d)) for _ in range(N)] for _ in range(N)]
    for i, c in enumerate(basis):
        for a in range(N):
            for b in range(N):
                if c[b] == 0:
                    continue
                n = list(c); amp = np.sqrt(n[b]); n[b] -= 1; amp *= np.sqrt(n[a] + 1); n[a] += 1
                E[a][b][idx[tuple(n)], i] += amp
    return basis, E


def sym_rep_chain(N, k, check=True):
    """Integrable SU(N) chain in Sym^k: h = sum_j H_{k-j} P_{(2k-j, j)}, ordered by the two-site
    Casimir X = sum_ab E_ab (x) E_ba (largest X on the fully symmetric component). For k = 1 this
    is (1 + P)/2. Returns (h, d, charges) with U(1)^(N-1) colour charges."""
    basis, E = sym_rep(N, k)
    d = len(basis)
    Xop = sum(np.kron(E[a][b], E[b][a]) for a in range(N) for b in range(N))
    ev = np.linalg.eigvalsh(Xop)
    vals = []
    for x in sorted(ev, reverse=True):
        if not vals or abs(x - vals[-1]) > 1e-8:
            vals.append(x)
    assert len(vals) == k + 1, vals
    h = np.zeros((d * d, d * d))
    for j, v in enumerate(vals):
        p = np.eye(d * d)
        for w in vals:
            if w != v:
                p = p @ (Xop - w * np.eye(d * d)) / (v - w)
        Hk = sum(1.0 / l for l in range(1, k - j + 1))
        h += Hk * p
    charges = np.array([c[:N - 1] for c in basis], dtype=int)
    return h.astype(complex), d, charges
