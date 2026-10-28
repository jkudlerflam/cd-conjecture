"""Closed-book finite-dimensional and numerical checks, numpy/scipy only."""
import numpy as np
from scipy.integrate import quad
from scipy.special import digamma


def comm(a, b):
    return a @ b - b @ a


def su_generators(n):
    out = []
    for i in range(n):
        for j in range(i + 1, n):
            a = np.zeros((n, n), dtype=complex)
            a[i, j] = a[j, i] = 0.5
            out.append(a)
            a = np.zeros((n, n), dtype=complex)
            a[i, j], a[j, i] = -0.5j, 0.5j
            out.append(a)
    for m in range(1, n):
        a = np.zeros(n)
        a[:m], a[m] = 1, -m
        out.append(np.diag(a) / np.sqrt(2 * m * (m + 1)))
    return out


for n in (2, 3, 4):
    ts = su_generators(n)
    gram = np.array([[np.trace(a @ b) for b in ts] for a in ts])
    assert np.max(np.abs(gram - np.eye(n*n-1)/2)) < 1e-13
    p = np.zeros((n*n, n*n))
    for i in range(n):
        for j in range(n):
            p[j*n+i, i*n+j] = 1
    d_op = sum(-comm(np.kron(a, np.eye(n)),
                     comm(np.kron(a, np.eye(n)), p)) for a in ts)
    error = np.max(np.abs(d_op - (np.eye(n*n) - n*p)))
    assert error < 1e-13
    e = 1 + 2/n * (digamma(1/n) + np.euler_gamma)
    d_expect = 1 - n*e
    c2 = (n*n-1)/(2*n)
    print(f"SU({n}): Fierz/commutator error {error:.3g}; "
          f"e0={e:.12f}, D={d_expect:.12f}, "
          f"f-sum ratio={(n*n-1)/(n*d_expect):.12f}, "
          f"operator ratio={d_expect/(4*c2):.12f}, final ratio=0.5")
    if n == 2:
        singlet = np.array([0, 1, -1, 0]) / np.sqrt(2)
        saturation = np.vdot(singlet, d_op @ singlet).real
        assert abs(saturation - 4*c2) < 1e-13
        print(f"SU(2) isolated singlet: D={saturation:.12f}=4 C2 ||P||")

# Exact spectral f-sum on a six-site SU(2) permutation ring.
length = 6
dim = 2**length
h = np.zeros((dim, dim))
for state in range(dim):
    for x in range(length):
        y = (x + 1) % length
        flipped = state
        if ((state >> x) & 1) != ((state >> y) & 1):
            flipped ^= (1 << x) | (1 << y)
        h[flipped, state] += 1
energies, vectors = np.linalg.eigh(h)
ground = vectors[:, 0]
e_bond = energies[0]/length
for m in (1, 2):
    q = 2*np.pi*m/length
    rho = np.array([sum(np.exp(-1j*q*x) * (((state >> x) & 1)-0.5)
                       for x in range(length)) for state in range(dim)])
    weights = np.abs(vectors.conj().T @ (rho * ground))**2 / length
    moment = np.dot(energies - energies[0], weights)
    predicted = (1-np.cos(q)) * (1-2*e_bond)/3
    error = abs(moment-predicted)
    assert error < 1e-12
    print(f"L=6, SU(2), q={q:.9f}: Lehmann={moment:.12f}, "
          f"bond={predicted:.12f}, error={error:.3g}")

# Exact thermodynamic particle-hole integrals, t=1.
for kf in (0.3, np.pi/2, 2.4):
    for q in (0.1, 0.001):
        moment = quad(lambda p: (-2*np.cos(p+q)+2*np.cos(p))/(2*np.pi),
                      kf-q, kf, epsabs=1e-14)[0]
        exact = 4*np.sin(kf)*np.sin(q/2)**2/np.pi
        assert abs(moment-exact) < 1e-13
        current = (2*np.sin(kf))*q*q/(2*np.pi)
        print(f"U(1), kF={kf:.9f}, q={q:g}: exact error={abs(moment-exact):.3g}, "
              f"moment/current={moment/current:.12f}")
