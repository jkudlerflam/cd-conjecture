#!/usr/bin/env python3
"""C172 (Claude Opus 5.5 agent, 2026-11-08): ring ED for the two-kink test of PROOF section 4.

H = sum_j [ sigma_j^T J sigma_{j+1} + B . sigma_j ] on a ring of L qubits (Pauli matrices).
Lowest eigenvalues by Lanczos with a matrix-free matvec; lattice momentum of each eigenvector from <psi|T|psi>.

Kink gap on odd rings:   D_K(L) = E0(L) - [E0(L-1) + E0(L+1)]/2          (L odd)
Kink-even gap on even rings: D_2(L) = E_2(L) - [E_0(L) + E_1(L)]/2      (L even; E_0, E_1 the Neel doublet)
Binding test: b(L) = D_2(L) - 2 D_K(L +- 1). Unbound kinks: b > 0, falling like 1/L^2. Bound pair: b -> negative const.

usage: python3 ed_kink.py 'a=0.9,Dz=2.0,g=0.3,hx=0.2,hz=0.05,Dx=0.03' Lmin Lmax [nev]
"""
import sys, json
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh


def fam(a=0.8, Dz=1.0, g=0.0, gxy=0.0, G=0.0, Gp=None, Dx=0.0, Dy=0.0, hx=0.0, hy=0.0, hz=0.0):
    """Identical to codex/work_C165/classical.fam (copied, not imported)."""
    if Gp is None:
        Gp = G
    J = np.array([[-np.cos(a) + g, np.sin(a) + gxy, G],
                  [-np.sin(a) + gxy, -np.cos(a) - g, Dx],
                  [Gp, -Dx, Dz]], float)
    J[0, 2] += -Dy; J[2, 0] += Dy
    return J, np.array([hx, hy, hz], float)


PAULI = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]], complex),
         np.array([[1, 0], [0, -1]], complex)]


def bond_matrix(J, B):
    """4x4 bond matrix in basis |b_j b_{j+1}> with b = 0 for up (sigma^z = +1). Field put on site j."""
    h = np.zeros((4, 4), complex)
    for a in range(3):
        for b in range(3):
            h += J[a, b] * np.kron(PAULI[a], PAULI[b])
        h += B[a] * np.kron(PAULI[a], np.eye(2))
    return h


def make_matvec(L, h):
    N = 1 << L
    s = np.arange(N, dtype=np.int64)
    terms = []
    for j in range(L):
        k = (j + 1) % L
        bj = (s >> j) & 1; bk = (s >> k) & 1
        idx_in = 2 * bj + bk
        for i in range(4):
            mask = np.nonzero(idx_in == i)[0]
            for o in range(4):
                amp = h[o, i]
                if abs(amp) < 1e-15:
                    continue
                flip = ((((o >> 1) ^ (i >> 1)) & 1) << j) | ((((o & 1) ^ (i & 1)) & 1) << k)
                terms.append((mask, mask ^ flip, amp))
    def mv(v):
        v = v.ravel()
        out = np.zeros(N, complex)
        for src, dst, amp in terms:
            out[dst] += amp * v[src]   # dst is a permutation of src within each term, no repeats
        return out
    return LinearOperator((N, N), matvec=mv, dtype=complex)


def translate(v, L):
    """(T v)(s) = v(rot^{-1} s): T moves site j to j+1."""
    N = 1 << L
    s = np.arange(N, dtype=np.int64)
    rot = ((s << 1) | (s >> (L - 1))) & (N - 1)  # site j -> j+1
    out = np.empty_like(v)
    out[rot] = v
    return out


def spectrum(J, B, L, nev=8):
    h = bond_matrix(J, B)
    A = make_matvec(L, h)
    rng = np.random.default_rng(1)
    v0 = rng.standard_normal(1 << L) + 1j * rng.standard_normal(1 << L)
    w, V = eigsh(A, k=nev, which='SA', v0=v0, tol=1e-10)
    o = np.argsort(w); w = w[o]; V = V[:, o]
    ks = []
    for n in range(nev):
        ov = np.vdot(V[:, n], translate(V[:, n], L))
        ks.append((float(np.angle(ov)) % (2 * np.pi), float(abs(ov))))
    return w, ks


if __name__ == '__main__':
    kw = eval('dict(' + sys.argv[1] + ')')
    Lmin, Lmax = int(sys.argv[2]), int(sys.argv[3])
    nev = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    J, B = fam(**kw)
    res = {}
    for L in range(Lmin, Lmax + 1):
        w, ks = spectrum(J, B, L, nev)
        res[L] = dict(E=[float(x) for x in w], k=ks)
        print(L, ' '.join(f'{x:.6f}' for x in w[:6]), ' k/pi:', ' '.join(f'{k[0]/np.pi:.3f}({k[1]:.2f})' for k in ks[:6]), flush=True)
    out = {'kw': sys.argv[1], 'res': res}
    # gaps
    rows = []
    for L in range(Lmin + 1, Lmax):
        if L % 2 == 1 and (L - 1) in res and (L + 1) in res:
            DK = res[L]['E'][0] - 0.5 * (res[L - 1]['E'][0] + res[L + 1]['E'][0])
            rows.append(('K', L, DK))
    for L in range(Lmin, Lmax + 1):
        if L % 2 == 0:
            E = res[L]['E']
            rows.append(('2', L, E[2] - 0.5 * (E[0] + E[1]), E[1] - E[0]))
    out['gaps'] = rows
    for r in rows:
        print(*r)
    tag = sys.argv[5] if len(sys.argv) > 5 else 'run'
    json.dump(out, open(f'ed_{tag}.json', 'w'), indent=1)
