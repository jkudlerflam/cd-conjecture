"""Leading (order g^2) quadratic symbol of C32's mobile 28-letter core, from
one-carrier graphs, without a 2^n logical space or a scalar subtraction.

Geometry (C32 PROOFS section 1). One excited carrier per cell lives on the
2n+2 vertices
    0 = A,  1..n = d_0..d_{n-1},  n+1..2n = x_0..x_{n-1},  2n+1 = B.
Every hop enters as its positive edge matrix, so each edge of amplitude c adds
|c| to both endpoint diagonals. On-site m everywhere, plus D on A, x, B.
After the Jordan-Wigner gauge (neutral token at data site j carries
Z_0...Z_{j-1}) and after moving the n neutral express phases e^{i phi} into the
B source, the neutral graph is real except for the loop edge d_{n-1}-B, whose
phase becomes psi_0 - n phi. Charged carriers have real express edges and loop
phase psi_c.

The pair formula is C33's equation (17) (work_C33/PROOF.md section 8 and
fast_loop.py), which holds for any carrier graph with this species/phase
structure at n phi = pi/2 and port phases (1, i, -i): the left cell's source
sits at B, the right cell's at A. validate() checks it against a literal
order-g^2 extraction in C32's complete excited pair space.

numpy and scipy only.
"""
import os
import sys

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu


def carrier(n, m=1.0, D=100.0, td=0.0, loopphase=0.0):
    S = 2 * n + 2
    H = np.zeros((S, S), complex)

    def edge(i, j, c):
        H[i, i] += abs(c)
        H[j, j] += abs(c)
        H[j, i] += c
        H[i, j] += np.conj(c)

    edge(0, 1, -1.0)                       # A_s R <-> A_0 d_s (left link)
    for j in range(n - 1):
        edge(1 + j, 2 + j, -1.0)           # data rail
    edge(1, n + 1, -1.0)                   # switch d_0 <-> x_0
    for j in range(n - 1):
        edge(n + 1 + j, n + 2 + j, -1.0)   # express rail (phases gauged)
    edge(2 * n, 2 * n + 1, -1.0)           # x_{n-1} -> B
    if td:
        edge(n, 2 * n + 1, -td * np.exp(1j * loopphase))  # loop d_{n-1} -> B
    H += m * np.eye(S)
    for v in [0] + list(range(n + 1, 2 * n + 2)):
        H[v, v] += D
    return sp.csc_matrix(H)


def leading(n, m=1.0, D=100.0, td=0.0, psi0=0.9, psic=0.0):
    """Return (T, M) of the order-g^2 symbol, visibility 1/(n+2)^2 included.
    T[i, j] moves a fermion from data site j of the right cell to data site i
    of the left cell; M is the sum of the two one-cell matrices."""
    phi = np.pi / (2 * n)
    S = 2 * n + 2
    H0 = carrier(n, m, D, td, psi0 - n * phi)
    Hc = carrier(n, m, D, td, psic)
    I = sp.identity(S, format='csc')

    def factor(A, B):
        return splu((sp.kron(A, I) + sp.kron(I, B)).tocsc())

    lu00, lucc, lu0c, luc0 = factor(H0, H0), factor(Hc, Hc), factor(H0, Hc), factor(Hc, H0)
    b = np.zeros(S * S, complex)
    b[(S - 1) * S + 0] = 1.0                # left carrier at B, right carrier at A
    x = lu00.solve(b).reshape(S, S)
    y = lucc.solve(b).reshape(S, S)

    def restrict(v):
        U = np.zeros((S * S, n), complex)
        V = U.copy()
        for j in range(n):
            U[(j + 1) * S:(j + 2) * S, j] = v[j + 1, :]
            V[(j + 1)::S, j] = v[:, j + 1]
        return U, V

    U, V = restrict(x)
    W, Z = restrict(y)
    RU, RV, RW, RZ = luc0.solve(U), lu0c.solve(V), lu0c.solve(W), luc0.solve(Z)
    vis = 1.0 / (n + 2) ** 2
    T = -2 * vis * np.real((Z.conj().T @ RU).T + W.conj().T @ RV)
    M = -2j * vis * np.imag(U.conj().T @ RU + W.conj().T @ RW + V.conj().T @ RV + Z.conj().T @ RZ)
    return T, M


def literal(n, m=1.0, D=2.0, td=0.0):
    """Order-g^2 coefficient K = -Z^dag R0 C R0 C R0 Z in C32's complete
    excited pair space (C32 check_asymptotic.py, extended to the onsite block
    and to a nonzero loop). C32's mobile.py fixes psi_0 = 0.9, psi_c = 0."""
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(here, '..', '..', 'codex', 'work_C32'))
    sys.dont_write_bytecode = True
    import mobile as mob
    c0 = mob.core(n, g=0, mass=m, express=D, loop=td, phase=np.pi / (2 * n))
    c1 = mob.core(n, g=1, mass=m, express=D, loop=td, phase=np.pi / (2 * n))
    ex = c0['ex']
    H = sp.csc_matrix(c0['H'][np.ix_(ex, ex)])
    C = sp.csc_matrix((c1['H'] - c0['H'])[np.ix_(ex, ex)])
    I = sp.identity(H.shape[0], format='csc')
    HH = (sp.kron(H, I) + sp.kron(I, H)).tocsc()
    CC = (sp.kron(C, I) + sp.kron(I, C)).tocsc()
    Zs = sum(v * np.kron((c0['F']['R'][fr] @ c0['G'])[ex], (c0['F']['L'][fl] @ c0['G'])[ex])
             for fr, fl, v in mob.TERMS)
    lu = splu(HH)
    X = lu.solve(Zs)
    X = lu.solve(CC @ X)
    X = lu.solve(CC @ X)
    K = -Zs.conj().T @ X
    ids = [2 ** (2 * n - 1 - j) for j in range(2 * n)]
    T = K[np.ix_(ids[:n], ids[n:])]
    Mfull = K[np.ix_(ids, ids)] - K[0, 0] * np.eye(2 * n)
    M = Mfull[:n, :n] + Mfull[n:, n:]
    return T, M, float(np.max(abs(K - K.conj().T)))


def validate():
    out = []
    for n, D, td in [(1, 2.0, 0.0), (1, 2.0, 0.05), (2, 2.0, 0.0), (2, 2.0, 0.05), (3, 2.0, 0.02)]:
        Tf, Mf = leading(n, m=1.0, D=D, td=td, psi0=0.9, psic=0.0)
        Tl, Ml, herm = literal(n, m=1.0, D=D, td=td)
        sc = np.max(abs(Tl))
        out.append(dict(n=n, D=D, td=td, T_scale=float(sc),
                        T_rel_err=float(np.max(abs(Tf - Tl)) / sc),
                        M_over_T=float(np.max(abs(Ml)) / sc),
                        M_err_over_T=float(np.max(abs(Mf - Ml)) / sc),
                        literal_hermiticity=herm))
    return out


if __name__ == '__main__':
    import json
    print(json.dumps(validate(), indent=1))
