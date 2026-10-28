#!/usr/bin/env python3
"""The route-A demonstrator bond on C^10, and its cell.

Local space C^10 = P (+) B.
  P = C^2, letters 0,1: one port qubit p.
  B = C^8, letters 2..9: three qubits (sL, c, sR), letter = 2 + 4 sL + 2 c + sR.
      sL, sR are the partners of the left and right ports, c is the logical qubit.
Qubit basis |0> (Z = +1), |1> (Z = -1); s+ = |0><1|.

U(1) charge  Q = sum_ports Z/2 - sum_partners Z/2 + sum_logical Z/2.
Pair states of (p, s): Phi+ = (|00>+|11>)/sqrt2 (charge 0), Phi- = (|00>-|11>)/sqrt2 (0),
|01> (charge +1), |10> (charge -1).

Port check on a (port, partner, logical) triple:
  H_pair = D0 |Phi-><Phi-| + D1 (|01><01| + |10><10|)
           + kap ( |Phi-><01| (x) s+_c  -  |Phi-><10| (x) s-_c  + h.c. )
It annihilates Phi+ (x) anything and is >= 0 iff D0 D1 >= kap^2.

Bond h on C^10 (x) C^10 (left site first):
  K n_B n_B                                   forbids BB
  -(mu/2)(n_B (x) 1 + 1 (x) n_B)              chemical potential of B
  [P,B]: H_pair(p = left site, s = sL, c) (x) 1_sR     left-port check
  [B,P]: H_pair(p = right site, s = sR, c) (x) 1_sL    right-port check
  [P,P]: lam [ Jp (s+ s- + s- s+) + Jz Z Z ]          port-port coupling

Cell = (P B P) = (pL, sL, c, sR, pR), dimension 32, ordering pL sL c sR pR (pL most significant).
Claude (Opus 5.5), 2026-10-27, for contrib/routeA_demonstrator.
"""
import numpy as np

I2 = np.eye(2)
SP = np.array([[0., 1.], [0., 0.]])   # |0><1|, raises Z
SM = SP.T.copy()
Z = np.diag([1., -1.])
X = SP + SM


def kron(*a):
    out = np.array([[1.]])
    for x in a:
        out = np.kron(out, x)
    return out


def pair_states():
    e = np.eye(4)                      # (p, s) -> index 2p + s
    phip = (e[0] + e[3]) / np.sqrt(2)
    phim = (e[0] - e[3]) / np.sqrt(2)
    return phip, phim, e[1], e[2]      # Phi+, Phi-, |01>, |10>


def h_pair(D0=1.0, D1=1.0, kap=0.5):
    """Operator on (p, s, c), ordering p s c, dimension 8."""
    phip, phim, v01, v10 = pair_states()
    P = lambda a, b: np.outer(a, b)
    H = np.kron(D0 * P(phim, phim) + D1 * (P(v01, v01) + P(v10, v10)), I2)
    T = kap * (np.kron(P(phim, v01), SP) - np.kron(P(phim, v10), SM))
    return H + T + T.T


def perm_op(op, n, perm):
    """op acts on qubits in order (0..n-1); return it acting with qubit k placed at position perm[k]."""
    T = op.reshape([2] * (2 * n))
    inv = np.argsort(perm)
    T = T.transpose(list(inv) + [n + i for i in inv])
    return T.reshape(2 ** n, 2 ** n)


def bond(D0=1.0, D1=1.0, kap=0.5, lam=0.2, Jp=1.0, Jz=1.0, mu=1.0, K=4.0):
    d = 10
    h = np.zeros((d * d, d * d))
    nP = np.diag([1.] * 2 + [0.] * 8)
    nB = np.diag([0.] * 2 + [1.] * 8)
    h += K * np.kron(nB, nB) - 0.5 * mu * (np.kron(nB, np.eye(d)) + np.kron(np.eye(d), nB))
    hp = h_pair(D0, D1, kap)                                    # (p, s, c)
    # left check: physical order (p | sL c sR); hp acts on (p, sL, c), identity on sR
    left = np.kron(hp, I2)                                      # order p sL c sR already
    # right check: physical order (sL c sR | p); hp on (p, sR, c)
    right = perm_op(np.kron(hp, I2), 4, [3, 2, 1, 0])           # qubits (p,s,c,x) -> positions (3,2,1,0)
    # careful: we want (p, s=sR, c) with sL the spectator; physical order sL c sR p.
    # perm maps operator qubit k to physical position perm[k]: p->3, s->2, c->1, spectator->0.
    V = lam * (Jp * (np.kron(SP, SM) + np.kron(SM, SP)) + Jz * np.kron(Z, Z))
    # embed: P letters are 0,1 ; B letters 2..9
    iP = np.arange(2); iB = np.arange(2, 10)

    def put(blockrows, blockcols, op):
        idx = np.array([a * d + b for a in blockrows[0] for b in blockrows[1]])
        h[np.ix_(idx, idx)] += op

    put((iP, iB), None, left)
    put((iB, iP), None, right)
    put((iP, iP), None, V)
    return h


def charge_site():
    """Diagonal charge of each of the 10 letters."""
    q = np.zeros(10)
    q[0], q[1] = 0.5, -0.5
    for sL in range(2):
        for c in range(2):
            for sR in range(2):
                zs = lambda b: 1 - 2 * b
                q[2 + 4 * sL + 2 * c + sR] = -zs(sL) / 2 + zs(c) / 2 - zs(sR) / 2
    return q


# ----------------------------------------------------------------------------- the cell (P B P)
def cell_ops(D0=1.0, D1=1.0, kap=0.5):
    """Cell Hamiltonian H0 on (pL, sL, c, sR, pR) and the port operators."""
    hp = h_pair(D0, D1, kap)                                      # (p, s, c)
    HL = perm_op(np.kron(hp, kron(I2, I2)), 5, [0, 1, 2, 3, 4])   # (pL, sL, c) then sR, pR
    HR = perm_op(np.kron(hp, kron(I2, I2)), 5, [4, 3, 2, 0, 1])   # p->pR(4), s->sR(3), c->2, spectators 0,1
    H0 = HL + HR
    ports = dict(
        L=dict(sp=kron(SP, I2, I2, I2, I2), sm=kron(SM, I2, I2, I2, I2), z=kron(Z, I2, I2, I2, I2)),
        R=dict(sp=kron(I2, I2, I2, I2, SP), sm=kron(I2, I2, I2, I2, SM), z=kron(I2, I2, I2, I2, Z)))
    phip = pair_states()[0]
    code = np.zeros((32, 2))
    for c in range(2):
        ec = I2[c]
        # order pL sL c sR pR ; Phi+ on (pL,sL) and on (pR,sR): (sR,pR) Bell is symmetric
        code[:, c] = np.kron(np.kron(phip, ec), phip)
    logical = dict(z=code @ np.diag([1., -1.]) @ code.T)
    qc = np.zeros(32)
    for i in range(32):
        b = [(i >> (4 - k)) & 1 for k in range(5)]
        zz = [1 - 2 * x for x in b]
        qc[i] = (zz[0] - zz[1] + zz[2] - zz[3] + zz[4]) / 2
    return H0, ports, code, qc


def inter_cell(lam=0.2, Jp=1.0, Jz=1.0, D0=1.0, D1=1.0, kap=0.5):
    """Two cells: H = H0 (x) 1 + 1 (x) H0 + lam V(pR of cell 1, pL of cell 2)."""
    H0, ports, code, qc = cell_ops(D0, D1, kap)
    I = np.eye(32)
    V = (Jp * (np.kron(ports['R']['sp'], ports['L']['sm']) + np.kron(ports['R']['sm'], ports['L']['sp']))
         + Jz * np.kron(ports['R']['z'], ports['L']['z']))
    return np.kron(H0, I) + np.kron(I, H0), V, code, qc


if __name__ == '__main__':
    # sanity: bond is symmetric, charge conserving
    h = bond()
    assert np.allclose(h, h.T)
    q = charge_site()
    Q2 = np.kron(np.diag(q), np.eye(10)) + np.kron(np.eye(10), np.diag(q))
    print('bond charge commutator', np.abs(h @ Q2 - Q2 @ h).max())
    H0, ports, code, qc = cell_ops()
    print('cell hermitian', np.allclose(H0, H0.T), 'charge comm', np.abs(H0 @ np.diag(qc) - np.diag(qc) @ H0).max())
    w = np.linalg.eigvalsh(H0)
    print('cell spectrum', np.round(w[:8], 6))
    print('code energy', np.abs(H0 @ code).max())
