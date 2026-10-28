#!/usr/bin/env python3
"""Two-flavour route-A demonstrator: a 13-letter bond with exact U(1) x U(1).

Local space C^13 = P (+) B.
  P = C^3 (x) C^3, letters 0..8 = 3x + y.  x is the port qutrit, y its partner.  Charges: x = 0, 1, 2
      carries 0, e1, e2; y carries minus that.  The port and partner sit in the same letter.
  B = C^2 (x) C^2, letters 9..12 = 9 + 2 n1 + n2: two logical modes (occupations n1, n2), charge n1 e1 + n2 e2.
Bell states on (x, y):  Phi_L = sum_x |x x>/sqrt3  (left ports),  Phi_R = sum_x w^x |x x>/sqrt3, w = e^{2 pi i/3}
(right ports).  They are orthogonal, so a P shared by two B's is frustrated.
Neutral syndromes N^s_a = (D_a (x) 1) Phi_s with D_1 = sqrt(3/2) diag(1,-1,0), D_2 = diag(1,1,-2)/sqrt2
(tr D_a D_b / 3 = delta_ab).  Charged syndromes |a b>, a != b.

Port check on (x, y, n1, n2), s = L or R (dimension 36):
  H_s = DN Pi_N^s + sum_{a != b} DC_{ab} |ab><ab|
        + sum_a kap_a ( |N^s_a><a0| (x) c_a^dag  +  sg |N^s_a><0a| (x) c_a  + h.c. )
The charge +e_a syndrome |a0> converts to N_a while flavour a is created, the charge -e_a syndrome |0a>
while it is annihilated.  H_s annihilates Phi_s (x) anything and is >= 0 if DC DN >= kap^2.
Bond h on C^13 (x) C^13:
  K n_B n_B  -  (mu/2)(n_B 1 + 1 n_B)
  (P,B): H_L(x, y of the P ; n1 n2 of the B)      (B,P): H_R(x, y of the P ; n1 n2 of the B)
  (P,P): lam [ Jp sum_a (|a><0| (x) |0><a| + h.c.)  +  Jz sum_a D_a (x) D_a ]  on (x, x')
  B on-site: -(h2/2)(...) logical field on flavour 2 (optional, a logical-scale chemical potential).
Claude (Opus 5.5), 2026-10-28, for contrib/routeA_demonstrator_n2.
"""
import numpy as np

W3 = np.exp(2j * np.pi / 3)
D1 = np.sqrt(1.5) * np.diag([1., -1., 0.])
D2 = np.diag([1., 1., -2.]) / np.sqrt(2)
DS = [D1, D2]
I2, I3 = np.eye(2), np.eye(3)
SP = np.array([[0., 0.], [1., 0.]])     # c^dag on an occupation qubit: |1><0|  (basis n = 0, 1)
SM = SP.T.copy()
NOP = np.diag([0., 1.])


def ket(*idx, dims):
    v = np.zeros(int(np.prod(dims)), complex)
    i = 0
    for k, dd in zip(idx, dims):
        i = i * dd + k
    v[i] = 1
    return v


def bell(side):
    ph = [1, 1, 1] if side == 'L' else [1, W3, W3 ** 2]
    return sum(ph[x] * ket(x, x, dims=(3, 3)) for x in range(3)) / np.sqrt(3)


def neutral(side):
    Ph = bell(side)
    return [np.kron(D, I3) @ Ph for D in DS]


def logical_ops():
    """c_a, c_a^dag on (n1, n2) (spin language; flavours are separately conserved, no JW string needed)."""
    c = [np.kron(SM, I2), np.kron(I2, SM)]
    return c, [x.T for x in c]


def h_pair(side, DN=1.0, DC=1.0, kap=(0.6, 0.6), sg=-1.0, DCd=(0.0, 0.0)):
    """Check on (x, y, n1, n2); dimension 36.  DCd[a] shifts the energy of |0a> (charge -e_a) relative
    to |a0>, which breaks particle-hole symmetry of flavour a and gives it a logical field at second order."""
    N = neutral(side)
    PN = sum(np.outer(v, v.conj()) for v in N)
    H = DN * np.kron(PN, np.eye(4))
    for a in range(3):
        for b in range(3):
            if a != b:
                e = DC
                if b != 0 and a == 0:
                    e = DC + DCd[b - 1]
                v = ket(a, b, dims=(3, 3))
                H += e * np.kron(np.outer(v, v), np.eye(4))
    c, cd = logical_ops()
    for a in (0, 1):
        x = a + 1
        up = np.kron(np.outer(N[a], ket(x, 0, dims=(3, 3)).conj()), cd[a])
        dn = np.kron(np.outer(N[a], ket(0, x, dims=(3, 3)).conj()), c[a])
        T = kap[a] * (up + sg * dn)
        H += T + T.conj().T
    return H


def port_V(Jp=1.0, Jz=1.0, gauge=(0.0, 0.0)):
    """gauge = (phi1, phi2): the exact unitary U = prod_j exp(i sum_a phi_a j N_a,j) maps H to the same bond with
    the flavour-a port hopping multiplied by exp(-i phi_a); used only to choose a convenient frame for VUMPS."""
    V = np.zeros((9, 9), complex)
    for a in (1, 2):
        up = np.zeros((3, 3)); up[a, 0] = 1
        ph = np.exp(-1j * gauge[a - 1])
        V += Jp * (ph * np.kron(up, up.T) + np.conj(ph) * np.kron(up.T, up))
    for D in DS:
        V += Jz * np.kron(D, D)
    return V


def charges():
    """(N1, N2) of each of the 13 letters."""
    q = np.zeros((13, 2))
    for x in range(3):
        for y in range(3):
            v = np.zeros(2)
            if x: v[x - 1] += 1
            if y: v[y - 1] -= 1
            q[3 * x + y] = v
    for n1 in range(2):
        for n2 in range(2):
            q[9 + 2 * n1 + n2] = (n1, n2)
    return q


def bond(lam=0.1, mu=0.5, K=2.0, DN=1.0, DC=1.0, kap=(0.6, 0.6), sg=-1.0, DCd=(0.0, 0.0), Jp=1.0, Jz=1.0, h2=0.0,
         gauge=(0.0, 0.0)):
    d = 13
    h = np.zeros((d * d, d * d), complex)
    nB = np.diag([0.] * 9 + [1.] * 4)
    Id = np.eye(d)
    h += K * np.kron(nB, nB) - 0.5 * mu * (np.kron(nB, Id) + np.kron(Id, nB))
    # logical field on flavour 2, split over the two bonds of each B
    fB = np.zeros((d, d)); fB[9:, 9:] = -h2 * np.kron(I2, NOP)
    h += 0.5 * (np.kron(fB, Id) + np.kron(Id, fB))
    iP, iB = np.arange(9), np.arange(9, 13)
    HL = h_pair('L', DN, DC, kap, sg, DCd)          # (x,y | n1,n2): P letter then B letter
    HR = h_pair('R', DN, DC, kap, sg, DCd)          # (x,y | n1,n2); physical order is B then P
    HRp = HR.reshape(9, 4, 9, 4).transpose(1, 0, 3, 2).reshape(36, 36)
    V = lam * port_V(Jp, Jz, gauge)

    def put(rows, cols, op):
        idx = np.array([a * d + b for a in rows for b in cols])
        h[np.ix_(idx, idx)] += op
    put(iP, iB, HL)
    put(iB, iP, HRp)
    # V acts on the x factors of two P letters, identity on their y factors
    VP = np.einsum('acbd,ef,gh->aecgbfdh', V.reshape(3, 3, 3, 3), I3, I3).reshape(81, 81)
    put(iP, iP, VP)
    return h


# ----------------------------------------------------------------------------- the cell (P B P), 324 states
def cell_ops(DN=1.0, DC=1.0, kap=(0.6, 0.6), sg=-1.0, DCd=(0.0, 0.0), h2=0.0):
    """Order (xL yL | n1 n2 | xR yR), dims (9, 4, 9)."""
    HL = h_pair('L', DN, DC, kap, sg, DCd)           # (xL yL n1 n2)
    HR = h_pair('R', DN, DC, kap, sg, DCd)           # (xR yR n1 n2)
    H0 = np.kron(HL, np.eye(9))
    # explicit: HR acts on (B, P_R) with B in the middle position
    T = np.zeros((9, 4, 9, 9, 4, 9), complex)
    HRb = HR.reshape(9, 4, 9, 4)                     # [pR, b, pR', b']
    for i in range(9):
        T[i, :, :, i, :, :] += HRb.transpose(1, 0, 3, 2)  # [b, pR, b', pR']
    H0 += T.reshape(324, 324)
    H0 -= h2 * np.kron(np.kron(np.eye(9), np.kron(I2, NOP)), np.eye(9))
    xop = lambda M: np.kron(M, I3)                   # operator on x of a P letter (9 x 9)
    ports = {}
    for side in ('L', 'R'):
        ops = {}
        for a in (1, 2):
            up = np.zeros((3, 3)); up[a, 0] = 1
            ops[f'up{a}'] = xop(up); ops[f'dn{a}'] = xop(up.T)
        for k, D in enumerate(DS):
            ops[f'D{k + 1}'] = xop(D)
        if side == 'L':
            ports['L'] = {k: np.kron(v, np.eye(36)) for k, v in ops.items()}
        else:
            ports['R'] = {k: np.kron(np.eye(36), v) for k, v in ops.items()}
    code = np.zeros((324, 4), complex)
    for l in range(4):
        code[:, l] = np.kron(np.kron(bell('L'), np.eye(4)[l]), bell('R'))
    q = charges()
    qc = np.zeros((324, 2))
    for i in range(324):
        a, rest = divmod(i, 36)
        b, c = divmod(rest, 9)
        qc[i] = q[a] + q[9 + b] + q[c]
    return H0, ports, code, qc


def bond_terms(Jp=1.0, Jz=1.0):
    """Inter-cell coupling as a list (coef, opR_name, opL_name)."""
    t = []
    for a in (1, 2):
        t += [(Jp, f'up{a}', f'dn{a}'), (Jp, f'dn{a}', f'up{a}')]
    t += [(Jz, 'D1', 'D1'), (Jz, 'D2', 'D2')]
    return t


if __name__ == '__main__':
    h = bond()
    print('bond hermitian', np.abs(h - h.conj().T).max())
    q = charges()
    for k in range(2):
        Q = np.kron(np.diag(q[:, k]), np.eye(13)) + np.kron(np.eye(13), np.diag(q[:, k]))
        print('charge', k + 1, 'commutator', np.abs(h @ Q - Q @ h).max())
    H0, ports, code, qc = cell_ops()
    w = np.linalg.eigvalsh(H0)
    print('cell spectrum', np.round(w[:10], 6), ' code residual', np.abs(H0 @ code).max())
    print('min eigenvalue of each check', np.linalg.eigvalsh(h_pair('L'))[:3], np.linalg.eigvalsh(h_pair('R'))[:3])
