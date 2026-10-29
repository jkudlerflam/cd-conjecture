#!/usr/bin/env python3
"""Bond families for the attack.  Every function returns (h, d, charges) with h the
d^2 x d^2 Hermitian bond (site i the most significant factor) and charges a d x nq
integer array of conserved additive site charges (nq may be 0).  Fermionic models are
written in the graded-permutation (Jordan-Wigner-free) form, which on a ring equals the
fermionic model up to the boundary sign, irrelevant for the Casimir coefficient."""
import numpy as np

X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1.0 + 0j, -1])


def spin(S):
    m = np.arange(S, -S - 1, -1)
    sp = np.zeros((len(m), len(m)), complex)
    for i in range(1, len(m)):
        sp[i - 1, i] = np.sqrt(S * (S + 1) - m[i] * (m[i] + 1))
    return (sp + sp.conj().T) / 2, (sp - sp.conj().T) / 2j, np.diag(m).astype(complex)


def swap(d):
    return np.eye(d * d).reshape(d, d, d, d).transpose(1, 0, 2, 3).reshape(d * d, d * d).astype(complex)


def graded_swap(par):
    d = len(par)
    P = np.zeros((d * d, d * d), complex)
    for a in range(d):
        for b in range(d):
            P[b * d + a, a * d + b] = (-1) ** (par[a] * par[b])
    return P


def onsite(o, d):
    return (np.kron(o, np.eye(d)) + np.kron(np.eye(d), o)) / 2


def ss(S):
    return sum(np.kron(a, a) for a in spin(S))


def proj(d, vec):
    v = np.asarray(vec, complex).reshape(-1); v = v / np.linalg.norm(v)
    return np.outer(v, v.conj())


# ------------------------------------------------------------ d = 2
def heisenberg():
    return ss(0.5), 2, np.array([[1], [0]])


def xxz(delta):
    sx, sy, sz = spin(0.5)
    return np.kron(sx, sx) + np.kron(sy, sy) + delta * np.kron(sz, sz), 2, np.array([[1], [0]])


# ------------------------------------------------------------ d = 3
def tj(mu, sign=-1):
    """Supersymmetric t-J: states (hole, up, down), hole bosonic.  h = sign*P_g - mu*(n_e)
    distributed on the bond.  sign=-1 is the physical t-J (t=1, J=2)."""
    P = graded_swap([0, 1, 1])
    ne = np.diag([0.0, 1, 1])
    h = sign * P - mu * onsite(ne, 3)
    return h, 3, np.array([[0, 0], [1, 1], [1, -1]])


def osp12(theta, mu=0.0):
    """OSp(1|2) invariant circle cos(theta) P_g + sin(theta) K_g on (boson, f1, f2), K_g the
    projector on the graded singlet |bb> + |f1 f2> - |f2 f1>."""
    P = graded_swap([0, 1, 1])
    s = np.zeros(9, complex); s[0 * 3 + 0] = 1; s[1 * 3 + 2] = 1; s[2 * 3 + 1] = -1
    K = 3 * proj(3, s)
    h = np.cos(theta) * P + np.sin(theta) * K - mu * onsite(np.diag([0.0, 1, 1]), 3)
    return h, 3, np.array([[0], [1], [-1]])


def blbq(theta):
    x = ss(1)
    return np.cos(theta) * x + np.sin(theta) * x @ x, 3, np.array([[1], [0], [-1]])


def spin1_xxz_d(delta, D):
    sx, sy, sz = spin(1)
    h = np.kron(sx, sx) + np.kron(sy, sy) + delta * np.kron(sz, sz) + D * onsite(sz @ sz, 3)
    return h, 3, np.array([[1], [0], [-1]])


# ------------------------------------------------------------ d = 4
def hubbard(U, mu, t=1.0):
    """Two fermion species per site, states (0, up, down, updown).  Graded hopping for each
    species, U n_up n_dn on site, -mu n."""
    # single species hopping as graded permutation on (0,1): P_g(1|1) = 1 - n_i - n_j + hop
    Pg = graded_swap([0, 1]).real
    hop1 = Pg - np.eye(4) + np.kron(np.diag([0, 1.0]), np.eye(2)) + np.kron(np.eye(2), np.diag([0, 1.0]))
    # species ordering: site state = (n_up, n_dn) -> index 2*n_up + n_dn ; bond index (a_i, a_j)
    I2 = np.eye(2)
    # up hopping acts on (up_i, up_j), identity on (dn_i, dn_j); fermionic sign from dn's between: handled by
    # graded structure: treat the four states with parity n_up+n_dn; build with full graded permutation algebra.
    # Simplest correct route: construct via Jordan-Wigner on the two-site, two-orbital problem (orbital order up_i, dn_i, up_j, dn_j).
    def cdag(n, k):  # JW creation operator k of n modes
        op = np.array([[1.0]])
        for m in range(n):
            if m < k:
                op = np.kron(op, np.diag([1.0, -1.0]))
            elif m == k:
                op = np.kron(op, np.array([[0, 0], [1.0, 0]]))  # |0> -> |1>, basis (0,1)
            else:
                op = np.kron(op, np.eye(2))
        return op
    n = 4
    c = [cdag(n, k).T for k in range(n)]
    cd = [cdag(n, k) for k in range(n)]
    num = [cd[k] @ c[k] for k in range(n)]
    h = -t * (cd[0] @ c[2] + cd[2] @ c[0] + cd[1] @ c[3] + cd[3] @ c[1])
    h += U / 2 * (num[0] @ num[1] + num[2] @ num[3])
    h -= mu / 2 * sum(num)
    # modes are (up_i, dn_i, up_j, dn_j) with basis index 8*n_upi + 4*n_dni + 2*n_upj + n_dnj = 4*(2 n_upi + n_dni) + (2 n_upj + n_dnj): site-major, good.
    return h, 4, np.array([[0, 0], [1, 0], [0, 1], [1, 1]])


def kugel_khomskii(x, y, J=1.0):
    """Spin-orbital h = J (S.S + x)(T.T + y), S and T spin-1/2, site state = spin x orbital."""
    s = ss(0.5)
    I4 = np.eye(4)
    # reorder: site index = 2*spin + orbital; bond basis (s_i, o_i, s_j, o_j) -> kron(spin bond, orbital bond) has order (s_i,s_j,o_i,o_j)
    SS = np.kron(s, np.eye(4)).reshape(2, 2, 2, 2, 2, 2, 2, 2)  # (si',sj',oi',oj', si,sj,oi,oj)
    TT = np.kron(np.eye(4), s).reshape(2, 2, 2, 2, 2, 2, 2, 2)
    perm = (0, 2, 1, 3, 4, 6, 5, 7)  # -> (si',oi',sj',oj', si,oi,sj,oj)
    SS = SS.transpose(perm).reshape(16, 16); TT = TT.transpose(perm).reshape(16, 16)
    h = J * (SS + x * np.eye(16)) @ (TT + y * np.eye(16))
    return h, 4, np.array([[1, 1], [1, -1], [-1, 1], [-1, -1]])


def two_chain_coupled(lam, mode='zz'):
    """Two spin-1/2 Heisenberg chains on one site: S.S' + T.T' + lam * coupling."""
    s = ss(0.5)
    SS = np.kron(s, np.eye(4)).reshape(2, 2, 2, 2, 2, 2, 2, 2)
    TT = np.kron(np.eye(4), s).reshape(2, 2, 2, 2, 2, 2, 2, 2)
    perm = (0, 2, 1, 3, 4, 6, 5, 7)
    SS = SS.transpose(perm).reshape(16, 16); TT = TT.transpose(perm).reshape(16, 16)
    sx, sy, sz = spin(0.5)
    # cross coupling S_i . T_{i+1} + T_i . S_{i+1}
    def op2(a, b, c, e):  # a on spin_i, b on orb_i, c on spin_j, e on orb_j
        return np.kron(np.kron(a, b), np.kron(c, e))
    I2 = np.eye(2)
    cross = sum(op2(a, I2, I2, a) + op2(I2, a, a, I2) for a in (sx, sy, sz))
    rung = onsite(sum(np.kron(a, a) for a in (sx, sy, sz)), 4)
    if mode == 'cross':
        h = SS + TT + lam * cross
        return h, 4, np.array([[2], [0], [0], [-2]])  # only the total S^z + T^z is conserved
    elif mode == 'rung':
        h = SS + TT + lam * rung
    else:
        h = SS + TT + lam * SS @ TT
    return h, 4, np.array([[1, 1], [1, -1], [-1, 1], [-1, -1]])


def graded_perm(par, mus, sign=1):
    """gl(N|M) chain: sign * P_g - sum_a mu_a n_a on the bond.  par: parities."""
    d = len(par)
    h = sign * graded_swap(par)
    for a, m in enumerate(mus):
        na = np.zeros((d, d)); na[a, a] = 1
        h = h - m * onsite(na, d)
    charges = np.eye(d, dtype=int)[:, 1:]  # independent species numbers (drop one)
    return h, d, charges
