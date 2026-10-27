"""Port-idle Feynman-Kitaev cell versus its target gadget (Claude's independent check).

Registers of one cell, in this order: clock t = 0..T, port p (qubit), partner s (qubit),
logical l (qubit), ancilla c (qubit).  Gates act on (l, c) only, so the port is never
touched by the computation.  Penalties h_t act on (p, s, c) at clock time t.

Target gadget: H_A = sum_t W_t^dag h_t W_t restricted to c = 0, acting on (p, s, l).
Its kernel is the code G|i> = |Omega>_{ps} |i>_l (dimension 2).

Checks, as functions of the clock coupling J_c:
  (1) the simulator has an exactly two-fold zero-energy ground space (frustration-free part),
  (2) port compression G~^dag F G~ (exactly zero for the port-idle cell),
  (3) the matrix-valued boundary resolvent R_ab(z) = G^dag F_a (H + z)^-1 F_b G,
  (4) the two-cell second-order exchange K2 (spectral formula, all excited states).
The same quantities are computed for a variant in which ONE gate acts on (p, c).

Run: python3 port_idle_check.py   (numpy, scipy; a few seconds)
"""
import numpy as np
from functools import reduce

rng = np.random.default_rng(7)
I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1.0, -1.0]).astype(complex)
PAULI = [X, Y, Z]


def kron(*ops):
    return reduce(np.kron, ops)


def haar(n):
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    q, r = np.linalg.qr(a)
    return q * (np.diag(r) / abs(np.diag(r)))


def rand_herm(n, scale=1.0):
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return scale * (a + a.conj().T) / 2


# input space ordering: p, s, l, c  (dim 16)
omega = np.zeros(4, complex)
omega[0] = omega[3] = 1 / np.sqrt(2)           # |Omega>_{ps}
P_om = np.outer(omega, omega.conj())
Q_om = np.eye(4) - P_om                        # syndrome projector on the pair
syn = np.linalg.eigh(Q_om)[1][:, 1:]           # 4x3 orthonormal syndrome frame


def build(T, touch_port=False, g=0.05):
    """Return gates U_t (16x16 on input space) and penalties h_t (16x16), t = 0..len(U).

    Port-idle: every gate acts on (l, c).  Touched variant: at the middle of the schedule a
    random gate V on (p, c) is applied and immediately undone (V then V^dag), with zero
    penalty at those two clock times, so the code is still exactly annihilated but the
    port is acted on during one clock time (C9-style dilution).
    """
    U, tags = [], []
    for t in range(1, T + 1):
        U.append(kron(I2, I2, haar(4))); tags.append('lc')
        if touch_port and t == T // 2:
            u = haar(4).reshape(2, 2, 2, 2)                 # (p',c',p,c)
            V = np.einsum('acbd,ef,gh->aegcbfhd', u, I2, I2).reshape(16, 16)
            U.append(V); tags.append('touch')
            U.append(V.conj().T); tags.append('undo')
    W = [np.eye(16, dtype=complex)]
    for u in U:
        W.append(u @ W[-1])
    h = [kron(Q_om, I2, I2)]                                 # t = 0: enforce |Omega>_{ps}
    for tag in tags:
        if tag in ('touch', 'undo'):
            h.append(np.zeros((16, 16), complex))
        else:
            K = syn @ rand_herm(3) @ syn.conj().T            # K|Omega> = 0
            h.append(g * kron(K, I2, rand_herm(2)))          # (p,s) x l x c
    return U, W, h


def target(W, h):
    """H_A on inputs with c = 0: indices (p, s, l)."""
    Heff = sum(Wt.conj().T @ ht @ Wt for Wt, ht in zip(W, h))
    sel = [i for i in range(16) if i % 2 == 0]               # c = 0
    return Heff[np.ix_(sel, sel)]


def simulator(U, h, Jc):
    T = len(U)
    D = 16 * (T + 1)
    H = np.zeros((D, D), complex)
    blk = lambda t: slice(16 * t, 16 * (t + 1))
    for t in range(1, T + 1):
        H[blk(t), blk(t)] += Jc / 2 * np.eye(16)
        H[blk(t - 1), blk(t - 1)] += Jc / 2 * np.eye(16)
        H[blk(t), blk(t - 1)] += -Jc / 2 * U[t - 1]
        H[blk(t - 1), blk(t)] += -Jc / 2 * U[t - 1].conj().T
    H[blk(0), blk(0)] += Jc * kron(I2, I2, I2, np.diag([0.0, 1.0]))   # ancilla init
    for t in range(T + 1):
        H[blk(t), blk(t)] += (T + 1) * h[t]
    return H


def code_target():
    G = np.zeros((8, 2), complex)        # (p,s,l) ordering, l last
    for i in range(2):
        e = np.zeros(2); e[i] = 1
        G[:, i] = np.kron(omega, e)
    return G


def code_sim(W, T):
    G = np.zeros((16 * (T + 1), 2), complex)
    for i in range(2):
        psi = kron(omega, np.eye(2)[i], np.array([1.0, 0.0]))
        G[:, i] = np.concatenate([W[t] @ psi for t in range(T + 1)]) / np.sqrt(T + 1)
    return G


def port_ops(T=None):
    if T is None:
        return [kron(P, I2, I2) for P in PAULI]                       # on (p,s,l)
    return [np.kron(np.eye(T + 1), kron(P, I2, I2, I2)) for P in PAULI]


def resolvent(H, G, F, z):
    R = np.linalg.solve(H + z * np.eye(len(H)), np.hstack([f @ G for f in F]))
    return np.array([[ (F[a] @ G).conj().T @ R[:, 2 * b:2 * b + 2] for b in range(3)]
                     for a in range(3)])                               # 3x3 blocks of 2x2


def K2(H, G, F, tol=1e-9):
    E, V = np.linalg.eigh(H)
    ex = E > tol
    E, V = E[ex], V[:, ex]
    A = np.array([V.conj().T @ (f @ G) for f in F])                    # (a, e, i)
    # M[a,b,e] (2x2) = A[a,e,:]^* A[b,e,:]
    Mae = np.einsum('aei,bej->abeij', A.conj(), A)
    C = 1.0 / (E[:, None] + E[None, :])
    K = -np.einsum('ef,abeij,abfkl->ikjl', C, Mae, Mae)
    return K.reshape(4, 4)


def run(touch, T=6, Jcs=(1e2, 1e3, 1e4, 1e5, 1e6)):
    U, W, h = build(T, touch_port=touch)
    T = len(U)
    HA = target(W, h)
    GA = code_target()
    FA = port_ops()
    eA = np.linalg.eigvalsh(HA)
    K2A = K2(HA, GA, FA)
    zs = [0.5, 1.0, 2.0]
    RA = [resolvent(HA, GA, FA, z) for z in zs]
    out = []
    for Jc in Jcs:
        H = simulator(U, h, Jc)
        E = np.linalg.eigvalsh(H)
        Gs = code_sim(W, T)
        Fs = port_ops(T)
        kern = np.linalg.norm(H @ Gs)
        comp = max(np.linalg.norm(Gs.conj().T @ f @ Gs) for f in Fs)
        Rerr = max(np.max(abs(resolvent(H, Gs, Fs, z) - RA[k])) for k, z in enumerate(zs))
        K2s = K2(H, Gs, Fs)
        out.append((Jc, E[:3], kern, comp, Rerr, np.linalg.norm(K2s - K2A)))
    return eA, K2A, out


if __name__ == '__main__':
    for touch in (False, True):
        eA, K2A, out = run(touch)
        print('\n=== port touched by one gate' if touch else '\n=== port-idle cell')
        print('target spectrum (lowest 4):', np.round(eA[:4], 6), '  ||K2_target|| =', round(float(np.linalg.norm(K2A)), 4))
        print(' J_c   E0,E1,E2 of simulator           ||H G~||   max|G~^dag F G~|  max|dR(z)|   ||dK2||')
        for Jc, E3, kern, comp, Rerr, dK in out:
            print(f'{Jc:6g}  {np.round(E3, 6)}  {kern:9.2e}  {comp:12.2e}  {Rerr:10.3e}  {dK:10.3e}')
