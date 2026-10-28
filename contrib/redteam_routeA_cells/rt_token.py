"""Leading order in the conversion amplitude g: the hopping matrix T from the
data-free single-token problem (positions x species), derived independently.

Gauge: neutral hops carry Z on the departure qubit; with W(d_j)=W(x_j)=Z_0..Z_{j-1}
and W(M)=P every neutral hop becomes scalar, conversions become JW fermions
(0->p: g c_j, 0->m: g c_j^dag), and the right-port neutral source carries P.
Then, with H2 = -P V (H0-E0)^{-1} V P,

  T_ij = -alpha^4 g^2 e^{i phi} [Amp1(i,j) - Amp2(i,j)]

  Amp1 = <M d0| R(m,p) [C_B R(m,0) C_A + C_A R(0,p) C_B] R(0,0) |M d0>
  Amp2 = <M d0| R(0,0) [C_B R(0,m) C_A + C_A R(p,0) C_B] R(p,m) |M d0>

R(s,s') = (K_s x 1 + 1 x K_s')^{-1} on (positions of A) x (positions of B),
C_A = |d_i><d_i| on A, C_B = |d_j><d_j| on B.  Channel 2 carries the minus
sign because P_A is applied after c_A^dag.
"""
import numpy as np
from rt_bond import P


def K_species(n, p, s):
    """Token Hamiltonian of one species on positions d_0..d_{n-1}, x_0..x_{n-1}, M,
    measured from the vacuum-walker ground energy."""
    m = 2 * n + 1
    K = np.zeros((m, m), complex)
    q = {'0': 0, 'p': 1, 'm': -1}[s]
    E0 = -2 * p['t'] * np.cos(np.pi / (n + 2))
    for j in range(n):
        K[j, j] = p['Delta'] + p['delta'] * q - E0
        K[n + j, n + j] = p['Delta'] + p['delta'] * q + p['Dx'] - E0
    K[2 * n, 2 * n] = p['Delta'] + p['delta'] * q + p['DM'] - E0

    def hop(to, fr, c):
        K[to, fr] += c
        K[fr, to] += np.conj(c)
    for j in range(n - 1):
        hop(j + 1, j, -p['t'] * np.exp(1j * p['ph_stub'][s]))
        hop(n + j + 1, n + j, -p['tx'] * np.exp(1j * p['ph_x'][s]))
    hop(2 * n, 2 * n - 1, -p['tx'] * np.exp(1j * p['ph_x'][s]))
    if p['td'] != 0:
        hop(2 * n, n - 1, -p['td'] * np.exp(1j * p['ph_d'][s]))
    hop(n, 0, -p['tsw'] * np.exp(1j * p['ph_sw'][s]))
    return K


def T_token(n, p, stub_only=False):
    """O(g^2) hopping matrix.  stub_only=True gives the Gram-limit reference in
    which both ports sit at d_0 of a bare stub (Gamma_ij = int f_i f_j)."""
    Ks = {s: K_species(n, p, s) for s in ('0', 'p', 'm')}
    m = 2 * n + 1
    if stub_only:
        Ks = {s: Ks[s][:n, :n] for s in Ks}
        m = n
    I = np.eye(m)
    R = {}
    for s in ('0', 'p', 'm'):
        for s2 in ('0', 'p', 'm'):
            R[(s, s2)] = np.linalg.inv(np.kron(Ks[s], I) + np.kron(I, Ks[s2]))
    srcA = np.zeros(m); srcA[0 if stub_only else 2 * n] = 1.0
    srcB = np.zeros(m); srcB[0] = 1.0
    v = np.kron(srcA, srcB)
    alpha = np.sqrt(2 / (n + 2)) * np.sin(np.pi / (n + 2))
    T = np.zeros((n, n), complex)
    for i in range(n):
        CA = np.kron(np.diag(np.eye(m)[i]), I)
        for j in range(n):
            CB = np.kron(I, np.diag(np.eye(m)[j]))
            a1 = v @ R[('m', 'p')] @ (CB @ R[('m', '0')] @ CA + CA @ R[('0', 'p')] @ CB) @ R[('0', '0')] @ v
            a2 = v @ R[('0', '0')] @ (CB @ R[('0', 'm')] @ CA + CA @ R[('p', '0')] @ CB) @ R[('p', 'm')] @ v
            if stub_only:
                T[i, j] = a1      # all species equal on the bare stub at delta = 0
            else:
                T[i, j] = -alpha ** 4 * p['g'] ** 2 * np.exp(1j * p['phi']) * (a1 - a2)
    return T


def rs_eigs(T):
    sc = np.abs(T).max()
    R = T * np.exp(-1j * np.angle(T[0, 0]))
    return np.linalg.eigvalsh((R + R.conj().T) / 2) / sc


def crawford(T, nth=7201):
    sc = np.abs(T).max()
    best, arg = -np.inf, None
    for th in np.linspace(0, 2 * np.pi, nth, endpoint=False):
        Rt = T * np.exp(-1j * th)
        e = np.linalg.eigvalsh((Rt + Rt.conj().T) / 2)[0] / sc
        if e > best:
            best, arg = e, th
    return best, arg
