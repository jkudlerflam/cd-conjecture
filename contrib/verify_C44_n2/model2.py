"""Projected n = 2 chain in cell coordinates and its exactly solvable references.

Each cell is one site with local index 4*l + z: l the length bit (0: short
length n+3, 1: long length n+4) and z = z0 + 2 z1 the two payload modes
(flavour f is bit f; Jordan-Wigner order inside a cell f = 0, 1, as carried by
the rebuilt bond). The bond is

  h = sum_{l l'} P_l P_l' (x) kappa_(l+1)(l'+1) + V n_l n_l' - j_l (|10><01| + h.c.)_len,

with kappa_pq the traceless part of the s/delta -> 0 bond K_pq (mpkbar.py),
the linear length terms of the scalar parts absorbed into the length
chemical potential, and V = c_22 - c_21 - c_12 + c_11 the length density
interaction. The uniform payload hopping phase is removed by the exact gauge
c_j -> e^{i j phi} c_j, chosen to centre the two Fermi pockets at k = 0 and pi.
Energies are in units of v_p, the mean payload Fermi velocity of the quadratic
part. j_l is set so that the free length velocity 2 j_l sin(pi nu) equals v_p.

Models
  M0  decoupled free reference: quadratic payload part + free length XX chain, c = 3 exactly
  M1  scheduled point (full projected bond)
  M2  M1 with every traceless length-payload term multiplied by STRESS
  C0  M0 with the length mode frozen (all cells short): c = 2 exactly
  C1  M1 with the length mode frozen: control, expect c = 2
"""
import os
import numpy as np
import mpmath as mp
from ferm import decompose, quad_operator
from scan import bands

HERE = os.path.dirname(os.path.abspath(__file__))
SCHED = dict(rho='0.008', g='0.22', w='1.8', a='3.5', eta='2.0', v='0.3', D='1.0', nu=0.375, stress=1000.0)


def load(rho, g, w, a, eta, v, D):
    z = np.load(os.path.join(HERE, 'kbar_n2_rho%s_g%s_w%s_a%s_eta%s_v%s_D%s.npz' % (rho, g, w, a, eta, v, D)))
    K = {k: z['K' + k] for k in ('11', '12', '21', '22')}
    return K, float(mp.mpf(str(z['V']))) * 1e12


def gauge(K, phi):
    Vg = np.kron(np.eye(4), np.diag([np.exp(1j * phi * bin(z).count('1')) for z in range(4)]))
    return {k: Vg.conj().T @ v @ Vg for k, v in K.items()}


def payload_data(K11):
    sc, h, rest = decompose(K11, 2)
    T = h[:2, 2:]
    E, fp = bands(T)
    return T, fp, rest


def embed(K):
    """sum_{l1 l2} |l1 l2><l1 l2| (x) K_(l1+1)(l2+1) in the basis (l1 z1)(l2 z2)."""
    h2 = np.zeros((64, 64), complex)
    for l1 in (0, 1):
        for l2 in (0, 1):
            Pl = np.zeros((4, 4))
            Pl[2 * l1 + l2, 2 * l1 + l2] = 1
            Kp = K['%d%d' % (l1 + 1, l2 + 1)]
            h2 += np.kron(Pl, Kp).reshape(2, 2, 4, 4, 2, 2, 4, 4).transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(64, 64)
    return h2


def len_ops():
    nl = np.kron(np.diag([0., 1.]), np.eye(4))
    hop = np.zeros((4, 4))
    hop[2, 1] = hop[1, 2] = 1.0
    # (l1 l2)(z1 z2) -> (l1 z1)(l2 z2)
    H = np.kron(hop, np.eye(16)).reshape(2, 2, 4, 4, 2, 2, 4, 4)
    H = H.transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(64, 64)
    nn = np.kron(nl, nl)
    return nl, H, nn


def build(name, rho=None, g=None, nu=None, stress=None):
    rho = rho or SCHED['rho']
    g = g or SCHED['g']
    w, a, eta, v, D = (SCHED[k] for k in ('w', 'a', 'eta', 'v', 'D'))
    nu = SCHED['nu'] if nu is None else nu
    stress = SCHED['stress'] if stress is None else stress
    K, V = load(rho, g, w, a, eta, v, D)
    T0, fp0, _ = payload_data(K['11'])
    # pair the zeros of each band and take the pocket centres
    byb = {}
    for b, k, v in fp0:
        byb.setdefault(b, []).append(k)
    cen = []
    for b, kk in byb.items():
        a, c = sorted(kk)
        cen.append((a + c) / 2 if c - a < np.pi else (a + c) / 2 + np.pi)
    # gauge so that one pocket is centred at 0 (the other then sits at pi)
    phi = min(cen, key=lambda x: abs(np.angle(np.exp(1j * x))))
    phi = float(np.angle(np.exp(1j * phi)))
    K = gauge(K, phi)
    T, fp, _ = payload_data(K['11'])
    vp = float(np.mean([abs(x[2]) for x in fp]))
    K = {k: v / vp for k, v in K.items()}
    V = V / vp
    T = T / vp
    jl = 1.0 / (2 * np.sin(np.pi * nu))
    kap = -2 * jl * np.cos(np.pi * nu)
    nl, Hhop, nn = len_ops()
    nz = np.kron(np.eye(2), np.diag([0., 1., 1., 2.]))
    if name in ('M0', 'C0'):
        q = quad_operator(np.block([[np.zeros((2, 2)), T], [T.conj().T, np.zeros((2, 2))]]), 2)
        Kuse = {k: q for k in K}
        Vuse = 0.0
    else:
        Kuse = dict(K)
        Vuse = V
        if name == 'M2':
            Kuse = {k: K['11'] + stress * (v - K['11']) for k, v in K.items()}
    h2 = embed(Kuse) - jl * Hhop + Vuse * nn
    h1 = -kap * nl
    if name in ('C0', 'C1'):
        h1 = 5.0 * nl                     # length band empty: every cell short
    mixed = max(float(np.abs(Kuse[k] - Kuse['11']).max()) for k in Kuse)
    info = dict(name=name, rho=float(rho), g=float(g), w=float(w), a=float(a), eta=float(eta), v_bus=float(v), D=float(D), nu=nu, jl=jl, kappa=kap, V=V, vp_raw=vp * 1e-12,
                gauge_phi=phi, fermi=[(int(b), k, v / vp) for b, k, v in fp], mixed_max=mixed,
                stress=stress if name == 'M2' else 1.0, T=[[T[i, j].real, T[i, j].imag] for i in range(2) for j in range(2)])
    return h2, h1, info, dict(nl=nl, nz=nz), T


def free_profile_payload(L, T):
    """Exact open-chain entanglement profile of the quadratic payload ground state."""
    H = np.zeros((2 * L, 2 * L), complex)
    for i in range(L - 1):
        H[2 * i:2 * i + 2, 2 * i + 2:2 * i + 4] = T
        H[2 * i + 2:2 * i + 4, 2 * i:2 * i + 2] = T.conj().T
    e, v = np.linalg.eigh(H)
    occ = v[:, e < 0]
    C = occ @ occ.conj().T
    prof = []
    for l in range(1, L):
        ev = np.linalg.eigvalsh(C[:2 * l, :2 * l])
        ev = ev[(ev > 1e-15) & (ev < 1 - 1e-15)]
        prof.append(float(-(ev * np.log(ev) + (1 - ev) * np.log(1 - ev)).sum()))
    return np.array(prof), float((e < 0).sum() / L)


def free_profile_xx(L, hop, onsite):
    H = np.zeros((L, L))
    for i in range(L - 1):
        H[i, i + 1] = H[i + 1, i] = hop
    H += onsite * np.eye(L)
    e, v = np.linalg.eigh(H)
    occ = v[:, e < 0]
    C = occ @ occ.T
    prof = []
    for l in range(1, L):
        ev = np.linalg.eigvalsh(C[:l, :l])
        ev = ev[(ev > 1e-15) & (ev < 1 - 1e-15)]
        prof.append(float(-(ev * np.log(ev) + (1 - ev) * np.log(1 - ev)).sum()))
    return np.array(prof), float((e < 0).sum() / L)


def exact_profile(name, L, info, T):
    Sp, _ = free_profile_payload(L, T)
    if name == 'C0':
        return Sp
    Sl, _ = free_profile_xx(L, -info['jl'], -info['kappa'])
    return Sp + Sl
