"""Projected n=1 model of the wrapped chain in cell coordinates.

Each cell is one site with a length bit l (0: length b=4, 1: length b+1=5)
and the payload qubit z; local index 2*l+z. The bond is

  h = sum_{l,l'} |l l'><l l'| (x) [K_(l+1)(l'+1) - c_11]  -  j_l (|10><01|+h.c.)_len (x) I,

with K_pq the rebuilt second-order source bond (units lambda^2/s) and the
first-order marker hopping j_l = u_mov beta_n s/lambda^2. Onsite
-kappa l - mu_p z sets the two densities. Everything is divided by the
payload hopping |t| = |<01|K_11|10>| so that energies are O(1).
"""
import numpy as np


def payload_hop(K):
    return K[1, 2]


def build(Kpq, jl=1.0, nu=0.35, mu_p=0.0, stress=1.0, drop_zz=False, decouple=False, V_override=None):
    antiherm = max(float(np.abs(v - v.conj().T).max()) for v in Kpq.values())
    Kpq = {k: (v + v.conj().T) / 2 for k, v in Kpq.items()}
    t_raw = payload_hop(Kpq['11'])
    # remove the uniform payload hopping phase by the exact gauge c_j -> e^{i j phi} c_j
    # (open or infinite chains); on every bond it acts as V^dag K V, V = 1 x diag(1, e^{i phi})
    phi = np.angle(t_raw) + np.pi
    Vg = np.kron(np.eye(2), np.diag([1.0, np.exp(1j * phi)]))
    Kpq = {k: Vg.conj().T @ v @ Vg for k, v in Kpq.items()}
    t = payload_hop(Kpq['11'])
    sc = abs(t)
    c = {k: np.trace(v).real / 4 for k, v in Kpq.items()}
    c11 = c['11']
    Vraw = c['22'] - c['21'] - c['12'] + c['11']
    # keep only the density interaction among the scalar parts; constant and
    # linear length terms are removed exactly (they are a chemical potential)
    lin = {'11': 0.0, '21': c['21'] - c11, '12': c['12'] - c11, '22': c['21'] - c11 + c['12'] - c11}
    elin_raw = (c['21'] - c11) + (c['12'] - c11)
    K = {k: (v - (c11 + lin[k]) * np.eye(4)) / sc for k, v in Kpq.items()}
    if decouple:
        K = {k: K['11'] for k in K}
    if stress != 1.0:
        base = K['11']
        new = {}
        for k, v in K.items():
            dlt = v - base
            sca = np.trace(dlt).real / 4
            new[k] = base + sca * np.eye(4) + stress * (dlt - sca * np.eye(4))
        K = new
    if drop_zz:
        K = {k: v - np.diag(np.diag(v)) + np.trace(v).real / 4 * np.eye(4) for k, v in K.items()}
    if V_override is not None:
        V = (np.trace(K['22'] - K['21'] - K['12'] + K['11']).real) / 4
        K['22'] = K['22'] + (V_override - V) * np.eye(4)
    h2 = np.zeros((16, 16), complex)
    for l1 in (0, 1):
        for l2 in (0, 1):
            Pl = np.zeros((4, 4))
            Pl[2 * l1 + l2, 2 * l1 + l2] = 1
            Kp = K['%d%d' % (l1 + 1, l2 + 1)]
            # basis (l1 z1)(l2 z2): reorder from (l1 l2)(z1 z2)
            h2 += np.kron(Pl, Kp).reshape(2, 2, 2, 2, 2, 2, 2, 2).transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(16, 16)
    hop = np.zeros((4, 4))
    hop[2, 1] = hop[1, 2] = -jl                        # |10><01| + h.c. on (l1 l2)
    h2 += np.kron(hop, np.eye(4)).reshape(2, 2, 2, 2, 2, 2, 2, 2).transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(16, 16)
    # free hard-core boson filling nu: onsite coefficient of n_l is -mu = 2 j_l cos(pi nu)
    elin = 0.0
    kappa = -2 * jl * np.cos(np.pi * nu)
    nl = np.diag([0., 0., 1., 1.])
    nz = np.diag([0., 1., 0., 1.])
    h1 = -kappa * nl - mu_p * nz
    V = (np.trace(K['22'] - K['21'] - K['12'] + K['11']).real) / 4
    info = dict(t=complex(t), t_raw=complex(t_raw), gauge_phi=float(phi), scale=float(sc), c11=float(c11), kappa=float(kappa), elin_removed=float(elin_raw / sc), V=float(V), jl=jl, antiherm_raw=antiherm)
    return h2, h1, info, dict(nl=nl, nz=nz)


def kpq_from_rebuild(rho, n=1):
    from precise import kbar_bonds
    return kbar_bonds(n, rho)
