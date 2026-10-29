"""Families and two-parameter cuts for the floating-Ising search, plus calibration controls.

Bond convention of contrib/d2_classification/ed.py (and vumps.bond):
H = sum_j [s_j^T J s_{j+1} + f . s_j], Pauli matrices, J = diag(Jx,Jy,Jz) + A(D), A_ab = eps_abc D_c,
so D_z multiplies X_j Y_{j+1} - Y_j X_{j+1}.

Families
  Z2 (C67 eq. (1), C68 anchor): (Jx,Jy,Jz,Dz,fz) = (1+g, 1-g, u, D, b).  The C68 chiral XXZ anchor is
      (g,u,D,b) = (0,0.4,0.5,0).  Strict single-Z2 stratum for g != 0 and D b != 0.
  J67 (C67 eq. (13)): (Jx,Jy,u,D,b) = (-1,-0.5,0.25,D,-gf), the interacting Ising/chiral-liquid junction.
  N68 (C68 Section 4, N family, with the field strength decoupled from eta):
      S = diag(1+eta,1-eta,0.4), D = (eta,2eta,0.5), f = b (3,4,5)/sqrt(50).
  P1282 (C68 Section 7): J = J0 + s e_z e_z^T, f = f0 + t e_z, point 1282 of d2_followup_220_1282.

Written for contrib/d2_floating_ising (Claude Opus 5.5 subagent, 2026-10-29).
"""
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
P1282 = json.load(open(os.path.join(HERE, '..', 'd2_followup_220_1282', 'candidates_params.json')))['1282']
NHAT = np.array([3., 4., 5.]) / np.sqrt(50.)


def J_from(Jd, D):
    J = np.diag(np.asarray(Jd, float))
    Dx, Dy, Dz = D
    J[0, 1] += Dz; J[1, 0] -= Dz
    J[1, 2] += Dx; J[2, 1] -= Dx
    J[2, 0] += Dy; J[0, 2] -= Dy
    return J


def family(fam, p):
    """Return (J, f) for family fam with parameter dict p."""
    if fam == 'Z2':
        g, u, D, b = p['g'], p['u'], p['D'], p['b']
        return J_from([1 + g, 1 - g, u], [0, 0, D]), np.array([0, 0, b], float)
    if fam == 'J67':
        return J_from([-1., -0.5, 0.25], [0, 0, p['D']]), np.array([0, 0, -p['gf']], float)
    if fam == 'N68':
        e = p['eta']
        return J_from([1 + e, 1 - e, 0.4], [e, 2 * e, 0.5]), p['b'] * NHAT
    if fam == 'P1282':
        J = np.array(P1282['J'], float); J[2, 2] += p['s']
        f = np.array(P1282['field'], float); f[2] += p['t']
        return J, f
    raise KeyError(fam)


# two-parameter cuts: name -> (family, fixed params, (x name, x0, x1), (y name, y0, y1))
CUTS = {
    'A': ('Z2', dict(u=0.4, D=0.5), ('g', 0.0, 0.8), ('b', 0.0, 2.4)),
    'B': ('Z2', dict(g=0.3, D=0.5), ('u', -0.8, 1.2), ('b', 0.0, 2.4)),
    'C': ('Z2', dict(u=0.4, b=1.2), ('g', 0.0, 0.8), ('D', 0.0, 0.8)),
    'J': ('J67', dict(), ('D', 0.0, 0.4), ('gf', 1.6, 2.4)),
    'N': ('N68', dict(), ('eta', 0.0, 0.4), ('b', 0.0, 2.4)),
    'P': ('P1282', dict(), ('s', -1.0, 1.0), ('t', -1.0, 1.0)),
}


def cut_points(cut, n=17):
    fam, fixed, (xn, x0, x1), (yn, y0, y1) = CUTS[cut]
    pts = []
    for i, x in enumerate(np.linspace(x0, x1, n)):
        for j, y in enumerate(np.linspace(y0, y1, n)):
            p = dict(fixed); p[xn] = float(x); p[yn] = float(y)
            pts.append((f'{cut}_{i:02d}_{j:02d}', fam, p))
    return pts


def parse_name(name):
    """'Z2:g=0.3,u=0.4,D=0.5,b=1.2' -> (fam, p)."""
    fam, rest = name.split(':')
    p = {k: float(v) for k, v in (kv.split('=') for kv in rest.split(','))}
    return fam, p


def control_term(name):
    """(h, d, c_exact) for the calibration chains, all nearest neighbour (two-site term)."""
    import vumps as V
    if name == 'ising':
        return V.bond(np.diag([-1., 0, 0]), [0, 0, -1.]), 2, 0.5
    if name == 'xx':
        return V.bond(np.diag([1., 1., 0.]), [0, 0, 0]), 2, 1.0
    if name == 'isxx':
        # Ising (x) XX on a d = 4 site (Ising velocity 2, XX velocity 4), c = 3/2
        hi = V.bond(np.diag([-1., 0, 0]), [0, 0, -1.]); hx = V.bond(np.diag([1., 1., 0.]), [0, 0, 0])
        one = np.eye(4).reshape(2, 2, 2, 2)
        T = (np.einsum('abcd,efgh->aebfcgdh', hi.reshape(2, 2, 2, 2), one)
             + np.einsum('abcd,efgh->aebfcgdh', one, hx.reshape(2, 2, 2, 2)))
        return T.reshape(16, 16), 4, 1.5
    raise KeyError(name)
