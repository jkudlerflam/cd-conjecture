"""Local terms for the D2 candidates and the calibration controls (Pauli convention)."""
import json, os
import numpy as np
from ed import bond, PAULI, I2

HERE = os.path.dirname(os.path.abspath(__file__))
X, Y, Z = PAULI


def Rz(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.]])


def dmz(Jd, D):
    """diag(Jd) + DM along z, J_xy = D = -J_yx."""
    J = np.diag(np.asarray(Jd, float))
    J[0, 1] += D; J[1, 0] -= D
    return J


# name: (J, field, expected c, note)
CONTROLS = {
    'ising':      (np.diag([-1., 0, 0]), [0, 0, -1.], 0.5, 'critical transverse-field Ising'),
    'xx':         (np.diag([1., 1, 0]), [0, 0, 0], 1.0, 'XX'),
    'heis':       (np.diag([1., 1, 1]), [0, 0, 0], 1.0, 'Heisenberg AF (log corrections)'),
    'xxz05':      (np.diag([1., 1, .5]), [0, 0, 0], 1.0, 'XXZ Delta=0.5'),
    'xxzm05':     (np.diag([1., 1, -.5]), [0, 0, 0], 1.0, 'XXZ Delta=-0.5'),
    'xxzm08':     (np.diag([1., 1, -.8]), [0, 0, 0], 1.0, 'XXZ Delta=-0.8'),
    'xxzm095':    (np.diag([1., 1, -.95]), [0, 0, 0], 1.0, 'XXZ Delta=-0.95 (near ferromagnet)'),
    'xxzm099':    (np.diag([1., 1, -.99]), [0, 0, 0], 1.0, 'XXZ Delta=-0.99 (ferromagnetic crossover)'),
    'xxzm0999':   (np.diag([1., 1, -.999]), [0, 0, 0], 1.0, 'XXZ Delta=-0.999 (ferromagnetic crossover)'),
    'xxz_dm':     (dmz([1., 1, .5], .6), [0, 0, 0], 1.0, 'XXZ 0.5 + DM_z 0.6 = twisted XXZ'),
    'xxzm07_dm_h': (dmz([-1., -1, .7], .8) , [0, 0, .3], 1.0, 'FM-XY + AF zz + DM_z + h_z: twisted XXZ in field'),
    'fm_twist':   (-Rz(0.6), [0, 0, 0], None, 'isotropic ferromagnet in a spiral frame (z=2)'),
}


def cluster_ising(lam=1.0):
    """h = -X Z X + lam Y Y 1 on three sites (c = 3/2 at lam = 1)."""
    return -np.kron(np.kron(X, Z), X) + lam * np.kron(np.kron(Y, Y), I2)


def local_term(name):
    pars = json.load(open(os.path.join(HERE, 'candidates_params.json')))
    if name in pars:
        p = pars[name]
        return bond(p['J'], p['field']), {'J': p['J'], 'field': p['field'], 'kind': 'candidate'}
    if name in CONTROLS:
        J, f, c, note = CONTROLS[name]
        return bond(J, f), {'J': np.asarray(J).tolist(), 'field': list(f), 'kind': 'control', 'c': c, 'note': note}
    if name.startswith('cluster_ising'):
        lam = float(name.split('_')[-1]) if name.count('_') == 2 else 1.0
        return cluster_ising(lam), {'kind': 'control', 'c': 1.5 if lam == 1.0 else None,
                                    'note': f'cluster-Ising, three-site, lambda={lam}'}
    extra = os.path.join(HERE, 'extra_points.json')
    if os.path.exists(extra):
        ex = json.load(open(extra))
        if name in ex:
            p = ex[name]
            return bond(p['J'], p['field']), {'J': p['J'], 'field': p['field'], 'kind': 'extra', 'note': p.get('note', '')}
    raise KeyError(name)
