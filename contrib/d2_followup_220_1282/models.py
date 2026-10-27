"""Local terms for the two candidates and the calibration controls.

The candidates use the bond convention of contrib/d2_classification (ed.bond):
h = sum_ab J_ab s^a (x) s^b + sum_a f_a (s^a (x) 1 + 1 (x) s^a)/2, Pauli matrices,
so H = sum_j [s_j^T J s_{j+1} + f . s_j] on the infinite chain.
"""
import json, os
import numpy as np
from vumps import bond, PAULI, I2

HERE = os.path.dirname(os.path.abspath(__file__))
X, Y, Z = PAULI
PARS = os.path.join(HERE, 'candidates_params.json')


def dmz(Jd, D):
    J = np.diag(np.asarray(Jd, float)); J[0, 1] += D; J[1, 0] -= D
    return J


# name: (local term, range r, exact c, note)
def model(name):
    if name in ('220', '1282'):
        p = json.load(open(PARS))[name]
        return bond(p['J'], p['field']), 2, None, 'D2 candidate ' + name
    if name == 'ising':
        return bond(np.diag([-1., 0, 0]), [0, 0, -1.]), 2, 0.5, 'critical TFIM'
    if name == 'xx':
        return bond(np.diag([1., 1, 0]), [0, 0, 0]), 2, 1.0, 'XX'
    if name == 'xxz05':
        return bond(np.diag([1., 1, .5]), [0, 0, 0]), 2, 1.0, 'XXZ Delta=0.5'
    if name == 'xxzm05':
        return bond(np.diag([1., 1, -.5]), [0, 0, 0]), 2, 1.0, 'XXZ Delta=-0.5'
    if name == 'xxzm08':
        return bond(np.diag([1., 1, -.8]), [0, 0, 0]), 2, 1.0, 'XXZ Delta=-0.8'
    if name == 'xydm':
        # XY + DM_z + h_z: free fermions, no U(1), incommensurate Fermi points, c = 1
        return bond(dmz([1., .6, 0], .5), [0, 0, .6]), 2, 1.0, 'XY(1,0.6)+D_z 0.5+h_z 0.6 (free fermion, c=1)'
    if name == 'xxzm07_dm_h':
        return bond(dmz([-1., -1, .7], .8), [0, 0, .3]), 2, 1.0, 'twisted XXZ in a field (U(1), c=1)'
    if name == 'cluster_ising':
        # three-site, c = 3/2
        return -np.kron(np.kron(X, Z), X) + np.kron(np.kron(Y, Y), I2), 3, 1.5, 'cluster-Ising, c=3/2'
    if name == 'xx_nnn':
        # two decoupled XX chains (next-nearest neighbour), c = 2
        return np.kron(np.kron(X, I2), X) + np.kron(np.kron(Y, I2), Y), 3, 2.0, 'two decoupled XX chains, c=2'
    raise KeyError(name)
