"""Bonds used in contrib/d2_two_fluid.

A model is given by a spec string:
  'xx', 'xxz05', 'xxzm05', 'xxh'           calibrations at d = 2 (c = 1)
  'xx2', 'xx2s'                            two decoupled XX chains as one d = 4 site (c = 2);
                                           xx2 has unequal velocities and fields, xx2s identical chains
  'p:t,phi,V,mu,gL,|gR|,arg gR[,Delta,fL,fR]'   particle bond of twobody.bond (d = 2)
  'J:J11,...,J33,B1,B2,B3'                 Pauli bond of specs/_d2_background.md (d = 2)
Returns (h, d).
"""
import numpy as np
from twobody import bond, from_pauli, PAULI, I2

X, Y, Z = PAULI


def xx_bond(J=1.0, Jz=0.0, h=0.0):
    return J * (np.kron(X, X) + np.kron(Y, Y)) + Jz * np.kron(Z, Z) + h * 0.5 * (np.kron(Z, I2) + np.kron(I2, Z))


def two_chains(ha, hb):
    """d = 4 site = qubit a (x) qubit b; bond = ha on the a qubits + hb on the b qubits."""
    A = ha.reshape(2, 2, 2, 2)          # a1' a2' a1 a2
    B = hb.reshape(2, 2, 2, 2)
    Ia = np.eye(2)
    # site index (a, b); bond index (a1 b1, a2 b2)
    Ha = np.einsum('wxyz,mn,op->wmxoynzp', A, Ia, Ia)      # (a1',b1',a2',b2',a1,b1,a2,b2)
    Hb = np.einsum('mn,op,wxyz->mwoxnypz', Ia, Ia, B)
    return (Ha + Hb).reshape(16, 16)


def get(spec):
    if spec == 'xx':
        return xx_bond(), 2
    if spec == 'xxz05':
        return xx_bond(1.0, 0.5), 2
    if spec == 'xxzm05':
        return xx_bond(1.0, -0.5), 2
    if spec == 'xxh':
        return xx_bond(1.0, 0.0, 0.8), 2
    if spec == 'xx2':
        return two_chains(xx_bond(1.0, 0.0, 0.6), xx_bond(0.55, 0.0, -0.35)), 4
    if spec == 'xx2s':
        return two_chains(xx_bond(), xx_bond()), 4
    if spec.startswith('p:'):
        v = [float(x) for x in spec[2:].split(',')]
        t, phi, V, mu, gL, gRa, gRp = v[:7]
        Delta, fL, fR = (v[7:] + [0, 0, 0])[:3]
        return bond(t, phi, V, mu, gL, gRa * np.exp(1j * gRp), Delta, fL, fR), 2
    if spec.startswith('J:'):
        v = [float(x) for x in spec[2:].split(',')]
        return from_pauli(np.array(v[:9]).reshape(3, 3), v[9:12]), 2
    raise KeyError(spec)
